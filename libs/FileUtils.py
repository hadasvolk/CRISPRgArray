import sys
import os
import shutil
from Bio import SeqIO
import re
import boto3
import botocore
from subprocess import Popen, PIPE
import logging
import datetime

MainTempDir = "/bigdisk/tmp"
s3_resource = boto3.resource('s3')
s3_client = boto3.client('s3')

class BasicFileUtils:

    def precopyDir(src, tmpdest):
        tmpdir = MainTempDir
        if len(tmpdest) > 0:
            tmpdir = tmpdest

        if not src.startswith('s3:'):
            if not os.path.isdir(src):
                print ("Input dir {} does not exist. Exiting ...".format(src))
                logging.error("Input dir {} does not exist".format(src))
                sys.exit(2)
            shutil.copytree(src, tmpdest)
        else:
            BasicFileUtils.processCommand("aws s3 sync " + src + "/ " + tmpdest + " --quiet")

        return tmpdir


    def fileExist(src):
        s = re.split('\/',src)
        bucket=s[2]
        key="/".join(s[3:])
        try:
            s3_resource.Object(bucket, key).load()
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                logging.error("{} does not exist".format(key))
                sys.exit(2)
            else:
                print("Connection to s3 failed")
                raise


    def precopyFile(src, tmpdest, dest):
        tmpdir = MainTempDir
        if len(tmpdest) > 0:
            tmpdir = tmpdest
        head, tail = os.path.split(src)
        otmpfile = tmpdir + "/" + tail
        if not src.startswith('s3:'):
            if not os.path.isfile(src):
                print ("{} file {} does not exist. Exiting ...".format(dest, src))
                logging.error("{} file {} does not exist".format(dest, src))
                sys.exit(2)
            shutil.copyfile(src, otmpfile)

        else:
            s = re.split('\/',src)
            bucket=s[2]
            key="/".join(s[3:])
            try:
                s3_resource.Bucket(bucket).download_file(key, otmpfile)
            except botocore.exceptions.ClientError as e:
                if e.response['Error']['Code'] == "404":
                    if dest == 'input':
                        print("The object does not exist.")
                        logging.error("The object does not exist")
                        sys.exit(2)
                else:
                    raise

        return otmpfile

    def postcopy(srcdir, destdir):
        if not destdir.startswith("s3:"):
            if not os.path.isdir(destdir):
                os.mkdir(destdir)
            shutil.copytree(srcdir, destdir, symlinks=False, ignore=None, dirs_exist_ok=True)
        else:
            s = re.split('\/',destdir)
            bucket=s[2]
            dstkey="/".join(s[3:])
            if dstkey.endswith('/'):
                dstkey = dstkey[:-1]
            BasicFileUtils.processCommand("aws s3 sync {} {} --quiet".format(srcdir, destdir))
            logging.shutdown()


    def processCommand(command):
        try:
            p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
            p.wait()
            output, error = p.communicate()
            if p.returncode != 0:
                print("{} -> Failed".format(command, output, error))
                logging.error("{} \nfailed {} -> {}".format(command, output, error))
                sys.exit(2)
        except Exception as e:
            print("{} -> Failed".format(command))
            logging.error("{} \nFailed {}".format(command, e))
            sys.exit(2)
        # logging.info("Successfully: {}".format(command))



    def uploadDirectory(srcpath, bucketname, dstkey):
        for root,dirs,files in os.walk(srcpath):
            for d0 in dirs:
                BasicFileUtils.uploadDirectory(root+'/'+d0, bucketname, dstkey+'/'+d0)
            for f0 in files:
                srcname = "{}/{}".format(root, f0)
                dstname = "{}/{}".format(dstkey, f0)
                #print(root, f0, dstkey, f0)
                print("Copying {} to {}".format(srcname, dstname))
                s3_client.upload_file(srcname, bucketname, dstname)


    def clean(clist):
        for f in clist:
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)

    # Clearing FASTA files from temp directory
    def clearFa(src):
        for f in os.listdir(src):
            if f.endswith('.fa'):
                BasicFileUtils.processCommand("rm {}/{}".format(src, f))

    # Initalize log file
    def create_log(path, curDate, var):
        log_file = "{}/gRNAs_{}_step1_{}.log".format(path, var, curDate)
        logging.basicConfig(filename=log_file,
                            level=logging.INFO,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.info("gRNAs Production for {} Started!".format(var))


class ExcelUtils:

    def read_sheet1(xfile):

        df = pd.read_excel(xfile,sheet_name='Sheet1', index_col=None)
        Name=df['User'][0]
        Email=df['Email'][0]
        Service=df['Service'][0]
        return Name, Email, Service

    def to_csv(xfile, sheetname):
        df = pd.read_excel(xfile,sheet_name=sheetname, index_col=None)
#        xpath = os.path.dirname(os.path.abspath(xfile))
        xfilename = os.path.splitext(xfile)[0]
        csvfile = xfilename + '.csv'
        df.to_csv(csvfile, header=True, index=None)
        return csvfile

    def copy_sheet(xsrcfile, srcsheetname, xdstfile, dstsheetname):
        df = pd.read_excel(xsrcfile,sheet_name=srcsheetname, index_col=None)
        # This section is sample code that creates a worbook in the current directory with 3 worksheets
        #df = pd.DataFrame(np.random.randn(10, 3), columns=list('ABC'))
        writer = pd.ExcelWriter(xdstfile, engine='xlsxwriter')
        df.to_excel(writer, sheet_name=dstsheetname, index=False)
        writer.close()


class FastaUtils:

    def batch_iterator(iterator, batch_size):
        """Returns lists of length batch_size.

        This can be used on any iterator, for example to batch up
        SeqRecord objects from Bio.SeqIO.parse(...), or to batch
        Alignment objects from Bio.AlignIO.parse(...), or simply
        lines from a file handle.

        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.
        """
        entry = True  # Make sure we loop once
        while entry:
            batch = []
            while len(batch) < batch_size:
                try:
                    entry = next(iterator)
                except StopIteration:
                    entry = None
                if entry is None:
                    # End of file
                    break
                batch.append(entry)
            if batch:
                yield batch

    def splitFastq (infile, outdir, batch_size):
        record_iter = SeqIO.parse(open(infile),"fastq")
        for i, batch in enumerate(FastaUtils.batch_iterator(record_iter, batch_size)):
            filename = outdir + "/" + batch[0].id + ".fastq"
            #filename = outdir + "/chunk_%i.fastq" % (i + 1)
            with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fastq")
                print("Wrote %i records to %s" % (count, filename))

    def splitFasta (infile, outdir, batch_size):
        record_iter = SeqIO.parse(open(infile),"fasta")
        for i, batch in enumerate(FastaUtils.batch_iterator(record_iter, batch_size)):
            filename = outdir + "/" + batch[0].id + ".fasta"
            #filename = outdir + "/chunk_%i.fasta" % (i + 1)
            with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
                print("Wrote %i records to %s" % (count, filename))


    def validProteinSequence(infastafile):
        records = list(SeqIO.parse(infastafile, "fasta"))
        valid = 'ACTG'
        for letter in records[0].seq:
            if letter not in valid:
                return True
        return False


    def getFastaSeq(filename):
        ifile = open(filename, 'rU')
        seqs = []
        for record in SeqIO.parse(ifile, "fasta"):
            sequence = str(record.seq).upper()
            seqs.append(sequence)

        return seqs

    def multi2linefasta(indir, outdir, filelist):
#
# May accept single file in filelist
#
        for items in filelist:
            mfasta = outdir +"/"+re.sub('\..*','',items)+'_twoline.fasta'
            ifile = open(indir+'/'+items,'rU')
            with open(mfasta, 'w') as ofile:
                for record in SeqIO.parse(ifile, "fasta"):
                    sequence = str(record.seq)
                    ofile.write('>' + record.id + '\n' + sequence + '\n')
