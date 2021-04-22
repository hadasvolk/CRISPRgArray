import sys
import os
import glob
import os.path
from os import path
import pandas as pd
import numpy as np
import string
import re
import shutil
import logging
import mysql.connector as mariadb
import boto3
from botocore.exceptions import ClientError
import datetime

from FileUtils import BasicFileUtils
from FileUtils import FastaUtils
from SequenceUtils import NucleotideUtils
import cfg

curDate = datetime.datetime.now().strftime("%d-%m-%Y")

class gRNA_Utils:
    def __init__ (self, ifile, file_ext, odir, tmpdir):
        # file_ext: \.xlsx , \.fa, \.fasta
        tmpfile = BasicFileUtils.precopyFile(ifile, tmpdir, 'input')
        logging.info("Input file: {}".format(ifile))
        logging.info("Output directory: {}".format(odir))
        logging.info("Temp directory: {}".format(tmpdir))

        files=os.listdir(tmpdir)
        self.ifile = ["{}/{}".format(tmpdir, f) for f in files if re.search(file_ext, f)][0]
        self.wd = tmpdir
        self.outdir = odir
        if odir.endswith('/'):
            self.outdir = odir[:-1]
        self.transtab = str.maketrans('TAGCtagc', 'ATCGATCG')

        self.__checkDBs()
        self.__readInfoSheet()
        #Get input fielname
        filename = ifile.split('/')[-1]
        self.filename = filename.replace(file_ext, '')

        ofile = "".join((self.outdir, '/', self.filename, '_OUT_{}.xlsx'.format(curDate)))
        self.out = BasicFileUtils.precopyFile(ofile, self.wd, 'output')


    def __del__ (self):
        BasicFileUtils.postcopy(self.wd, self.outdir)
        os.remove(cfg.Hg38_Genome)
        # shutil.rmtree(self.wd)


    def __checkDBs(self):
        if not path.exists(cfg.Hg38_Genome):
            cmd = "aws s3 sync s3://emendobio-compute/ReferenceGenomes/UCSC {}".format(os.path.abspath(cfg.Hg38_Genome))
            BasicFileUtils.processCommand(cmd)

    def __readInfoSheet(self):
        # Valdiate sheet1 from input file and extract information
        self.sheet1_header = ["user", "email", "purpose", "service"]
        try:
            self.sheet1 = pd.read_excel(self.ifile, sheet_name='Sheet1')
            assert len(self.sheet1.index) != 0
        except AssertionError:
            print("\nPlease insert user detalis in sheet1\n")
            logging.error("Empty sheet1")
            sys.exit(2)
        except:
            print("{} unable to read sheet1".format(self.ifile))
            logging.error("Unable to read sheet1")
            sys.exit(2)
        self.sheet1 = self.sheet1.iloc[:,0:4]
        self.sheet1.columns = [col.strip() for col in self.sheet1.columns]
        self.sheet1.columns = [s.lower() for s in self.sheet1.columns]
        if not all([s in self.sheet1_header for s in self.sheet1.columns]):
            print("{}, expected header for sheet1: \n{}".format(self.ifile,
                                                            self.sheet1_header))
            logging.error("Expected header for sheet1: \n{}".format(self.sheet1_header))
            sys.exit(2)
        self.user = self.sheet1['user'][0]
        self.mail = self.sheet1['email'][0]
        logging.info("Successfully read sheet1")

        # Valdiate Info sheet from input file and extract information
        expected = ["GeneName", "Chromosome", "RefSeqTranscript",
            "GenomeAssembly", "GeneStartPos", "GeneEndPos", "Strand", "SpacerSizes",
            "PAMsize", "TestingGuides", "AdditionalPAMs"]
        self.info_header = [s.lower() for s in expected]
        try:
            self.info = pd.read_excel(self.ifile, sheet_name='Info')
            assert len(self.info.index) != 0
        except AssertionError:
            print("\nPlease insert info\n")
            logging.error("Empty Info sheet")
            sys.exit(2)
        except:
            print("{} unable to read Info sheet".format(self.ifile))
            logging.error("Unable to read Info sheet")
            sys.exit(2)

        self.info = self.info.iloc[:,0:11]
        self.info.columns = [col.strip() for col in self.info.columns]
        self.info.columns = [s.lower() for s in self.info.columns]
        if not all([s in self.info_header for s in self.info.columns]):
            print("{}, expected header for Info: \n{}".format(self.ifile,
                                                            self.info_header))
            logging.error("Expected header for Info: \n{}".format(self.info_header))
            sys.exit(2)
        self.GeneName = self.info['genename'][0]
        try:
            self.Chromosome = int(self.info['chromosome'][0])
        except:
            self.Chromosome = self.info['chromosome'][0].split('r')[-1]
        self.RefSeqTranscript = self.info['refseqtranscript'][0].strip()
        print("RefSeqTranscript: " + self.RefSeqTranscript)
        self.GenomeAssembly = self.info['genomeassembly'][0]
        if self.GenomeAssembly != "hg38":
            print("Software supports only hg38 build")
            logging.error("Software supports only hg38 build")
            sys.exit(2)
        try:
            self.GeneStartPos = int(self.info['genestartpos'][0])
            self.GeneEndPos = int(self.info['geneendpos'][0])
        except:
            print("\nUnable to read gene positions from Info sheet\n")
            logging.error("Unable to read gene positions from Info sheet")
            sys.exit(2)
        self.strand = self.info['strand'][0][1]
        if self.strand != "-" and self.strand != "+":
            print("undefined strand. Please use: \"-\" or \"+\"")
            logging.error("undefined strand. Please use: \"-\" or \"+\"")
            sys.exit(2)

        num_format = re.compile("^[\-]?[1-9][0-9]*\.?[0-9]+$")
        self.SpacerSizes = [int(x) for x in self.info['spacersizes'][0].split(',')]
        self.PAMsize = int(self.info['pamsize'][0])
        self.PAMs = [x for x in self.info['additionalpams'].tolist() if str(x) != 'nan']
        self.TestGuides = [x.upper() for x in self.info['testingguides'].tolist() if (not re.match(num_format,str(x)) and str(x) != 'nan')]
        logging.info("Successfully read Info sheet")
        logging.info("{} {} {} {} {} {} {} {} {}".format(self.RefSeqTranscript,
            self.GeneName, self.GenomeAssembly, self.Chromosome, self.GeneStartPos,
            self.GeneEndPos, self.strand, self.SpacerSizes, self.PAMsize))


    # Utility function to validate input data
    def valdSheet(self, sheet, input_header, expected, cols):
        try:
            id = pd.read_excel(self.ifile, sheet_name=sheet)
        except:
            print("{} unable to read {} sheet".format(self.ifile, sheet))
            logging.error("Unable to read {} sheet".format(sheet))
            sys.exit(2)
        if len(id.index) == 0:
            print("Empty {} table".format(sheet))
            logging.warning("Empty {} table".format(sheet))
            sys.exit(0)

        id = id.iloc[:,0:cols]
        id.columns = [col.strip() for col in id.columns]
        id.columns = [s.lower() for s in id.columns]
        if not all([s in input_header for s in id.columns]):
            print("{}, expected header for {}: \n{}".format(self.ifile, sheet, expected))
            logging.error("Expected header for {}: \n{}".format(sheet, expected))
            sys.exit(2)
        for col in list(id.columns):
            try:
                id[col] = id[col].str.strip()
            except:
                pass
        logging.info("Successfully read {} sheet".format(sheet))
        return id


    # Utility function to extract fasta sequence in a given chromose, start, stop
    # Input: chromose, start, stop
    # Output: DNA sequence, string
    def getFASTA(self, chr, start, end):
        # TwoBitToFa path
        self.getSeqApp = cfg.TwoBitToFa
        # hg38 genome location
        self.getSeqGenome = cfg.Hg38_Genome
        filename = "{}/{}_chr{}_{}_{}.fa".format(self.wd, self.GeneName, chr,
                                                 start, end)
        cmd = "{} {} -seq=chr{} -start={} -end={} {}".format(self.getSeqApp,
                                   self.getSeqGenome, chr, start, end, filename)
        BasicFileUtils.processCommand(cmd)
        return FastaUtils.getFastaSeq(filename)[0]


    # Utility function to extract exon and intron information in given gene
    # Output: Exon and Intron location in doctonary iterator
    def getExonInfo(self):
        # Helper function to exit saftely in case of error
        # Local call only
        def rollback():
            mariadb_connection.rollback()
            sys.exit(2)

        try:
            # MariaDB access user information
            mariadb_connection = mariadb.connect(user=cfg.User,
                                                 password=cfg.Password,
                                                 database=cfg.Database,
                                                 host=cfg.Host)
        except Exception as e:
            print("Unable to connect to hg38 database")
            logging.error("Unable to connect to hg38 database\n{}".format(e))
            rollback()
        # Global query template
        query = "SELECT txStart, txEnd, name2, chrom, exonStarts, exonEnds \
                 FROM   {} AS A\
                 WHERE  A.name=\"{}\" order by chrom;"
        try:
            # Open database connection and execute query in specific transcript
            # and gene
            cursor = mariadb_connection.cursor()
            cursor.execute(query.format(cfg.Table, self.RefSeqTranscript.split('.')[0]))
        except Exception as e:
            print("Failed to retrieve exons information")
            logging.error("Failed to retrieve exons information\n{}".format(e))
            rollback()
        # Convert query output to list from iteratior exit if not found
        raw_info = list(cursor)
        if not any(True for _ in raw_info):
            print("Database does not contain exons information on:\
                   \n{} {}".format(self.RefSeqTranscript, self.GeneName))
            logging.error("Database does not contain exons information on:\
                   \n{} {}".format(self.RefSeqTranscript, self.GeneName))
            rollback()
        start = int(raw_info[0][0])
        end = int(raw_info[0][1])
        r_s = range(start-2, start+2, 1)
        r_e = range(end-2, end+2, 1)
        if self.GeneStartPos not in r_s or self.GeneEndPos not in r_e:
            print("Gene Start/End Position does not match database values [{},{}]".format(start,end))
            logging.error("Gene Start/End Position does not match database value\n \
                           Start: {}/{}\nEnd: {}/{}".format(start, self.GeneStartPos,
                                                end, self.GeneEndPos))
            sys.exit(2)
        if self.GeneName != raw_info[0][2]:
            print("Gene name does not match transcript\nCheck info sheet")
            logging.error("Gene name does not match transcript")
            sys.exit(2)
        # Manipulate output
        chr = raw_info[0][3]
        exonStarts = raw_info[0][4].split(',')[:-1]
        exonStarts = [int(x) for x in exonStarts]
        exonEnds = raw_info[0][5].split(',')[:-1]
        exonEnds = [int(x) for x in exonEnds]
        exons = zip(exonStarts, exonEnds)
        # Verify chromose number from input sheet
        if chr != "chr" + str(self.Chromosome):
            print("Mistake in chromosome")
            logging.error("Mistake in chromosome")
            rollback()

        mariadb_connection.rollback()
        logging.info("Successfully pulled exon info from hg38 database")
        return exons


    # Produce separte lists of exon and intron coordinates and adjust by strand
    def exonIntron(self):
        self.exons_raw = self.getExonInfo()
        exons = list(self.exons_raw)
        introns = []
        for num, exon in enumerate(exons):
            try:
                introns.append((exon[1] + 1, exons[num+1][0] - 1))
            except:
                pass
        self.downstream = exons[0][0]
        self.upstream = exons[-1][1]
        if self.strand == "-":
            exons.reverse()
            introns.reverse()
        self.exons = exons
        self.introns = introns


    # Utility function to locate mutation location in exon or intron proximity
    def mutPosition(self, exons, introns, position):
        if position < self.downstream:
            close_exon = abs(position - self.downstream) + 1
            if self.strand == "-":
                return "downstream", close_exon
            return "upstream", close_exon
        elif position > self.upstream:
            close_exon = abs(position - self.downstream) + 1
            if self.strand == "-":
                return "upstream", close_exon
            return "downstream", close_exon

        num_exons = len(exons)
        num_introns = len(introns)
        # Check if mutation in exon
        for num, exon in enumerate(exons, start=1):
            if position in range(exon[0], exon[1]+1):
                mut_pos = "Exon {} of {}".format(num, num_exons)
                close_exon = 0
                return mut_pos, close_exon
        # Check if mutation in intron and proximity
        for num, intron in enumerate(introns, start=1):
            if position in range(intron[0], intron[1]+1):
                mut_pos = "Intron {} of {}".format(num, num_introns)
                close_exon = min(abs(position - intron[0]),
                                 abs(position - intron[1]))
                return mut_pos, (close_exon + 1)

        print("Unable to locate mutation position\n \
               chr{}:{}".format(self.Chromosome, position))
        logging.error("Unable to locate mutation position\n \
                       chr{}:{}".format(self.Chromosome, position))
        sys.exit(2)


    def UpdatePatentsList(self, df_out, patlist):
        for index, row in df_out['PatentsPrefix'].dropna().iteritems():
            arr=row.split(',')
            for pat in arr:
                patlist.append(pat)

        patlist = list(dict.fromkeys(patlist))
        res = sorted(patlist, reverse=True)
        print("res len:" + str(len(res)))
        print(' '.join(map(str, res)))
        return res


    def FindSequenceInPatents(self, df):
        gRNA_header="Spacer"
        gRNAlen_header="Spacer.length"
        df['index_col'] = df.index

        filters=NucleotideUtils.getAllSequences('NNNN')
        colnames=[gRNA_header, 'Patents', 'PatentsPrefix']

        df_result = None
        for f in filters:
            for l in range(20,23):
                df_filter=df[(df[gRNA_header].str.startswith(f)) & (df[gRNAlen_header]==l)]
                if not df_filter.empty:
                    lens_file=cfg.lens_prefix+f+"_"+str(l)+".txt.tsv"
                    df_lens=pd.read_csv(lens_file, sep='\t', index_col=None, names=colnames)
                    df_filter = pd.merge(df_filter, df_lens, how='left', on=[gRNA_header])
                    df_result = pd.concat([df_result,df_filter], ignore_index=True)

        df=df_result.sort_values('index_col')
        df[['Patents','PatentsPrefix']]=df[['Patents','PatentsPrefix']].fillna('')
        df.drop('index_col', axis='columns', inplace=True)
        return df


    # Generate excel output from dataframe
    def writeInputOutput(self, old_df, old_sheet, new_df, new_sheet, empty_sheet):
        ofile = "".join((self.outdir, '/', self.filename, '_OUT_{}.xlsx'.format(curDate)))
        tmpfile = BasicFileUtils.precopyFile(ofile, self.wd, 'output')
        # Build QA dataframe
        QA = pd.DataFrame(self.TestGuides, columns=["TestingGuides"])
        QA["ExistInNonFilteredLists"] = np.nan
        QA["ExistInFilteredLists"] = np.nan
        # Helper function to search guides in frames
        def fnc(guide, df):
            if guide in df.Spacer.tolist():
                return "TRUE"
            else:
                return "FALSE"
        QA.ExistInNonFilteredLists = QA.TestingGuides.apply(fnc, args=[new_df])

        writer = pd.ExcelWriter(self.out)
        self.sheet1.to_excel(writer, sheet_name = 'Sheet1', index=False)
        self.info.to_excel(writer, sheet_name = 'Info', index=False)
        # read all existing sheets and write them back
        try:
            xlsx = pd.ExcelFile(self.out)
            for sheet in xlsx.sheet_names:
                if sheet in [new_sheet, 'Sheet1', 'Info', 'QA']:
                    continue
                df = xlsx.parse(sheet_name=sheet)
                df.to_excel(writer, sheet_name=sheet, index=False)
                if sheet.endswith("_out"):
                    for guide in self.TestGuides:
                        if guide in df.Spacer.tolist():
                            QA.loc[QA.TestingGuides==guide, "ExistInNonFilteredLists"] = "TRUE"

        except:
            pass

#        pat_summary = ' OR '.join([r.replace('_', '') for r in patlist])
#        patlist.append(pat_summary)
#        dfpat = pd.DataFrame(patlist, columns=['Patents'])
#        dfpat.to_excel(writer, sheet_name='PatentsList', index=False)
        old_df.to_excel(writer, sheet_name=old_sheet, index=False)
        new_df.to_excel(writer, sheet_name=new_sheet, index=False)
        QA.to_excel(writer, sheet_name='QA', index=False)
        empty = pd.DataFrame()
        empty.to_excel(writer, sheet_name=empty_sheet, index=False)
        writer.save()
        writer.close()
        logging.info("Successfully created excel output")


    # CasOffTargets with 0 mismatches for gRNAs that target a SNP or mutation
    def casOff(self, df):
        # CasOffTargets script directory
        casoff = cfg.CasOffinder
        # Hg38 database directory
        db = cfg.Hg38_DB
        # Temporary input and output directories
        pathInput = os.path.join(self.wd, "CasOffinder_0mm")
        pathOutput = os.path.join(self.wd, "CasOffinder_0mm_out")
        try:
            os.mkdir(pathInput)
            os.mkdir(pathOutput)
        except Exception as e:
            print("Unable to make casoff directories")
            logging.error("Unable to make casoff directories\n{}".format(e))
            sys.exit(2)

        df["OffTargets.0.mm"] = int(0)
        # replace(df, 'U', 'T')
        for col in ["PAM", "Target"]:
            df[col] = df[col].apply(lambda x: x.replace('U', 'T'))
        # Producing input files by dividing target seq to PAM and Spacer size
        df2 = df.copy(deep=True)
        df2.drop_duplicates(subset='Target', inplace=True)
        while not df2.empty:
            pam = df2["PAM"].values[0]
            length = df2["Spacer.length"].values[0]
            temp = df2.loc[(df2["PAM"] == pam) & (df2["Spacer.length"] == length)]
            file = [db, "N"*length + pam]
            targs = temp["Target"].tolist()
            file += [x + " 0" for x in targs]
            name = "casOffinder_{}_input_{}.txt".format(length, pam)
            path = "/".join((pathInput, name))
            with open(path, 'w') as f:
                f.writelines("%s\n" % x for x in file)
            df2.drop(temp.index.values.tolist(), inplace=True)
        # Run CasOffTarget
        try:
            cwd = os.getcwd()
            os.chdir(casoff)
            print("Runinng CasOffinder, this might take awhile...")
            cmd = "./OffTarget.sh {} {} {} {}".format(pathInput, pathOutput,
                                               self.user, self.mail)
            print(cmd)
            BasicFileUtils.processCommand(cmd)
#            os.system(cmd)
            inputCounter = len(glob.glob1(pathInput,"*.txt"))
            outputCounter = len(glob.glob1(pathOutput,"*.txt"))
            if inputCounter != outputCounter:
                print ("Command: {} failed. Exiting ...".format(cmd))
                sys.exit(2)
                cmd = "./OffTarget_serial.sh {} {} {} {}".format(pathInput, pathOutput,
                                               self.user, self.mail)
                os.system(cmd)
            os.chdir(cwd)
        except Exception as e:
            print("Failed to execute casOffinder")
            logging.error("Failed to execute casOffinder\n{}".format(e))
            sys.exit(2)
        inputCounter = len(glob.glob1(pathInput,"*.txt"))
        outputCounter = len(glob.glob1(pathOutput,"*.txt"))
        if inputCounter != outputCounter:
            print("No of input files " + str(inputCounter))
            print("No of output files " + str(outputCounter))
            print("No output generated for CasOffinder\n{}".format(pathOutput))
            logging.error("No output generated for CasOffinder\n{}".format(pathOutput))
            sys.exit(2)
        # Process output, increment OffTargets.0.mm
        data = pd.DataFrame(columns=['Target', 'count'])
        for f in os.listdir(pathOutput):
            with open("/".join((pathOutput, f)), "r") as f:
                try:
                    data_tmp = pd.read_csv(f, sep="\t", header=None, index_col=False)
                except pd.errors.EmptyDataError:
                    continue
                data_tmp.columns = ["Target", "Chr", "Location", "Target2", "Strand", "#"]
                data_tmp = data_tmp.groupby(["Target"]).size().reset_index(name='count')
                data = data.append(data_tmp, ignore_index=True)
        df = df.join(data.set_index('Target'), on='Target')
        df = df.fillna(0)
        df["OffTargets.0.mm"] = df["count"]
        df.drop(columns=["count"], inplace=True)

        for col in ["PAM", "Target"]:
            df[col] = df[col].apply(lambda x: x.replace('T', 'U'))

        try:
            shutil.rmtree(pathInput)
            shutil.rmtree(pathOutput)
        except Exception as e:
            print("Unable to delete casoff directories")
            logging.error("Unable to delete casoff directories\n{}".format(e))
            sys.exit(2)
        logging.info("Successfully ran casOffinder")
        return df


    # Upload a file to an S3 bucket
    # :param file_name: File to upload
    # :param bucket: Bucket to upload to
    # :param object_name: S3 object name.
    # :return: True if file was uploaded, else False
    def uploadAWS(self, file_name):
        odir = self.outdir
        odir = odir.replace('s3://', '')
        bucket = odir.split('/')[0]
        object_name = "".join((odir.replace(bucket + '/', ''), '/',
                              self.filename,'_OUT.xlsx'))
        # If S3 object_name was not specified, use file_name
        if object_name is None:
            object_name = file_name
        # Upload the file
        s3_client = boto3.client('s3')
        try:
            response = s3_client.upload_file(file_name, bucket, object_name)
        except Exception as e:
            print("Unable to upload output file to s3 bucket\n{}".format(e))
            sys.exit(2)
