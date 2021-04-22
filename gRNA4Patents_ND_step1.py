import pandas as pd
import sys
import argparse
import tempfile
import os
import shutil
import logging
import datetime

curdir=os.path.dirname(os.path.realpath(__file__))
sys.path.append(curdir+'/libs')
assert sys.version_info >= (3, 6)

from gRNA_NonDiscriminatory import gRNA_NonDiscriminatory
from FileUtils import BasicFileUtils

curDate = datetime.datetime.now().strftime("%d-%m-%Y")

def main():
    # Reported to console only
    parser = argparse.ArgumentParser(description="gRNAs4IP_Step1",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-is3xlsx", required=True, help="Input xlsx file on s3")
    parser.add_argument("-s3odir", required=True, help="output s3 directory")

    args = parser.parse_args()

    ifile = args.is3xlsx
    if not ifile.startswith('s3://emendobio'):
        print ("Input file {} is not in emendo s3 bucket. \
                Exiting ...".format(ifile))
        sys.exit()
    if not (ifile.endswith('xlsx') or ifile.endswith('xls')):
        print ("Input file {} is not excel file. \
                Exiting ...".format(ifile))
        sys.exit()

    ofile = args.s3odir
    if not ifile.startswith('s3://emendobio'):
        print ("Output directiry {} is not in emendo s3 bucket. \
                Exiting ...".format(ofile))
        sys.exit()

    BasicFileUtils.fileExist(ifile)

    # Reported in consloe and in log file
    tmpdir = tempfile.mkdtemp(dir="/efs/tmp")
    os.chmod(tmpdir, 0o2775)
    BasicFileUtils.create_log(tmpdir, curDate, 'ND')

    print("Producing guides for NonDiscriminatory ...")

    gRNA_non_discriminatory = gRNA_NonDiscriminatory(ifile, '.xlsx', ofile, tmpdir)

    print("Finished succesfully !!!")

if __name__ == "__main__":
    main()
