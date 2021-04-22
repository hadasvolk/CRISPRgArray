import sys
import shutil
import logging
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from gRNA_Utils import gRNA_Utils
from FileUtils import BasicFileUtils, FastaUtils
from SequenceUtils import NucleotideUtils

"""
Class Guide RNA produces guides for Non-Discriminatory regions in input file
and inherent from genral gRNA_Utils class
Guides are generated for Spacer and PAM lenghts given by user input file.

Hadas Volkov 03/2020
"""

class gRNA_NonDiscriminatory(gRNA_Utils):
    def __init__ (self, ifile, file_ext, odir, tmpdir):
        # Apply inital gRNA class from super class gRNA_Utils
        logging.info("Init gRNA_NonDiscriminatory")
        super().__init__(ifile, file_ext, odir, tmpdir)
        # Acquire input from mutation sheet
        logging.info("Reading non-discriminatory sheet")
        self.readNDSheet()
        # Acquire exon and intron information for gene from 'refGene_38_copy' database
        logging.info("Pulling exon-intron information from database")
        self.exonIntron()
        # Build guides dataframe
        logging.info("Starting main gRNA process")
        nd_df = self.work()
        # casOffinder
        logging.info("Running CasOffinder")
        nd_df = self.adjustCasoff(self.casOff(nd_df.copy(deep=True)))
        # Generate excel output from dataframe
        logging.info("Writing excel output")
        self.writeInputOutput(self.nd, 'non-discriminatory_in', nd_df, 'non-discriminatory_out', 'non-discriminatory_filtered')


    def __del__ (self):
        # Remove short fasta sequences from temp directory
        try:
            BasicFileUtils.clearFa(self.wd)
            BasicFileUtils.postcopy(self.wd, self.outdir)
            shutil.rmtree(self.wd)
        except AttributeError:
            pass


    def readNDSheet(self):
        # Valdiate non-discriminatory sheet from input file and convert to dataframe
        expected = ["RegionName", "Chromosome", "RegionStartPosition",
            "RegionEndPosition", "ClosestExon"]
        self.input_header = [s.lower() for s in expected]
        self.output_header = ["gRNA.Name", "Source", "ND.Position", "Closest.Exon.dist",
            "Region.Name", "Spacer.length", "ND.5prime.dist", "Spacer", "PAM", "Target",
            "OffTargets.0.mm", "GC.percent", "Single.consecutive.base",
            "Double.consecutive.base", "Triple.consecutive.base"]

        self.nd = self.valdSheet('non-discriminatory', self.input_header, expected, 5)

        chr = [x for x in self.nd.chromosome.tolist() if str(x) != 'nan']
        for col in ["regionstartposition", "regionendposition", "closestexon"]:
            self.nd[col] = pd.to_numeric(self.nd[col])

        for index, row in self.nd.iterrows():
            # Sequence from start to end positionsq
            if row.regionstartposition >= row.regionendposition:
                logging.error("Mistake in row {} in region start/end positions".format(index))
                print("Mistake in row {} in region start/end positions".format(index+1))
                sys.exit(2)

        if self.Chromosome != int(chr[0]) or not all(x == chr[0] for x in chr):
            print("Mistake in chromosome {}".format(self.Chromosome))
            logging.error("Mistake in chromosome {}".format(self.Chromosome))
            sys.exit(2)


    def adjustCasoff(self, df):
        df["OffTargets.0.mm"] -= 1
        temp = df.loc[df["OffTargets.0.mm"] == -1]
        if not temp.empty:
            logging.warning("Guides not found by CasOffinder\n{}".format(temp))
            print("Guides not found by CasOffinder\n{}".format(temp))
        return df

    def work(self):
        # Main function to produce guides
        # Main dataframe to append guides
        df = pd.DataFrame(columns=self.output_header)

        for index, row in self.nd.iterrows():
            # Sequence from start to end positionsq
            regionStartPosition = row.regionstartposition-1
            fasta = self.getFASTA(self.Chromosome, regionStartPosition,
                                  row.regionendposition).replace('T', 'U')
            # Iterate ove spacers lenghts
            for spacer in self.SpacerSizes:
                logging.info("ND spacer length: {}".format(spacer))
                # Target length
                targ_len = spacer + self.PAMsize
                # Targets for foward strand
                targs = [fasta[i:i+targ_len] for i in range(len(fasta))]
                targs = [s for s in targs if len(s) == targ_len]
                spacers = [targ[:spacer] for targ in targs]
                pams = [targ[spacer:] for targ in targs]
                li = range(row.regionstartposition, row.regionendposition+1)
                tr = [self.mutPosition(self.exons, self.introns, i) for i in li]
                zipped = zip(range(len(targs)), spacers, pams, targs, li, tr)
                # Targets for reverse strand
                targs_rev = [str(NucleotideUtils.getReverseComplement(targ)) for targ in targs]
                targs_rev = [targ.replace('T', 'U') for targ in targs_rev]
                spacers_rev = [targ[:spacer] for targ in targs_rev]
                pams_rev = [targ[spacer:] for targ in targs_rev]
                li_rev = range(row.regionstartposition+targ_len-1, row.regionendposition+1)
                # tr_rev = [self.mutPosition(self.exons, self.introns, i) for i in li_rev]
                # tr_rev = [(x, (y - targ_len + 1)) for (x, y) in tr_rev_raw if not x.startswith('Exon')]
                zipped_rev = zip(range(len(targs_rev)), spacers_rev, pams_rev, targs_rev, li_rev, tr)
                # Generate current dataframe for foward and reverse
                init_cols = ["Region.Name", "Spacer", "PAM", "Target", "ND.5prime.dist", "Closest.Exon.dist"]
                cur = pd.DataFrame(zipped, columns=init_cols)
                cur_rev = pd.DataFrame(zipped_rev, columns=init_cols)
                # gRNA naming
                cur["Region.Name"] = cur.apply(lambda x: "{}:{}_{}_ND_fw_{}bp_{}".format(
                    self.Chromosome, x["ND.5prime.dist"], row.regionname, spacer,
                    x["Region.Name"]), axis=1)
                cur_rev["Region.Name"] = cur_rev.apply(lambda x: "{}:{}_{}_ND_rev_{}bp_{}".format(
                    self.Chromosome, x["ND.5prime.dist"]-targ_len+1,
                    row.regionname, spacer, x["Region.Name"]), axis=1)
                # Concat foward and reverse dataframes
                cur = cur.append(cur_rev, ignore_index=True)
                # Operations on both foward and reverse guides
                cur["ND.Position"], cur["Closest.Exon.dist"] = cur["Closest.Exon.dist"].str
                cur["GC.percent"] = cur["Target"].apply(NucleotideUtils.GCContent)
                cur["GC.percent"] = cur["GC.percent"].apply(lambda x: "{0:.2f}".format(x))
                repeats = ["Single.consecutive.base", "Double.consecutive.base",
                           "Triple.consecutive.base"]
                for idx, col in enumerate(repeats, start=1):
                    cur[col] = cur["Target"].apply(NucleotideUtils.getMaxConsecutiveRepeats, args=(idx,))
                    # cur[col] = cur[col].apply(lambda x: x-1)
                cur["Source"] = "ND"
                cur["Spacer.length"] = spacer
                cur["OffTargets.0.mm"] = 0
                cur["gRNA.Name"] = row.regionname + "__chr" + str(self.Chromosome) + ":" + str(row.regionstartposition) + "-" + str(row.regionendposition)
                cur = cur[self.output_header]
                # Append to main dataframe
                df = df.append(cur, ignore_index=True)

        df = NucleotideUtils.mapPams(df, self.PAMs, self.PAMsize)
        df.drop(columns=["ND.5prime.dist"], inplace=True)
        df = self.FindSequenceInPatents(df)
        return df
