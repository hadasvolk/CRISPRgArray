import sys
import shutil
import logging
import pandas as pd
import numpy as np

from gRNA_Utils import gRNA_Utils
from FileUtils import BasicFileUtils, FastaUtils
from SequenceUtils import NucleotideUtils

"""
Class Guide RNA produces guides for Mutation table in input file and inherent
from genral gRNA_Utils class
Guides are generated for Spacer and PAM lenghts given by user input file.

Hadas Volkov 03/2020
"""

class gRNA_Mutation(gRNA_Utils):
    def __init__ (self, ifile, file_ext, odir, tmpdir):
        # Apply inital gRNA class from super class gRNA_Utils
        logging.info("Init gRNA_Mutation")
        super().__init__(ifile, file_ext, odir, tmpdir)
        # Acquire input from mutation sheet
        logging.info("Reading Mutation sheet")
        self.readMutSheet()
        # Acquire exon and intron information for gene from 'refGene_38_copy' database
        logging.info("Pulling exon-intron information from database")
        self.exonIntron()
        # Build guides dataframe
        logging.info("Starting main gRNA process")
        mut_df = self.work()
        # casOffinder
        logging.info("Running CasOffinder")
        mut_df = self.casOff(mut_df.copy(deep=True))
        # Generate excel output from dataframe
        logging.info("Writing excel output")
        self.writeInputOutput(self.mut, 'Mutations_in', mut_df, 'Mutations_out', 'Mutations_filtered')


    def __del__ (self):
        # Remove short fasta sequences from temp directory
        try:
            BasicFileUtils.clearFa(self.wd)
            BasicFileUtils.postcopy(self.wd, self.outdir)
            shutil.rmtree(self.wd)
        except AttributeError:
            pass


    def readMutSheet(self):
        # Valdiate mutation sheet from input file and convert to dataframe
        expected = ["Chromosome", "Position", "rsID", "Reference", "Alternate",
            "Source", "TranscriptConsequence"]
        self.input_header = [s.lower() for s in expected]
        self.output_header = ["gRNA.Name", "Info", "Source", "Mutation.Position",
            "Closest.Exon", "MUT.Name", "Spacer.length", "5prime.dist", "Spacer",
            "PAM", "Target", "OffTargets.0.mm", "GC.percent", "Single.consecutive.base",
            "Double.consecutive.base", "Triple.consecutive.base"]

        self.mut = self.valdSheet('Mutations', self.input_header, expected, 7)


    def work(self):
        # Main function to produce guides
        # Main dataframe to append guides
        df = pd.DataFrame(columns=self.output_header)
        # Iteration over each mutation in input file
        for index, row in self.mut.iterrows():
            chrm = str(row.chromosome)
            chrm = chrm.split('r')[-1]
            pos = int(row.position)
            ref = str(row.reference)
            len_ref = len(ref)
            alt = str(row.alternate)
            len_alt = len(alt)
            # Acuire mutation position by exon or intron
            mutpos, close_exon = self.mutPosition(self.exons, self.introns, pos)

            if ref != self.getFASTA(chrm, pos - 1, pos + len_ref -1):
                print("Reference not found in hg38;\n{}".format(row.T))
                logging.error("Reference not found in hg38;\n{}".format(row.T))
                sys.exit(2)

            length = max(self.SpacerSizes) + self.PAMsize + max(len_alt, len_ref)

            try:
                fastaseq_before = self.getFASTA(chrm, pos - length, pos - 1)
                fastaseq_after = self.getFASTA(chrm, pos + len_ref - 1,
                                               pos + len_ref + length)
            except Exception as e:
                print("Unable to get sequence from hg38 fasta")
                logging.error("Unable to get sequence from hg38 fasta\n{}".format(e))
                sys.exit(2)

            seq = fastaseq_before + alt + fastaseq_after
            seq = seq.replace('T', 'U')
            mut_pos = len(fastaseq_before)

            for spacer in self.SpacerSizes:
                logging.info("Mut: {} {} {} {} {}".format(chrm, pos, ref, alt, spacer))
                len_target = spacer + self.PAMsize
                for i in range(len_target):
                    target = seq[(mut_pos - i):(mut_pos + len_target - i)]
                    dist_5prime = i
                    guide = target[:spacer]
                    pam = target[spacer:]
                    gc = "{0:.2f}".format(NucleotideUtils.GCContent(guide))
                    repeats = []
                    for repeat in range(1, 4):
                        repeats.append(NucleotideUtils.getMaxConsecutiveRepeats(guide, repeat))
                    name = "{}:{}_{}_{}".format(chrm, pos, ref, alt)
                    gname = "{}:{}_{}_{}_MUT_fw_{}bp_{}".format(chrm, pos, ref,
                                                               alt, spacer,
                                                               dist_5prime)
                    fw = [name, row.transcriptconsequence, row.source, mutpos, close_exon, gname,
                          spacer, dist_5prime, guide, pam, target, np.nan, gc]
                    fw = fw + repeats

                    target_rev = str(NucleotideUtils.getReverseComplement(target))
                    target_rev = target_rev.replace('T', 'U')
                    dist_5prime_rev = len(target_rev) - i - 1
                    guide_rev = target_rev[:spacer]
                    pam_rev = target_rev[spacer:]
                    gc_rev = "{0:.2f}".format(NucleotideUtils.GCContent(guide_rev))
                    repeats_rev = []
                    for repeat in range(1, 4):
                        repeats_rev.append(NucleotideUtils.getMaxConsecutiveRepeats(guide_rev, repeat))
                    gname_rev = "{}:{}_{}_{}_MUT_rev_{}bp_{}".format(chrm, pos,
                                                            ref, alt, spacer,
                                                            dist_5prime_rev)
                    rev = [name, row.transcriptconsequence, row.source, mutpos, close_exon,
                           gname_rev, spacer, dist_5prime_rev, guide_rev, pam_rev,
                           target_rev, np.nan, gc_rev]
                    rev = rev + repeats_rev

                    df_temp = pd.DataFrame([fw, rev], columns=self.output_header)
                    df = df.append(df_temp, ignore_index=True)

        df = NucleotideUtils.mapPams(df, self.PAMs, self.PAMsize)
        df = self.FindSequenceInPatents(df)
        return df
