import sys
import os
import shutil
import pandas as pd
import numpy as np
import tabix
import logging

from gRNA_Utils import gRNA_Utils
from FileUtils import BasicFileUtils, FastaUtils
from SequenceUtils import NucleotideUtils
import cfg

"""
Class Guide RNA produces guides for SNPs locations in input file and inherent
from genral gRNA_Utils class
Guides are generated for Spacer and PAM lenghts given by user input file.

Hadas Volkov 03/2020
"""

class gRNA_SNP(gRNA_Utils):
    def __init__ (self, ifile, file_ext, odir, tmpdir):
        # Apply inital gRNA class from super class gRNA_Utils
        logging.info("Init gRNA_SNP")
        super().__init__(ifile, file_ext, odir, tmpdir)
        # Acquire input from mutation sheet
        logging.info("Reading SNP sheet")
        self.readSnpSheet()
        # Acquire exon and intron information for gene from 'refGene_38_copy' database
        logging.info("Pulling exon-intron information from database")
        self.exonIntron()
        # Build guides dataframe
        logging.info("Starting main gRNA process")
        snp_df = self.work()
        # casOffinder
        logging.info("Running CasOffinder")
        snp_df = self.adjustCasoff(self.casOff(snp_df.copy(deep=True)))
        # Generate excel output from dataframe
        logging.info("Writing excel output")
        self.writeInputOutput(self.snp, 'SNPs_in', snp_df, 'SNPs_out', 'SNPs_filtered')


    def __del__ (self):
        # Remove short fasta sequences from temp directory
        try:
            BasicFileUtils.clearFa(self.wd)
            BasicFileUtils.postcopy(self.wd, self.outdir)
            shutil.rmtree(self.wd)
        except AttributeError:
            pass


    def readSnpSheet(self):
        # Valdiate mutation sheet from input file and convert to dataframe
        expected = ["Chromosome", "SNPStartPosition", "SNPEndPosition", "MinHetFreq"]
        self.input_header = [s.lower() for s in expected]
        self.output_header = ["gRNA.Name", "HetFreq V3", "HetFreq V2.1", "Source",
            "SNP.Position", "Closest.Exon", "SNP.Name", "Spacer.length",
            "SNP.5prime.dist", "Spacer", "PAM", "Target", "OffTargets.0.mm", "GC.percent",
            "Single.consecutive.base", "Double.consecutive.base", "Triple.consecutive.base"]

        self.snp = self.valdSheet('SNPs', self.input_header, expected, 4)

        chrom = [x for x in self.snp.chromosome.tolist() if str(x) != 'nan']
        if self.Chromosome != int(chrom[0]) or not all(x == chrom[0] for x in chrom):
            print("Mistake in chromosome {}".format(self.Chromosome))
            logging.error("Mistake in chromosome {}".format(self.Chromosome))
            sys.exit(2)

        self.Source = "gnomAD hg38"


    def adjustCasoff(self, df):
        df.loc[df['gRNA.Name'].str.contains('REF'), "OffTargets.0.mm"] -= 1
        temp = df.loc[df["OffTargets.0.mm"] == -1]
        if not temp.empty:
            logging.warning("Guides not found by CasOffinder\n{}".format(temp))
            print("Guides not found by CasOffinder\n{}".format(temp))
        return df


    def produceRecords(self, gnomad_dir, genome, Chrom, Start, End, hetFreq):
        try:
            tb = tabix.open(gnomad_dir + genome)
            records = tb.query("chr{}".format(Chrom), Start, End)
        except Exception as e:
            print("Unable to access genome database\n{}".format(e))
            logging.error("Unable to access genome database\n{}".format(e))
            sys.exit(2)
        v = dict()
        for record in records:
            infoDict = {}
            for r in range(7):
                infoDict[r] = record[r]
            infos = record[7].split(';')
            for i, info in enumerate(infos, 1):
                try:
                    key, value = info.split('=')
                except ValueError:
                    key = 'INFO{}'.format(i)
                    value = info
                infoDict[key] = value
            chrom, pos, rsID = (infoDict[0].split('r'))[1], int(infoDict[1]), infoDict[2]
            ref, alt = infoDict[3], infoDict[4]

            AC, AN, nhomalt = int(infoDict['AC']), int(infoDict['AN']), int(infoDict['nhomalt'])
            het_freq = (AC-2*nhomalt)/AN*2
            if het_freq < hetFreq:
                continue
            het_freq = "{0:.4f}".format(het_freq)

            # # VEP information
            # type = infoDict['variant_type']
            # vep = infoDict['vep'].split('|')
            # info = "{} {} {}".format(vep[3], vep[1], vep[6])

            v["{}{}{}{}".format(chrom, pos, ref, alt)] = [chrom, pos, rsID, ref,
                                                        alt, het_freq]
        return v

    # Main function to produce guides
    def work(self):
        # Main dataframe to append guides
        df = pd.DataFrame(columns=self.output_header)

        gnomad_dir_v3 = cfg.GnomAD_v3_Dir
        genome_v3 = cfg.GnomAD_v3_VCF.format(self.Chromosome)

        gnomad_dir_v2 = cfg.GnomAD_v21_Dir
        genome_v2 = cfg.GnomAD_v21_VCF.format(self.Chromosome)

        for index, row in self.snp.iterrows():
            Start = int(self.snp.loc[index, 'snpstartposition'])
            End = int(self.snp.loc[index, 'snpendposition'])
            hetFreq = float(self.snp.loc[index, 'minhetfreq'])
            v3 = self.produceRecords(gnomad_dir_v3, genome_v3, self.Chromosome,
                                     Start, End, hetFreq)
            v2 = self.produceRecords(gnomad_dir_v2, genome_v2, self.Chromosome,
                                     Start, End, hetFreq)

            merged = list(set(v2.keys()).union(set(v3.keys())))
            merged.sort()
            for snp in merged:
                if snp in v3.keys() and snp in v2.keys():
                    chrom, pos, ref, alt = v3[snp][0], v3[snp][1], v3[snp][3], v3[snp][4]
                    het_freq_v3, het_freq_v2 = v3[snp][5], v2[snp][5]
                    if v3[snp][2] != '.':
                        rsID = v3[snp][2]
                    else:
                        rsID = v2[snp][2]
                elif snp in v3.keys():
                    chrom, pos, ref, alt = v3[snp][0], v3[snp][1], v3[snp][3], v3[snp][4]
                    het_freq_v3, het_freq_v2, rsID = v3[snp][5], '--', v3[snp][2]
                else:
                    chrom, pos, ref, alt = v2[snp][0], v2[snp][1], v2[snp][3], v2[snp][4]
                    het_freq_v3, het_freq_v2, rsID = '--', v2[snp][5], v2[snp][2]

                logging.info("SNP: {} {} {} {} {}".format(rsID, chrom, pos, ref, alt))
                len_ref, len_alt = len(ref), len(alt)
                mutpos, close_exon = self.mutPosition(self.exons, self.introns, pos)

                if ref != self.getFASTA(chrom, pos - 1, pos + len_ref -1):
                    print("Reference not found in hg38;\n{} {} {} {}".format(chrom, pos, ref, alt))
                    logging.error("Reference not found in hg38;\n{} {} {} {}".format(chrom, pos, ref, alt))
                    sys.exit(2)

                length = max(self.SpacerSizes) + self.PAMsize + max(len_alt, len_ref)

                try:
                    fastaseq_before = self.getFASTA(chrom, pos - length, pos - 1)
                    fastaseq_after = self.getFASTA(chrom, pos + len_ref - 1,
                                                   pos + len_ref + length)
                except Exception as e:
                    print("Unable to get sequence from hg38 fasta\n{}".format(e))
                    logging.error("Unable to get sequence from hg38 fasta\n{}".format(e))
                    sys.exit(2)

                seq_ref = fastaseq_before + ref + fastaseq_after
                seq_alt = fastaseq_before + alt + fastaseq_after
                mut_pos = len(fastaseq_before)

                sequences = [(seq_ref, 'ref', 'REF'), (seq_alt, 'alt', 'SNP')]
                for seq, strand, r in sequences:
                    seq = seq.replace('T', 'U')
                    for spacer in self.SpacerSizes:
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
                            gname = "{}:{}_{}_{}_SNP_{}_fw_{}bp_{}".format(chrom, pos, ref,
                                    alt, strand, spacer, dist_5prime)
                            name = "{}:{}_{}_{}_{}_{}".format(chrom, pos, ref, alt,
                                                              rsID, r)
                            fw = [name, het_freq_v3, het_freq_v2, self.Source, mutpos, close_exon, gname,
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
                            gname_rev = "{}:{}_{}_{}_SNP_{}_rev_{}bp_{}".format(chrom, pos,
                                        ref, alt, strand, spacer, dist_5prime_rev)
                            rev = [name, het_freq_v3, het_freq_v2, self.Source, mutpos, close_exon,
                                   gname_rev, spacer, dist_5prime_rev, guide_rev, pam_rev,
                                   target_rev, np.nan, gc_rev]
                            rev = rev + repeats_rev

                            df_temp = pd.DataFrame([fw, rev], columns=self.output_header)
                            df = df.append(df_temp, ignore_index=True)
        df = NucleotideUtils.mapPams(df, self.PAMs, self.PAMsize)
        df = self.FindSequenceInPatents(df)
        return df
