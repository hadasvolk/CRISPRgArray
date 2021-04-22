import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO
import re
import itertools

class ProteinUtils:

    def validProteinSequence(infastafile):
        records = list(SeqIO.parse(infastafile, "fasta"))
        valid = 'ACTG'
        for letter in records[0].seq:
            if letter not in valid:
                return True
        return False

class NucleotideUtils:

    def validNucleotideSequence(infastafile):
        records = list(SeqIO.parse(infastafile, "fasta"))
        valid = 'ACTG'
        for letter in records[0].seq:
            if letter not in valid:
                return False
        return True


    def GCContent(seq):
        gc_count = seq.count("G") + seq.count("C")
        gc_fraction = float(gc_count) / len(seq)
        return 100 * gc_fraction


    def getMaxConsecutiveRepeats(seq, unitlen):
        final = 0
        for j in range(0, unitlen):
            substrs = [seq[i:i+unitlen] for i in range(j, len(seq), unitlen)]
            if len(substrs[-1]) != unitlen:
                substrs.pop()
            res = [list(i) for j, i in itertools.groupby(substrs)]
            res_ = []
            skip = False
            for s in res:
                for i in range(len(s[0])-1):
                    if s[0][i] == s[0][i+1]:
                        skip = True
                if not skip:
                    res_.append(len(s))
                skip = False
            try:
                result = max(res_)
                final = max(result, final)
            except:
                final = 0
        return final


    def getReverseComplement(seq):
        if isinstance(seq, str):
            seq = Seq(seq)
        return seq.reverse_complement()


    def mapPams(df, pams, pamsize):
        pamsDict = {}

        def fnc(pam, cpam):
            for num, nuc in enumerate(list(pam)):
                if nuc not in pamsDict[cpam][num]:
                    return 'FALSE'
            return 'TRUE'

        for pam in pams:
            nuc_list = []
            if len(pam) < pamsize:
                pam_m = pam + 'N'*(pamsize - len(pam))
            elif len(pam) >= pamsize:
                pam_m = pam[:pamsize]
            for nuc in pam_m:
                if nuc == 'N':
                    nuc_list.append(['A', 'T', 'C', 'G', 'U'])
                elif nuc == 'R':
                    nuc_list.append(['A', 'G'])
                elif nuc == 'Y':
                    nuc_list.append(['C', 'T', 'U'])
                elif nuc == 'W':
                    nuc_list.append(['A', 'T', 'U'])
                elif nuc == 'T':
                    nuc_list.append(['T', 'U'])
                else:
                    nuc_list.append([nuc])
            pamsDict[pam] = nuc_list

        for pam in pamsDict:
            df[pam] = df.PAM.apply(fnc, args=[pam])

        return df

    def getAllSequences(sequence):
        nuc_list = []
        for nuc in sequence:
            if nuc == 'N':
                nuc_list.append(['A', 'C', 'G', 'U'])
            elif nuc == 'R':
                nuc_list.append(['A', 'G'])
            elif nuc == 'Y':
                nuc_list.append(['C', 'U'])
            elif nuc == 'S':
                nuc_list.append(['C', 'G'])
            elif nuc == 'K':
                nuc_list.append(['G', 'U'])
            elif nuc == 'M':
                nuc_list.append(['A', 'C'])
            elif nuc == 'B':
                nuc_list.append(['C', 'G', 'U'])
            elif nuc == 'D':
                nuc_list.append(['A', 'G', 'U'])
            elif nuc == 'H':
                nuc_list.append(['A', 'C', 'U'])
            elif nuc == 'V':
                nuc_list.append(['A', 'C', 'G'])
            elif nuc == 'W':
                nuc_list.append(['A', 'U'])
            elif nuc == 'T':
                nuc_list.append(['U'])
            else:
                nuc_list.append([nuc])
        seq_lists = list(itertools.product(*nuc_list))
        return list(map(''.join, seq_lists))

