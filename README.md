# gRNA_IP_Step1
#
# Generating gRNAs for mutations, SNPs and non-discriminatory ranges
#

"""
gRNA4Patents step1 produce guides for Mutation, SNPs and non-discriminatory
tables in an input Excel file.
All listings refer to  the human reference genome GRCh38 assembly

Outputs excel sheet file for each table and QA for testing guides listed in Info sheet

Hadas Volkov - AppliedGenomics 04/2020
"""

Usage: ./gRNAs_Step1.sh -i [] -o [] -t [all/mut/snp/nd]
        -i Input xlsx file on s3
        -o output s3 directory
        -t Description of what gRNAs to produce

Mandatory Sheets and Headers for Input Excel File:
  Sheet1 - "user", "email", "purpose", "service"
  Info - "GeneName", "Chromosome", "RefSeqTranscript",
      "GenomeAssembly", "GeneStartPos", "GeneEndPos", "Strand", "SpacerSizes",
      "PAMsize", "TestingGuides", "AdditionalPAMs"
  Mutations - "Chromosome", "Position", "rsID", "Reference",
      "Alternate", "Source", "TranscriptConsequence"
  SNPs - "Chromosome", "SNPStartPosition", "SNPEndPos", "MinHetFreq"
  non-discriminatory - "RegionName", "Chromosome", "RegionStartPosition",
      "RegionEndPosition", "ClosestExon"

cfg.py - Paths and Headers format

Requriments:
  Python 3.6.10
    pandas>=1.0.1
    numpy>=1.18.1
    mysql-connector-python>=8.0.18
    boto3>=1.12.12
    botocore>=1.15.12
    pytabix>=0.0.2
    biopython>=1.76
