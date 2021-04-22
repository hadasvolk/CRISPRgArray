# CRISPRgArray 
#
# Generating gRNAs for mutations, SNPs and non-discriminatory ranges
#

"""<br/>
Produce guides for Mutation, SNPs and non-discriminatory
tables in an input Excel file.
All listings refer to  the human reference genome GRCh38 assembly

Outputs excel sheet file for each table and QA for testing guides listed in Info sheet

Hadas Volkov - AppliedGenomics 04/2020<br/>
"""

Usage: ./gRNAs_Step1.sh -i [] -o [] -t [all/mut/snp/nd]<br/>
        -i Input xlsx file on s3<br/>
        -o output s3 directory<br/>
        -t Description of what gRNAs to produce<br/>

Mandatory Sheets and Headers for Input Excel File:<br/>
  Sheet1 - "user", "email", "purpose", "service"<br/>
  Info - "GeneName", "Chromosome", "RefSeqTranscript",<br/>
      "GenomeAssembly", "GeneStartPos", "GeneEndPos", "Strand", "SpacerSizes",<br/>
      "PAMsize", "TestingGuides", "AdditionalPAMs"<br/>
  Mutations - "Chromosome", "Position", "rsID", "Reference",<br/>
      "Alternate", "Source", "TranscriptConsequence"<br/>
  SNPs - "Chromosome", "SNPStartPosition", "SNPEndPos", "MinHetFreq"<br/>
  non-discriminatory - "RegionName", "Chromosome", "RegionStartPosition",<br/>
      "RegionEndPosition", "ClosestExon"<br/>

cfg.py - Paths and Headers format<br/>

Requriments:<br/>
  Python 3.6.10<br/>
    pandas>=1.0.1<br/>
    numpy>=1.18.1<br/>
    mysql-connector-python>=8.0.18<br/>
    boto3>=1.12.12<br/>
    botocore>=1.15.12<br/>
    pytabix>=0.0.2<br/>
    biopython>=1.76<br/>
