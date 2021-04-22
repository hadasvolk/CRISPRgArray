
# gRNA4Patents paths configuration file
# TwoBitToFa path
TwoBitToFa = "/efs/emendo/dev/utils/bin/twoBitToFa"
# hg38 genome location - small enough to remain on EFS
Hg38_Genome = "/efs/database/UCSC/hg38.2bit"

# MariaDB access user information
#User = 'genome'
#Password = ''
#Database = 'hg38'
#Host = 'genome-mysql.cse.ucsc.edu'
#Table = 'refGene'
#User = 'root'
User = 'genome'
#Password = 'emendobio123'
Password = ''
#Database = 'emendobio'
Database = 'hg38'
#Host = 'devsql.cluster-cxu1aid2cajt.eu-west-1.rds.amazonaws.com'
Host = 'genome-mysql.cse.ucsc.edu'
#Table = 'refGene_38_220420'
Table = 'refGene'

# CasOffTargets script directory
CasOffinder = "/efs/emendo/dev/IP/OffTarget/"
# Hg38 database directory
Hg38_DB = "/bigdisk/database/human/hg38_casOffinder"
# Lens nucleotide file prefix
lens_prefix="/efs/database/Lens/Nucleotides/all__application__na-all/gRNA_patents_"

# GnomAD
GnomAD_v3_Dir = "/efs/database/gnomAD_hg38/genomes/gnomad_v3/"
GnomAD_v3_VCF = "gnomad.genomes.r3.0.sites.chr{}.vcf.bgz"
GnomAD_v21_Dir = "/efs/database/gnomAD_hg38/genomes/gnomad_v2_161019/"
GnomAD_v21_VCF = "gnomad.genomes.r2.1.1.sites.{}.liftover_grch38.vcf.bgz"
