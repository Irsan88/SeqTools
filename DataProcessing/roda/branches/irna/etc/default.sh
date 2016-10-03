# ---------------------------------------------------------------------------- #
# Scripts and tools
# ---------------------------------------------------------------------------- #
export DIR_SCRIPTS=/bioinformatics/programs

export SCRIPT_JAVA=/usr/bin/java
export SCRIPT_PYTHON=/usr/bin/python
export SCRIPT_BWA=$DIR_SCRIPTS/bwa-0.7.9a/bwa
export SCRIPT_SAMTOOLS=$DIR_SCRIPTS/samtools-0.1.18/samtools
export SCRIPT_GATK=$DIR_SCRIPTS/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar

export SCRIPT_PT=$DIR_SCRIPTS/picard-tools-1.115/
export SCRIPT_PT_SORTSAM=$SCRIPT_PT/SortSam.jar
export SCRIPT_PT_REORDER=$SCRIPT_PT/ReorderSam.jar
export SCRIPT_PT_MARKDUP=$SCRIPT_PT/MarkDuplicates.jar
export SCRIPT_PT_FIXMATE=$SCRIPT_PT/FixMateInformation.jar

export SCRIPT_AV=$DIR_SCRIPTS/annovar-Jul2014
export SCRIPT_AV_CONVERT=$SCRIPT_AV/convert2annovar.pl
export SCRIPT_AV_SUMMARIZE=$SCRIPT_AV/summarize_annovar.pl

export SCRIPT_BT=$DIR_SCRIPTS/bedtools2-2.20.1/
export SCRIPT_BT_COVBED=$SCRIPT_BT/bin/coverageBed

export SCRIPT_FIXAF=$DIR_BASE/bin/scripts/calc_AF_VCF.py
export SCRIPT_CALCTOWERS=$DIR_BASE/bin/scripts/calcTowers.py


# ---------------------------------------------------------------------------- #
# Data files
# ---------------------------------------------------------------------------- #
export DIR_DATA=/bioinformatics/library
export FILE_REFERENCE=$DIR_DATA/reference-genomes/hg19/bwa-0.7.9a-r786/hg19.fa
export FILE_DBSNP=$DIR_DATA/variation/dbsnp_138.hg19.vcf
export FILE_KNOWNINDELS=$DIR_DATA/variation/Mills_and_1000G_gold_standard.indels.hg19.vcf
export FILE_TARGETAREAS=""
export FILE_ANNOVARDB=$SCRIPT_AV/humandb



# ---------------------------------------------------------------------------- #
# Arguments used in all steps
# ---------------------------------------------------------------------------- #
export ARG_BWA_OPTIONS=""
export ARG_BWA_THREADS=8 # Bit pointless due to retarded SGE functionality
export ARG_GATK_THREADS=8 # Again, SGE stupidity but we try...
export ARG_AV_DBSNP=138
export ARG_CLEANUP=true

export ARG_GATK_BED=false # Use bed file in GATK steps to focus on target areas
export ARG_ST_DEPTH=false # Use sam tools to provide information on coverage depths (contis)
export ARG_BT_DEPTH=false # Use bed tools to provide information on coverage depths (exomes)

export RODA_ID=rd
export RODA_PIPELINE=snpindel
