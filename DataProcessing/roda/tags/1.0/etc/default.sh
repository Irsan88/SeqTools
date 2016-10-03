# ---------------------------------------------------------------------------- #
# Scripts and tools
# ---------------------------------------------------------------------------- #
export DIR_SCRIPTS=/illumina/diagnostics/bin

export SCRIPT_JAVA=$DIR_SCRIPTS/jre1.6.0_38/bin/java
export SCRIPT_PYTHON=python2.7
export SCRIPT_BWA=$DIR_SCRIPTS/bwa-0.6.2/bwa
export SCRIPT_SAMTOOLS=$DIR_SCRIPTS/samtools-0.1.18/samtools
export SCRIPT_GATK=$DIR_SCRIPTS/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar

export SCRIPT_PT=$DIR_SCRIPTS/picard-tools-1.81
export SCRIPT_PT_SORTSAM=$SCRIPT_PT/SortSam.jar
export SCRIPT_PT_REORDER=$SCRIPT_PT/ReorderSam.jar
export SCRIPT_PT_MARKDUP=$SCRIPT_PT/MarkDuplicates.jar
export SCRIPT_PT_FIXMATE=$SCRIPT_PT/FixMateInformation.jar

export SCRIPT_AV=$DIR_SCRIPTS/annovar_121101
export SCRIPT_AV_CONVERT=$SCRIPT_AV/convert2annovar.pl
export SCRIPT_AV_SUMMARIZE=$SCRIPT_AV/summarize_annovar.pl

export SCRIPT_BT=$DIR_SCRIPTS/bedtools-2.17.0/bin
export SCRIPT_BT_COVBED=$SCRIPT_BT/coverageBed



# ---------------------------------------------------------------------------- #
# Data files
# ---------------------------------------------------------------------------- #
export DIR_DATA=/illumina/diagnostics/lib
#export FILE_REFERENCE=$DIR_SCRIPTS/hg19/hg19.fa
#export FILE_REFERENCE=$DIR_DATA/hg19gatk/ucsc.hg19.fasta
export FILE_REFERENCE=$DIR_DATA/hg19daud/ucsc.daud.hg19.fasta
export FILE_DBSNP=$DIR_DATA/dbsnp/dbsnp_137.hg19.vcf
export FILE_KNOWNINDELS=$DIR_DATA/knownindels/Mills_and_1000G_gold_standard.indels.hg19.vcf
export FILE_TARGETAREAS=""
export FILE_ANNOVARDB=$DIR_DATA/humandb2012



# ---------------------------------------------------------------------------- #
# Arguments used in all steps
# ---------------------------------------------------------------------------- #
export ARG_BWA_THREADS=8 # Bit pointless due to retarded SGE functionality
export ARG_GATK_THREADS=8 # Again, SGE stupidity but we try...
export ARG_AV_DBSNP=137
export ARG_CLEANUP=true

export ARG_GATK_BED=false # Use bed file in GATK steps to focus on target areas
export ARG_ST_DEPTH=false # Use sam tools to provide information on coverage depths (contis)
export ARG_BT_DEPTH=false # Use bed tools to provide information on coverage depths (exomes)

export RODA_ID=rd
export RODA_PIPELINE=snpindel
