export RODA_ID=tr_unstr
export RODA_PIPELINE=dge-rmd

# remind the start dir of the output, it is overwritten in roda_single
export DIR_START=$DIR_OUTPUT

# scripts for RNA-seq pipeline
export SCRIPT_STAR=/home/irsan/programs/star/STAR_2.3.0e.Linux_x86_64/STAR
export SCRIPT_HTSEQ=/home/irsan/programs/HTSeq-0.5.4p1/scripts/htseq-count
export SCRIPT_RSCRIPT=/usr/bin/Rscript
export SCRIPT_PROCESS_HTSEQ=$DIR_BASE/bin/scripts/processHTseqs.R
export SCRIPT_PT_ADDREADGROUP=$SCRIPT_PT/AddOrReplaceReadGroups.jar
export SCRIPT_GATK=$DIR_SCRIPTS/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar
export SCRIPT_JAVA=$DIR_SCRIPTS/jre1.7.0_25/bin/java

# data files for RNA-seq
export INDEX_STAR_HG19=$DIR_DATA/hg19_star_ensGene_splice_db
export ENSEMBL_GTF_HG19=$DIR_DATA/Homo_sapiens.GRCh37.65_v2.gtf

# command line options
export STRAND_DIRECTION=no #yes, no or reverse
export AMOUNT_MISMATCHES=10
export CHIM_SEGMENT_MIN=500

# command line arguments for the variant calling part
export ARG_BT_DEPTH=true
export FILE_EXONBED=$DIR_DATA/bed/Homo_sapiens.GRCh37.65_v3.bed
export FILE_TARGETAREAS=$DIR_DATA/bed/Homo_sapiens.GRCh37.65_v3.bed

# Restrict RNAseq variant calling to Ensembl exons
export ARG_TARGETAREAS="-L $DIR_DATA/bed/Homo_sapiens.GRCh37.65_v3_offset1.bed"

# allow gatk to work with rna-seq reads
export GATK_RNASEQ="-U ALLOW_N_CIGAR_READS"

export ARG_CLEANUP=true

