export RODA_ID=tru_str
export RODA_PIPELINE=dge_snv

# remind the start dir of the output, it is overwritten in roda_single
export DIR_START=$DIR_OUTPUT

# scripts for RNA-seq pipeline
export SCRIPT_STAR=$DIR_SCRIPTS/STAR_2.3.0e.Linux_x86_64_static/STAR
export SCRIPT_HTSEQ=$DIR_SCRIPTS/HTSeq-0.6.1/scripts/htseq-count
export SCRIPT_RSCRIPT=/usr/bin/Rscript
export SCRIPT_PROCESS_HTSEQ=$DIR_BASE/bin/scripts/process_irna.R
export SCRIPT_PT_ADDREADGROUP=$SCRIPT_PT/AddOrReplaceReadGroups.jar
export SCRIPT_GATK=$DIR_SCRIPTS/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar
export SCRIPT_JAVA=/usr/bin/java

# data files for RNA-seq
export INDEX_STAR_HG19=$DIR_DATA/reference-genomes/hg19/STAR_2.3.0e/STAR-hg19-refGene20nov2014-sjdboverhang99/
export GTF=$DIR_DATA/gene-annotation/refGene-21Nov2014.gtf

# command line options
export STRAND_DIRECTION=reverse #yes, no or reverse
export AMOUNT_MISMATCHES=10
export CHIM_SEGMENT_MIN=1
export FEATURE_TYPE='exon'
export FILTER_INTRON_MOTIFS='None'
export HTSEQ_COUNT_MODE="union"

# command line arguments for the variant calling part
export ARG_BT_DEPTH=true
export FILE_EXONBED=$DIR_DATA/gene-annotation/refGene-21Nov2014-merged.bed
export FILE_TARGETAREAS=$DIR_DATA/gene-annotation/refGene-21Nov2014-merged-plusmin5.bed

# Restrict RNAseq variant calling to Ensembl exons
export ARG_TARGETAREAS="-L $DIR_DATA/gene-annotation/refGene-21Nov2014-merged-plusmin5.bed"

# allow gatk to work with rna-seq reads
export GATK_RNASEQ="-U ALLOW_N_CIGAR_READS"

export ARG_CLEANUP=true

