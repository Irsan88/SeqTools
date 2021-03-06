export RODA_ID=tru_unstr
export RODA_PIPELINE=dge

# scripts for RNA-seq pipeline
export SCRIPT_STAR=$DIR_SCRIPTS/STAR_2.3.0e.Linux_x86_64_static/STAR
export SCRIPT_HTSEQ=$DIR_SCRIPTS/HTSeq-0.6.1/scripts/htseq-count
export SCRIPT_PT_ADDREADGROUP=$SCRIPT_PT/AddOrReplaceReadGroups.jar

# data files for RNA-seq
export INDEX_STAR_HG19=$DIR_DATA/reference-genomes/DDX11-only/STAR
export ENSEMBL_GTF_HG19=$DIR_DATA/gene-annotation/ensembl-genes-rnaseq/Homo_sapiens.GRCh37.75.chr1-22XY.gtf

# command line options
export STRAND_DIRECTION=no #yes, no or reverse
export AMOUNT_MISMATCHES=10
export CHIM_SEGMENT_MIN=1
export FILTER_INTRON_MOTIFS='None'
export FEATURE_TYPE=exon

# clean up intermediate files
export ARG_CLEANUP=true

