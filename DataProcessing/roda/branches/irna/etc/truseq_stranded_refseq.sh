export RODA_ID=tru_str
export RODA_PIPELINE=dge

# scripts for RNA-seq pipeline
export SCRIPT_STAR=$DIR_SCRIPTS/STAR_2.3.0e.Linux_x86_64_static/STAR
export SCRIPT_HTSEQ=$DIR_SCRIPTS/HTSeq-0.6.1/scripts/htseq-count
export SCRIPT_PT_ADDREADGROUP=$SCRIPT_PT/AddOrReplaceReadGroups.jar

# data files for RNA-seq
export INDEX_STAR_HG19=$DIR_DATA/reference-genomes/hg19/STAR_2.3.0e/STAR-hg19-refGene20nov2014-sjdboverhang99
export GTF=$DIR_DATA/gene-annotation/refGene-21Nov2014.gtf

# RNA-seq options
export STRAND_DIRECTION=reverse #yes, no or reverse
export AMOUNT_MISMATCHES=10
export CHIM_SEGMENT_MIN=1
export FEATURE_TYPE='exon'
export FILTER_INTRON_MOTIFS='None'
export HTSEQ_COUNT_MODE="union"

# clean up intermediate files
export ARG_CLEANUP=true

