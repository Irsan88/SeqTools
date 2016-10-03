export RODA_ID=tru_str
export RODA_PIPELINE=dge_snv

# scripts for RNA-seq pipeline
export SCRIPT_STAR=/bioinformatics/programs/STAR_2.3.0e.Linux_x86_64_static/STAR
export SCRIPT_HTSEQ=/bioinformatics/programs/HTSeq-0.6.1/scripts/htseq-count
export SCRIPT_PT_ADDREADGROUP=$SCRIPT_PT/AddOrReplaceReadGroups.jar

# data files for RNA-seq
export INDEX_STAR_HG19=$DIR_DATA/reference-genomes/hg19/STAR_2.3.0e/STAR-hg19-ensembl75-sjdboverhang99/
export GTF=$DIR_DATA/ensemble-genes/Homo_sapiens.GRCh37.75.onlyCDSCodonUTR.chr1-22XY.gtf

# command line options
export STRAND_DIRECTION=reverse #yes, no or reverse
export AMOUNT_MISMATCHES=10
export CHIM_SEGMENT_MIN=1
export FEATURE_TYPE='exon'
export FILTER_INTRON_MOTIFS='None'
export HTSEQ_COUNT_MODE="union"

# command line arguments for the variant calling part
export ARG_BT_DEPTH=true
export FILE_EXONBED=$DIR_DATA/ensemble-genes/Homo_sapiens.GRCh37.75.onlyCDSCodonUTR.chr1-22XY.merged.bed
export FILE_TARGETAREAS=$DIR_DATA/ensemble-genes/Homo_sapiens.GRCh37.75.onlyCDSCodonUTR.chr1-22XY.plusmin5.merged.bed

# Restrict RNAseq variant calling to Ensembl exons
export ARG_TARGETAREAS="-L $DIR_DATA/ensemble-genes/Homo_sapiens.GRCh37.75.onlyCDSCodonUTR.chr1-22XY.plusmin5.merged.bed"

# allow gatk to work with rna-seq reads
export GATK_RNASEQ="-U ALLOW_N_CIGAR_READS"

export ARG_CLEANUP=true

