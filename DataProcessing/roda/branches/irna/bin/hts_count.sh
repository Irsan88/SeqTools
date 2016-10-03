#$ -pe singlenode 2

$SCRIPT_HTSEQ \
	--mode=$HTSEQ_COUNT_MODE \
	--stranded=$STRAND_DIRECTION \
	--minaqual=35 \
	--type=$FEATURE_TYPE \
	--idattr=gene_id \
	$FILE_OUTPUT\.sam \
	$GTF > $FILE_OUTPUT\.htSeqCounts_ensemblID.ssv

