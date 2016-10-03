#$ -pe singlenode 4

$SCRIPT_CUFFLINKS \
	--output-dir $DIR_OUTPUT/cufflinks \
	--GTF-guide $GTF \
	--num-threads 4 \
	--multi-read-correct \
	--library-type $CUFF_LIBRARY_TYPE \
	--mask-file $CUFF_rRNA_MASK \
	--frag-bias-correct $HG19_FA \
	$FILE_OUTPUT.rg.sort.sam


