# -pe singlenode 3

$SCRIPT_SAMTOOLS rmdup -s $FILE_OUTPUT.sort.bam - | \
	$SCRIPT_SAMTOOLS view - -q 1 | \
	$SCRIPT_PYTHON $SCRIPT_WC_CONSAM \
		$FILE_OUTPUT.pickle

