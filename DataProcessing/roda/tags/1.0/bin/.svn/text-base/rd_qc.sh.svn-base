

ARG_TARGETAREAS_DEPTH=""
if [ ! $FILE_TARGETAREAS = "" ]; then
	ARG_TARGETAREAS_DEPTH="-b $FILE_TARGETAREAS"
fi

if $ARG_ST_DEPTH; then
	$SCRIPT_SAMTOOLS depth $ARG_TARGETAREAS_DEPTH -q 10 $FILE_OUTPUT.recal.bam > $FILE_OUTPUT.stat.recal.bam.depth
fi

if $ARG_BT_DEPTH; then
	$SCRIPT_BT_COVBED -abam $FILE_OUTPUT.recal.bam $ARG_TARGETAREAS_DEPTH -hist | awk '$1 == "all"' > $FILE_OUTPUT.stat.recal.bam.hist
fi

FILE_BLAAT=$FILE_TARGETAREAS
if [ ! $FILE_EXONBED = "" ]; then
	FILE_BLAAT=$FILE_EXONBED
fi

if [ ! $FILE_BLAAT = "" ]; then
	$SCRIPT_SAMTOOLS view -b -L $FILE_BLAAT $FILE_OUTPUT.recal.bam -o $FILE_OUTPUT.target.recal.bam
	$SCRIPT_SAMTOOLS flagstat $FILE_OUTPUT.target.recal.bam > $FILE_OUTPUT.stat.target.recal.bam.flagstat 
	rm $FILE_OUTPUT.target.recal.bam
fi
