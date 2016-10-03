#$ -pe singlenode 2

$SCRIPT_JAVA -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_GATK \
	-I $FILE_OUTPUT.mark.bam \
	-R $FILE_REFERENCE \
	-T IndelRealigner \
	-targetIntervals $FILE_OUTPUT.mark.list \
	-known $FILE_KNOWNINDELS \
	-L $FILE_EXONBED \
	$GATK_RNASEQ \
	-o $FILE_OUTPUT.rlgn.bam

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.mark.ba*
	rm $FILE_OUTPUT.mark.list
	#rm $FILE_OUTPUT.mark.metrics
fi

# Does NOT support parallelisation except for Scatter-Gather using the GATK queue implementation
