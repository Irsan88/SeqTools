#$ -pe singlenode 2

$SCRIPT_JAVA -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_MARKDUP \
	INPUT=$FILE_OUTPUT.sort.bam \
	OUTPUT=$FILE_OUTPUT.mark.bam \
	METRICS_FILE=$FILE_OUTPUT.mark.metrics \
	CREATE_INDEX=true \
	TMP_DIR=$DIR_OUTPUT \
	VALIDATION_STRINGENCY=LENIENT

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.sort.bam
	rm $FILE_OUTPUT.sort.bai
	rm $FILE_OUTPUT.mark.metrics # Did not provide the right information (dafuq...)
fi

# Generate some info for the statistics later on
$SCRIPT_SAMTOOLS flagstat $FILE_OUTPUT.mark.bam > $FILE_OUTPUT.stat.mark.bam.flagstat
echo 'In mark.bam mapped with MQ of 0: ' >> $FILE_OUTPUT.stat.mark.bam.flagstat
$SCRIPT_SAMTOOLS view -F 4 $FILE_OUTPUT.mark.bam | awk '$5=="0" ' | wc -l >> $FILE_OUTPUT.stat.mark.bam.flagstat

