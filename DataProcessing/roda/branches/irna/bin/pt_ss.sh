#$ -pe singlenode 2

$SCRIPT_JAVA -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_SORTSAM \
	SO=coordinate \
	INPUT=$FILE_OUTPUT.reordered.bam \
	OUTPUT=$FILE_OUTPUT.reordered.sorted.bam \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$DIR_OUTPUT \
	CREATE_INDEX=true

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.reordered.ba*
fi