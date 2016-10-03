

$SCRIPT_JAVA -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_FIXMATE \
	INPUT=$FILE_OUTPUT.rlgn.bam \
	OUTPUT=$FILE_OUTPUT.rlgnfmi.bam \
	SORT_ORDER=coordinate \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$DIR_OUTPUT \
	CREATE_INDEX=true

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.rlgn.bam
	rm $FILE_OUTPUT.rlgn.bai
fi
