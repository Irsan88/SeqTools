#$ -pe singlenode 2

$SCRIPT_JAVA -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_SORTSAM \
	SO=coordinate \
	INPUT=$FILE_OUTPUT.rg.sam \
	OUTPUT=$FILE_OUTPUT.rg.sort.sam \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$DIR_OUTPUT \
	CREATE_INDEX=true

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.rg.sa*
fi
