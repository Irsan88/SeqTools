#$ -pe singlenode 2

$SCRIPT_JAVA -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_SORTSAM \
	SO=queryname \
	INPUT=$FILE_OUTPUT.rg.sam \
	OUTPUT=$FILE_OUTPUT.sam \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$DIR_OUTPUT \
	CREATE_INDEX=false

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.rg.sa*
fi
