#$ -pe singlenode 4

$SCRIPT_JAVA -Xmx16g -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_MARKDUP \
	INPUT=$FILE_OUTPUT.rg.sort.sam \
	OUTPUT=$FILE_OUTPUT.nodupl.sam \
	REMOVE_DUPLICATES=true \
	METRICS_FILE=$FILE_OUTPUT.mark.metrics \
	CREATE_INDEX=false \
	TMP_DIR=$DIR_OUTPUT \
	VALIDATION_STRINGENCY=LENIENT

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.rg.sort.sa*
	rm $FILE_OUTPUT.mark.metrics # Did not provide the right information (dafuq...)
fi

