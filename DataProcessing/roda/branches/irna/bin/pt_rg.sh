#$ -pe singlenode 2

$SCRIPT_JAVA -Djava.io.tmpdir=/tmp \
	-jar $SCRIPT_PT_ADDREADGROUP \
	INPUT=$FILE_OUTPUT.Aligned.out.sam \
	OUTPUT=$FILE_OUTPUT.rg.sam \
	RGLB=RNAseq \
	RGPL=illumina \
	RGPU=platform1 \
	RGSM=$NAME_SAMPLE \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$DIR_OUTPUT \
	CREATE_INDEX=true

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.Aligned.out.sa*
fi
