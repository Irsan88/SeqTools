

$SCRIPT_BWA sampe \
	-r "@RG\tID:idname\tLB:libname\tSM:$NAME_SAMPLE\tPL:ILLUMINA" \
	$FILE_REFERENCE \
	${FILE_OUTPUT}_R1.sai \
	${FILE_OUTPUT}_R2.sai \
	$1 \
	$2 \
	> $FILE_OUTPUT.sam

if $ARG_CLEANUP; then
	rm ${FILE_OUTPUT}_R1.sai
	rm ${FILE_OUTPUT}_R2.sai
fi
