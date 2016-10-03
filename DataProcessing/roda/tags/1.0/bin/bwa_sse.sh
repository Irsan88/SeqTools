

$SCRIPT_BWA samse \
	-r "@RG\tID:idname\tLB:libname\tSM:$NAME_SAMPLE\tPL:ILLUMINA" \
	$FILE_REFERENCE \
	${FILE_OUTPUT}.sai \
	-n -1 \
	$1 \
	> $FILE_OUTPUT.sam

if $ARG_CLEANUP; then
	rm ${FILE_OUTPUT}.sai
fi
