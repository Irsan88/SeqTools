
$SCRIPT_JAVA -Xmx4g -jar $SCRIPT_GATK \
	-T PrintReads \
	-I $FILE_OUTPUT.reordered.sorted.bam \
	-o $FILE_OUTPUT.reordered.sorted.mapQ60.bam \
	-R $FILE_REFERENCE \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	$GATK_RNASEQ

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.reordered.sorted.ba*
fi
