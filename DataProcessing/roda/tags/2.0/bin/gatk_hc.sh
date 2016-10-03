
$SCRIPT_JAVA -Xmx4g -jar $SCRIPT_GATK \
	-R $FILE_REFERENCE \
	-T HaplotypeCaller \
	-I $FILE_OUTPUT.recal.bam \
	-D $FILE_DBSNP \
	-o $FILE_OUTPUT.snps.vcf \
	-stand_call_conf 50.0 \
	-stand_emit_conf 10.0 \
	$ARG_TARGETAREAS \
	$ARG_HCAREAS

#if $ARG_CLEANUP; then
#	rm $FILE_OUTPUT.snps.metrics
#fi

# UG Commands:
#	-glm BOTH \
#	-dcov 1000 \
#	-A DepthOfCoverage \
#	-A AlleleBalance \
#	-nt 1 \
#	-nct 8 #$ARG_GATK_THREADS
