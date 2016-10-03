#$ -pe singlenode 8

$SCRIPT_JAVA -Xmx4g -jar $SCRIPT_GATK \
	-glm BOTH \
	-R $FILE_REFERENCE \
	-T UnifiedGenotyper \
	-I $FILE_OUTPUT.recal.bam \
	-D $FILE_DBSNP \
	-o $FILE_OUTPUT.snps.vcf \
	-metrics $FILE_OUTPUT.stat.snps.metrics \
	-stand_call_conf 50.0 \
	-stand_emit_conf 10.0 \
	-dcov 1000 \
	-A DepthOfCoverage \
	-A AlleleBalance \
	$ARG_TARGETAREAS \
	-nt 1 \
	-nct 8 #$ARG_GATK_THREADS

#if $ARG_CLEANUP; then
#	rm $FILE_OUTPUT.snps.metrics
#fi

# Does support parallelisation
# Lower paralellisation: Too many open files error otherwise. (dafuq...)
