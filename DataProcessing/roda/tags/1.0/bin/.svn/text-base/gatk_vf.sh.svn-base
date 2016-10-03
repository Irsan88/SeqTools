

$SCRIPT_JAVA -Xmx4g -jar $SCRIPT_GATK \
	-R $FILE_REFERENCE \
	-T VariantFiltration \
	-V $FILE_OUTPUT.snps.vcf \
	-o $FILE_OUTPUT.snps.filt.vcf \
	$ARG_TARGETAREAS \
	--clusterWindowSize 10 \
	--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
	--filterName "HARD_TO_VALIDATE" \
	--filterExpression "DP < 5 " \
	--filterName "LowCoverage" \
	--filterExpression "QUAL < 30.0 " \
	--filterName "VeryLowQual" \
	--filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
	--filterName "LowQual" \
	--filterExpression "QD < 1.5 " \
	--filterName "LowQD" \
	--filterExpression "SB > -10.0 " \
	--filterName "StrandBias"

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.snps.vcf
	rm $FILE_OUTPUT.snps.vcf.idx
fi

# SB -> FS en 200 als value
