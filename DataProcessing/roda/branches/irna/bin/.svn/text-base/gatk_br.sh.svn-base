#$ -pe singlenode 4

$SCRIPT_JAVA -Xmx4g -jar $SCRIPT_GATK \
	-l INFO \
	-R $FILE_REFERENCE \
	-I $FILE_OUTPUT.rlgnfmi.bam \
	-T BaseRecalibrator \
	-knownSites $FILE_DBSNP \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov ContextCovariate \
	$ARG_TARGETAREAS \
	-nct 4 \
	-o $FILE_OUTPUT.recal.grp


# Does NOT support parallelisation nt, does support nct
# specify to use 8, really use about 300%... back to max 4 then
# Removed --disable_indel_quals
