#$ -pe singlenode 4

$SCRIPT_JAVA -Xmx4g -jar $SCRIPT_GATK \
	-T RealignerTargetCreator \
	-R $FILE_REFERENCE \
	-o $FILE_OUTPUT.mark.list \
	-I $FILE_OUTPUT.mark.bam \
	--known $FILE_KNOWNINDELS \
	$ARG_TARGETAREAS \
	-nt 4

# Does support parallelisation nt
# Specify nt 8, takes 235% cpu, down to 4
