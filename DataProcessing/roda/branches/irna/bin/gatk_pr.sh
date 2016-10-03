#$ -pe singlenode 4

$SCRIPT_JAVA \
	-Djava.io.tmpdir=$DIR_OUTPUT \
	-Xmx4g -jar $SCRIPT_GATK \
	-T PrintReads \
	-R $FILE_REFERENCE \
	-I $FILE_OUTPUT.rlgnfmi.bam \
	-BQSR $FILE_OUTPUT.recal.grp \
	$ARG_TARGETAREAS \
	$GATK_RNASEQ \
	-nct 4 \
	-o $FILE_OUTPUT.recal.bam

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.rlgnfmi.ba*
	rm $FILE_OUTPUT.recal.grp
fi

# Does NOT support parallelisation nt, does support nct
# specify 8, uses 300%, back to 4
