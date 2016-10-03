
# Add argument to numerous GATK Steps to only look into regions we are interested in
export ARG_TARGETAREAS=""
if $ARG_GATK_BED; then
#if [ ! $FILE_TARGETAREAS = "" ]; then
	export ARG_TARGETAREAS="-L $FILE_TARGETAREAS"
fi
export ARG_HCAREAS=""
if $ARG_BT_DEPTH; then
	export ARG_HCAREAS="-L $FILE_TARGETAREAS"
fi

# Obtain the other end fastq by replacing R1 in the filename
INPUT=$NAME_SAMPLE
INPUT_R1=$FILE_INPUT #$DIR_SAMPLE/*R1*.fastq*
INPUT_R2=${INPUT_R1/_R1_/_R2_} # Replace R1 attempt one: _R1_L001.fastq
INPUT_R2=${INPUT_R2/_R1./_R2.} # Replace R1 attempt two: _R1.fastq

# Align both ends
STEP=bwa_aln
#echo '
qsub \
	-N ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${STEP}1 \
	-V \
	-e $DIR_SGELOG/${STEP}1.e${STAMP} \
	-o $DIR_SGELOG/${STEP}1.o${STAMP} \
	$DIR_BASE/bin/$STEP.sh \
		$INPUT_R1 \
		${FILE_OUTPUT}_R1

qsub \
	-N ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${STEP}2 \
	-V \
	-e $DIR_SGELOG/${STEP}2.e${STAMP} \
	-o $DIR_SGELOG/${STEP}2.o${STAMP} \
	$DIR_BASE/bin/$STEP.sh \
		$INPUT_R2 \
		${FILE_OUTPUT}_R2
PREVSTEP=$STEP
#'

# Wait for ALN1 and ALN2 to finish,
# Combine alignments
# CAUTION: special case here, more arguments passed into echo qsub script than usual
STEP=bwa_spe
qsub -hold_jid ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${PREVSTEP}1,${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${PREVSTEP}2 \
	-V \
	-N ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${STEP} \
	-e $DIR_SGELOG/${STEP}.e${STAMP} \
	-o $DIR_SGELOG/${STEP}.o${STAMP} \
	$DIR_BASE/bin/$STEP.sh \
		$INPUT_R1 \
		$INPUT_R2
#'
PREVSTEP=$STEP

PIPELINE=('pt_ss_exome');

runpipe $PIPELINE
