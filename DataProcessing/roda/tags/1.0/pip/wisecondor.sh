INPUT=$NAME_SAMPLE
INPUT_R1=$FILE_INPUT #$DIR_SAMPLE/*R1*.fastq*

STEP=bwa_aln
qsub \
	-N ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${STEP} \
	-V \
	-e $DIR_SGELOG/${STEP}1.e${STAMP} \
	-o $DIR_SGELOG/${STEP}1.o${STAMP} \
	$DIR_BASE/bin/$STEP.sh \
		$INPUT_R1 \
		${FILE_OUTPUT}
PREVSTEP=$STEP

STEP=bwa_sse
qsub -hold_jid ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${PREVSTEP} \
	-V \
	-N ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${STEP} \
	-e $DIR_SGELOG/${STEP}.e${STAMP} \
	-o $DIR_SGELOG/${STEP}.o${STAMP} \
	$DIR_BASE/bin/$STEP.sh \
		$INPUT_R1
PREVSTEP=$STEP

PIPELINE=('pt_ss' 'wc_rf' 'wc_cs');

runpipe $PIPELINE
