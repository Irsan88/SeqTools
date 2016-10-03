#$ -pe singlenode 2

# Obtain the other end fastq by replacing R1 in the filename
INPUT=$NAME_SAMPLE
INPUT_R1=$FILE_INPUT #$DIR_SAMPLE/*R1*.fastq*
INPUT_R2=${INPUT_R1/_R1_/_R2_} # Replace R1 attempt one: _R1_L001.fastq
INPUT_R2=${INPUT_R2/_R1./_R2.} # Replace R1 attempt two: _R1.fastq
# check if R2 file exists, if not, assume it is single-end data
if [ ! -f $INPUT_R2 ]; then
    INPUT_R2=""
fi
# INPUT_R2="" # is a hack to force single-end

echo "Fastq-files used for $NAME_SAMPLE: $INPUT_R1 and $INPUT_R2"


$SCRIPT_STAR \
	--genomeDir $INDEX_STAR_HG19 \
	--readFilesIn $INPUT_R1 $INPUT_R2 \
	--runThreadN 8 \
	--runMode alignReads \
	--readFilesCommand zcat \
	--sjdbScore 2 \
	--outFilterMismatchNmax $AMOUNT_MISMATCHES \
	--chimSegmentMin $CHIM_SEGMENT_MIN \
	--outFileNamePrefix $FILE_OUTPUT.


	

