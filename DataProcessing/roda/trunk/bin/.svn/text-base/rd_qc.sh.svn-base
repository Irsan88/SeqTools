#$ -pe singlenode 2

# Set general bed (bed + extra region)
ARG_TARGETAREAS_DEPTH=""
if [ ! $FILE_TARGETAREAS = "" ]; then
	ARG_TARGETAREAS_DEPTH="$FILE_TARGETAREAS"
fi

# Calculate depth file for targeted kits
if $ARG_ST_DEPTH; then
	$SCRIPT_SAMTOOLS depth -b $ARG_TARGETAREAS_DEPTH -q 10 $FILE_OUTPUT.recal.bam > $FILE_OUTPUT.stat.recal.bam.depth
fi

# Calculate depth distribution for exome kit - also calculate horizontal coverage
if $ARG_BT_DEPTH; then
	$SCRIPT_SAMTOOLS mpileup -l $FILE_EXONBED -q 1 $FILE_OUTPUT.recal.bam | awk '{if($4 >= 1) {cnt1++}; if($4 >= 10) {cnt10++}; if($4 >= 30) {cnt30++}; if($4 >= 50) {cnt50++}}; END {print cnt1; print cnt10; print cnt30; print cnt50}' > $FILE_OUTPUT.stat.recal.bam.hist
	echo "64190747" >> $FILE_OUTPUT.stat.recal.bam.hist # Nr of bases in bed region exome, can be calculated with commented line right below this one
	#$SCRIPT_BT_COVBED -abam $FILE_OUTPUT.recal.bam -b $FILE_EXONBED -d | wc -l >> $FILE_OUTPUT.stat.recal.bam.hist
	$SCRIPT_BT_COVBED -abam $FILE_OUTPUT.recal.bam -b $FILE_EXONBED -hist | awk '$1 == "all"' >> $FILE_OUTPUT.stat.recal.bam.hist
fi

# If there is a bed file available, load in here
FILE_BLAAT=$FILE_TARGETAREAS
if [ ! $FILE_EXONBED = "" ]; then
	FILE_BLAAT=$FILE_EXONBED
fi

# If bed file: perform calculations on the target.
if [ ! $FILE_BLAAT = "" ]; then
	$SCRIPT_SAMTOOLS view -b -L $FILE_BLAAT -q 1 $FILE_OUTPUT.recal.bam -o $FILE_OUTPUT.target.recal.bam
	$SCRIPT_SAMTOOLS flagstat $FILE_OUTPUT.target.recal.bam > $FILE_OUTPUT.stat.target.recal.bam.flagstat
	echo 'Mapped reads in pair mapped on RS: ' $($SCRIPT_SAMTOOLS view $FILE_OUTPUT.target.recal.bam -f 0x0002 -f 0x0010 -F 0x0400 -c) >> $FILE_OUTPUT.stat.target.recal.bam.flagstat 
	echo 'Mapped reads in pair mapped on FS: ' $($SCRIPT_SAMTOOLS view $FILE_OUTPUT.target.recal.bam -f 0x0002 -F 0x0010 -F 0x0400 -c) >> $FILE_OUTPUT.stat.target.recal.bam.flagstat
	rm $FILE_OUTPUT.target.recal.bam
fi
