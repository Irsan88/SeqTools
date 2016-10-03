
$SCRIPT_SAMTOOLS flagstat $FILE_OUTPUT.sort.bam > $FILE_OUTPUT.sort.bam.flagstat
echo 'In sort.bam mapped with MQ of 0: ' >> $FILE_OUTPUT.sort.bam.flagstat 
$SCRIPT_SAMTOOLS view -F 4 $FILE_OUTPUT.sort.bam | awk '$5=="0" ' | wc -l >> $FILE_OUTPUT.sort.bam.flagstat
echo 'Number of reads left after RETRO filter: ' >> $FILE_OUTPUT.sort.bam.flagstat
python $SCRIPT_CALCTOWERS $FILE_OUTPUT.pickle >>$FILE_OUTPUT.sort.bam.flagstat
