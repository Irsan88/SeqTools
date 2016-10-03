#! /bin/bash 

#run as ./runSamtoolsMpileup.sh [input regions file] [.bam file] [output.txt file]

#/illumina/diagnostics/bin/samtools-0.1.18/samtools index $2
cat $1 | \
while read line
do
	region="$(echo $line | awk '{print $1":"$2"-"$2}')"
    /illumina/diagnostics/bin/samtools-0.1.18/samtools mpileup -f /illumina/diagnostics/lib/hg19daud/ucsc.daud.hg19.fasta -q 1 -Q 13 -r $region $2 >> $3 2> /dev/null #| awk ' $4 >= 30 {print $0}' >> $3 
	#echo "$region done"
done
