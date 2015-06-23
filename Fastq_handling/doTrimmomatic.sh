#!/bin/bash

# Usage: ./doTrimmomatic.sh /path/to/input/directory /path/to/output/directory trimmomatic.jar adapters.fa

# prerequisites: 
# - make sure all R1 and R2 files of each sample are in provided input dir
# - the _ character is used to separte sample name from the file name prefix
#   e.g. VU834-10M_ATCGT_R1.fastq.gz will have sample name VU834-10M
# - make sure provided output directory exists

# Author: E.I. Kooi
# Date: 10 July 2014

# To-Do:
# - check if trimmomatic is installed
# - check command line arguments

# Fetch command line arguments
inputDir=$1
outputDir=$2
trimJar=$3
adapters=$4

# loop over all samples in the provided input dir 
# and perform trimmomatic cleaning
for sample in $(ls $inputDir | cut -f1 -d'_' | sort -u)
do
	# fetch the names for the forward and reverse fastq-files
	forward=$(readlink -f $inputDir/$sample\_*R1*.gz)
	reverse=$(readlink -f $inputDir/$sample\_*R2*.gz)
	java -classpath $trimJar org.usadellab.trimmomatic.TrimmomaticPE \
		-threads 8 \
		-phred33 \
		$forward \
		$reverse \
		$outputDir/$sample\_R1_paired_clean.fastq \
		$outputDir/$sample\_R1_notpaired_clean.fastq \
		$outputDir/$sample\_R2_paired_clean.fastq \
		$outputDir/$sample\_R2_notpaired_clean.fastq \
		ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
	
	# clean up singleton files
	rm $outputDir/$sample\_R1_notpaired_clean.fastq 
	rm $outputDir/$sample\_R2_notpaired_clean.fastq
	
	# zip fastq files
	gzip $outputDir/$sample\_R1_paired_clean.fastq
	gzip $outputDir/$sample\_R2_paired_clean.fastq
done
