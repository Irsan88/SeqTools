#!/bin/bash

# Usage: ...
# Arguments:	1) bam-file (compatible with bam2bed)
#		2) windows where counts are needed in bed-file
#		3) mapping quality threshold
# Output: one column text file with counts for each window
# Author: E.I. Kooi
# Date: 24 Nov 2013
# Dependencies: bedops suite

# usage: bam2counts.sh bamfile.bam windows.bed 35 > answer.count

bam=$1
windows=$2
mapQ=$3

# fetch the sample name from the bamfile
# this will be the prefix of the output-file
name=$(basename $bam)
name=${name/\.bam/}

# convert bam to sam for parsing
samtools view $bam | \
# only keep reads that:
# -are not too close to chromosome position 0
# -have a mapping quality bigger than provided mapQ
awk -v q="$mapQ" 'BEGIN{FS="\t"}{if($4 > 200 && $5 > q)print $3,$4,$4+100}' | \
# sort resulting bed-file
sort-bed - | \
# count the amount of mapped reads
# overlapping with the features in the
# provided windows bed-file
bedmap --count $windows -



