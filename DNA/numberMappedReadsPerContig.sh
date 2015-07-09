#!/bin/bash

# check dependencies (samtools)
# check command line arguments
# provide help message

# author: E.I. Kooi
# date: 9 July 2015
# rationale: mapped cancer exomeseq data to viral genomes
# to see if there are viral traces in tumors. Now I want to
# see how much viral traces there are and which virals
# Each virus is one record in the reference genome fasta
# so 1 virus =~ 1 contig =~ 1 chromosome

bamInput=$1
mappingQualityThreshold=$2

# -F is negative filter, 4 means read is unmapped
# so -F 4 means dont return unmapped reads
samtools view -F 4 $bamInput | \
# filter reads with mapping quality (5th column) above (>) threshold
awk -v awkMapQ="$mappingQualityThreshold" 'BEGIN{FS="\t";OFS="\t"}{if($5 > awkMapQ)print $3}' - | \
# how much reads per contig
sort | uniq -c | \
# remove leading white spaces
sed 's/^\s\+//g' | \
# make tab-delimited (is tricky: contig names should not contain whitespaces)
# I believe that bwa already removes white spaces from contig names
# because the initial viral contig names contained white spaces but after
# mapping they didn't
sed 's/\s/\t/g' | \
# sort by number, large to small
sort -k1,1 -n -r # print to stdout

