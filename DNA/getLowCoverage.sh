#!/bin/bash

###################################################
#
# INITIALIZATION
#
###################################################

bam=$1
bed=$2
t=$3

# check if usage is correct
# usage: getLowCoverage.sh alignments.bam exons.bed refGene_hg19.gtf 6

###################################################
#
# MAIN PROGRAM EXECUTION
#
###################################################

# fetch the sample name from the bamfile
# this will be the prefix of the output-file
name=$(basename $bam)
name=${name/\.bam/}
# first melt the bed-file
bedtools makewindows -b $bed -w 1 | sort-bed - > melted.bed
# then use bedtools to get coverage per base
bedtools coverage -abam $bam -b melted.bed | \
# output contains: chr, start, end, coverage, ..., ..., ...
awk -F "\t" -v T="$t" '{if($4 < T)print $1,$2,$3}' | \
# then collapse adjacent bases
bedops --merge - > lowCoverageRegions.bed
# annotate the resulting bed
bedmap --delim "\t" --echo --echo-map-id-uniq \
lowCoverageRegions.bed \
/home/irsan/references/beds/RefSeq_geneSymbols.bed > $name\_lowCoverageRegions_annotated.bed


###################################################
#
# FINISHING
#
###################################################
# clean up
rm melted.bed
rm lowCoverageRegions.bed
