# Goal: LOH-analysis
# Arguments: 
#	1) bam-file 
#	2) reference genome
#	3) bed-file with snps 
#	4) threshold for the coverage that a base 
#	   position should have to be considered 
# Output: one column text file with B allele frequency for each SNP position
# Author: E.I. Kooi
# Date: 24 Nov 2013
# Dependencies: samtools

# Usage: bam2baf.sh file.bam hg19.fa fileWithSNPPositions.bed 10 > BAF.bed

bam=$1
ref=$2
snps=$3
cov=$4

# get directory where this file is stored in so 
# you can use mpileup2baf.pl
dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

samtools mpileup \
	-f $ref \
	-l $snps \
	-q 35 \
	-Q 20 \
	$bam | \
$dir/mpileup2baf.pl \
	--min-reads $cov 
