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

samtools mpileup \
	-f $ref \
	-l $snps \
	-q 35 \
	-Q 20 \
	$bam | \
mpileup2baf.pl \
	--min-reads $cov 
