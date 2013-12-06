# 1)Prepare a directory containing the FastQC-outputs
# of each fastq-file. Only include the FastQC-directory
# output, not the zip file.
#
# make sure htmldoc is installed
#
# [user@machine ~]$ bash fastqcToPdf.sh inputDir_containingFastQCs/
#
# Author: E.I. Kooi
# Date: 6 Dec 2013


inputDir=$1
outputDir='./FastQC_PDFs'

for f in $(ls $inputDir)
do
	htmldoc \
		--webpage \
		--browserwidth 800 \
		--fontsize 7 \
		-f $outputDir/$f.pdf \
		$inputDir/$f/fastqc_report.html	
done
