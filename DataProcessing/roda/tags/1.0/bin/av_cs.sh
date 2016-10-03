# Make the file human readable by splitting the final column in several columns
awk -F'[:\t]' '{for(i=1;i<=NF;i++)if(i<=7||i>=14)printf$i"\t";printf$8"\t";print""}' \
	$FILE_OUTPUT.snps.filt.vcf \
		> $FILE_OUTPUT.snps.2av.vcf


$SCRIPT_AV_CONVERT --format vcf4 \
	--includeinfo \
	$FILE_OUTPUT.snps.2av.vcf | \
	awk '{for(i=1;i<=NF;i++)if(i<=5||i>=11)printf$i"\t";print""}' \
		> $FILE_OUTPUT.snps.av
#		awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 "\t" $18 "\t" $20 "\t" $16}' \


$SCRIPT_AV_SUMMARIZE --buildver hg19 \
	--ver1000g 1000g2012feb \
	--verdbsnp $ARG_AV_DBSNP \
	$FILE_OUTPUT.snps.av \
	$FILE_ANNOVARDB \
	-outfile $FILE_OUTPUT.snps.av

# Fix Annovars header output
sed -i 's/Otherinfo/Quality,Filter,Genotype,DepthPerAllele,DepthOfCoverage,GenotypeQuality,PhredLikelyhood,Info/' \
	$FILE_OUTPUT.snps.av.exome_summary.csv
sed -i 's/Otherinfo/Quality,Filter,Genotype,DepthPerAllele,DepthOfCoverage,GenotypeQuality,PhredLikelyhood,Info/' \
	$FILE_OUTPUT.snps.av.genome_summary.csv
#sed -i 's/;//g' $FILE_OUTPUT.snps.av.exome_summary.csv
#sed -i 's/;//g' $FILE_OUTPUT.snps.av.genome_summary.csv

if $ARG_CLEANUP; then
	rm $FILE_OUTPUT.snps.2av.vcf
	rm $FILE_OUTPUT.snps.av
	rm $FILE_OUTPUT.snps.av.hg*
	rm $FILE_OUTPUT.snps.av.*variant_function
	rm $FILE_OUTPUT.snps.av.log
	rm $DIR_OUTPUT/snappy*.so
fi
