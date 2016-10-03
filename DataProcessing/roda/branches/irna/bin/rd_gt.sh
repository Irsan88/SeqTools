

function gendertest()
{
	TOTAL=$(cat $1 | awk '$1=="chrX"{ print $0 }' | awk '$7=="PASS"{ print $(NF) }' | awk -F ":" '{ print $1 }')
	#TOTAL=$(cat $1 | awk '$1=="chrX"{print $(NF) }' | awk -F ":" '{ print $1 }')
	#echo "$TOTAL"
	HOMOZYGOTE=$(echo "$TOTAL" | grep -c '1/1')
	HETEROZYGOTE=$(echo "$TOTAL" | grep -c '0/1')
	MULTIALLELIC=$(echo "$TOTAL" | grep -c '1/2')

	TOTPERC=$(echo "$TOTAL" | wc -l)
	HOMOZ_RATE=$(echo 'scale=2;'$HOMOZYGOTE'*100/'$TOTPERC | bc)
	HETEROZ_RATE=$(echo 'scale=2;'$HETEROZYGOTE'*100/'$TOTPERC | bc)
	RATIO=$(echo 'scale=3;'$HOMOZYGOTE'/'$HETEROZYGOTE | bc)

	echo RESULTS
	if [[ $(echo "if ("$RATIO" < 2) 1 else 0" | bc) -eq 1 ]]; then
		echo "Test subject is FEMALE"
	elif [[ $(echo "if ("$RATIO" > 3) 1 else 0" | bc) -eq 1 ]]; then
		echo "Test subject is MALE"
	else
		echo "Cannot confidently determine test subjects gender."
	fi
	echo
	echo -e "Homozygosity/heterozygosity ratio:\t"$RATIO
	echo -e "Homozygosity rate:\t"$HOMOZ_RATE" %"
	echo -e "Heterozygosity rate:\t"$HETEROZ_RATE" %"
	echo
	echo -e "Homozygote SNPs:\t"$HOMOZYGOTE
	echo -e "Heterozygote SNPs:\t"$HETEROZYGOTE
	echo -e "Multi-allelic SNPs:\t"$MULTIALLELIC
	echo
	echo -e "Error:\t"$(echo "$TOTAL" | grep -c -v -e "1/1" -e "0/1" -e "1/2")
}

gendertest $FILE_OUTPUT.snps.filt.vcf > $FILE_OUTPUT.stat.gender

# If there is a bed file available, load in here
FILE_BED=$FILE_TARGETAREAS
if [ ! $FILE_EXONBED = "" ]; then
	FILE_BED=$FILE_EXONBED
fi

# If bed file: perform calculations on the target. Also create XHMM input file
if [ ! $FILE_BED = "" ]; then	
	$SCRIPT_JAVA -Xmx3072m -jar $SCRIPT_GATK \
	-T DepthOfCoverage \
	-I $FILE_OUTPUT.recal.bam \
	-L $FILE_BED -R $FILE_REFERENCE \
	-U ALLOW_N_CIGAR_READS \
	-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase \
	--omitLocusTable --minBaseQuality 0 --minMappingQuality 20 \
	--start 1 --stop 5000 --nBins 200 --includeRefNSites -o $FILE_OUTPUT.stat.XHMM
fi
