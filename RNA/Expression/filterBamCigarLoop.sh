for f in $(find bams/ -name '*.bam')
do
	name=$(basename $f)
	echo $name
	samtools view -h $f > $name.sam
	./filterBamCigar.pl $name.sam > $name.splice
	java -jar /illumina/diagnostics/bin/picard-tools-1.81/SortSam.jar INPUT=$name.splice OUTPUT=$name.splice.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT
	htseq-count \
		--mode=union \
		--stranded=no \
		--minaqual=35 \
		--type=exon \
		--idattr=gene_id \
		$name.splice.sam \
		/illumina/diagnostics/lib/Homo_sapiens.GRCh37.65_v2.gtf  > $name.htseq.ssv
	rm $name.sam $name.splice
	samtools view -Sbh $name.splice.sam > /media/irsan/thromboseq/splicedOnly/$name.splice.bam
	rm $name.splice.sam
done 
