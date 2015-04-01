# this script assumes that there is a directory called ./bams 
# where (links of) all the bam files are
# It converts these bams to sams, sorts them for HTSeq input, 
# and HTseq does the read summarization

i=1
for f in $(find bams/ -name '*.bam')
do
	name=$(basename $f)
	echo "Sample $i:$name"
	samtools view -h $f > $name.sam
	./filterBamCigar.pl $name.sam > $name.splice
	java -jar /bioinformatics/programs/picard-tools-1.115/SortSam.jar INPUT=$name.splice OUTPUT=$name.splice.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT
	htseq-count \
		--mode=union \
		--stranded=no \
		--minaqual=35 \
		--type=exon \
		--idattr=gene_id \
		$name.splice.sam \
		/bioinformatics/library/gene-annotation/ensembl-genes-rnaseq/Homo_sapiens.GRCh37.75.chr1-22XY.gtf  > $name.htseq.ssv
	rm $name.sam $name.splice $name.splice.sam
	i=$((i+1))
done 
