
for SAMPLE in $1/*.conf*
do
	NAME_SAMPLE=$(basename "$SAMPLE")
	EXT_SAMPLE=${NAME_SAMPLE#*.}
	NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
	echo -e '\tProcessing:\t' $NAME_SAMPLE
	python ./new_qc_v6.py $1/$NAME_SAMPLE -general ~/Workspace/bedfiles/contis_general_Alessandra_130311.bed -exon ~/Workspace/bedfiles/contis_exon_Alessandra_130423.bed -panel ~/Workspace/bedfiles/130627_contis_genePanels\
		#> $1/$NAME_SAMPLE.statlog
	rm $1/*.png
	echo '' # Clear a line to keep things readable
done

