
for SAMPLE in $1/*.gender
do
	echo $SAMPLE
	NAME_SAMPLE=$(basename "$SAMPLE")
	echo $NAME_SAMPLE
	EXT_SAMPLE=${NAME_SAMPLE#*.}
	NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
	echo -e '\tAdding:\t' $NAME_SAMPLE
	#python ./new_qc.py $1/$NAME_SAMPLE \
	#	> $1/$NAME_SAMPLE.statlog
	#echo $1/$NAME_SAMPLE.statlog
	#echo ./new_qc_v6.py $1/$NAME_SAMPLE -generalbed ~/Workspace/bedfiles/SeqCap_EZ_Exome_v3_capture_450bpExtra_noOverlap_130307.bed -exon ~/Workspace/bedfiles/exome_SeqCap_EZ_v3_capture_minWhiteLines_130307.bed -hist
	python ./new_qc_v6.py $1/$NAME_SAMPLE -generalbed ~/Workspace/bedfiles/SeqCap_EZ_Exome_v3_capture_450bpExtra_130716_sortmerge.bed -exon ~/Workspace/bedfiles/SeqCap_EZ_Exome_v3_capture_minWhiteLines_sorted.bed -hist > $1/$NAME_SAMPLE.statlog
	rm $1/*.png
	echo '' # Clear a line to keep things readable
done

