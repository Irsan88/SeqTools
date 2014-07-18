# bedops should be installed
# input should be in DNAcopy-format
input=$1
TMP_DIR="./tmp"
if [[ -e $TMP_DIR ]]; then
    rm -r $TMP_DIR
fi
mkdir $TMP_DIR

# get all positions of segments 
sed '1d' $input | \
cut -f2-4 | \
sort-bed - > $TMP_DIR/segment-regions

# get partitioned regions of all segments
bedops --partition $TMP_DIR/segment-regions > $TMP_DIR/reduced-regions

# get sample names 
sed '1d' $input | \
cut -f1 | \
sort | \
uniq > $TMP_DIR/samples

# convert each sample to bed file
for sample in $(cat $TMP_DIR/samples)
do
	awk -v s="$sample" 'BEGIN{FS="\t";OFS="\t"}{if(NR>1 && $1==s){print $2,$3,$4,$1,$6}}' $input | \
	sort-bed - | \
 	bedmap --mean $TMP_DIR/reduced-regions - | sed 's:NAN:NA:g' > $TMP_DIR/$sample.mappedBack
done

# collect all data in matrix
# header line
echo -e "chr\tstart\tend" > $TMP_DIR/h1
cat $TMP_DIR/samples | tr '\n' '\t' | sed 's:\t$::g' > $TMP_DIR/h2
paste $TMP_DIR/h1 $TMP_DIR/h2 > $TMP_DIR/h3
paste $TMP_DIR/reduced-regions $TMP_DIR/*.mappedBack > reduced-matrix.tsv
cat $TMP_DIR/h3 reduced-matrix.tsv > tmp-file && mv tmp-file reduced-matrix.tsv

# clean up
rm -r ./tmp
