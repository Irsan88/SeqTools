#!/bin/bash

chim=$1
genesBed=$2

# also see https://groups.google.com/forum/#!msg/rna-star/HUxFCaHSX6c/iSudPgceUXkJ
# for explanation of columns of $chim file
awk 'BEGIN{FS="\t";OFS="\t"}{if($1!="chrM" && $4!="chrM" && $7>0 && $8+$9<=5)print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | \
sort | \
uniq -c | \
sort -k1,1rn \
sed 's/^\s\+//g' > $chim\.unique
 

