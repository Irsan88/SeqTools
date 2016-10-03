#!/bin/bash

input=$1

# remove header
sed '1d' $input | \
# first get rid of the double information
grep -v 'CN values' | \
# column 1 is GISTIC peak id
# column 3 is regions limits chrX:start-end(probe1-probeN)
# column 6 is q-value, 7 is residual q-value
cut -f 1,3,7 | \
# reformat
sed 's/\s\+Peak\s\+/Peak/g' | \
sed 's:(.*)::g' | \
sed 's/e-/keepthis/g' | \
sed 's/[:-]/\t/g' | \
sed 's/keepthis/e-/g' | \
awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$1,$5}'
