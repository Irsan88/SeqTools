inFile=$1
# remove the CN value matrix info
grep -v 'CN values' $inFile | \
# take wide region limits, peak ID, cytoband and residual q-value
awk  'BEGIN{FS="\t"}{print $3,$1,$2,$7}' | \
# remove the probes (number:number)
sed 's/(probes [0-9]*:[0-9]*)//g' | \
# change double white space characters to tab 
sed 's/\s\+/\t/g' | \
# remove header
sed '1d' | \
# ad Amplification and Peak in one column
sed 's/\tPeak\t/Peak/g' | \
# split chr:start-end in columns, but preserve 1.45e-5
sed 's/[:]/\t/g' | \
# make sure numbers like 1.45e-5 are preserved
sed 's/e-/e@/g' | sed 's/-/\t/g' | sed 's/e@/e-/g' | \
# merge AmplificationPeak23	1q32.2
awk 'BEGIN{FS="\t";OFS=""}{print $1,"\t",$2,"\t",$3,"\t",$4,"-",$5,"\t",$6}' | \
# rename AmplificationPeak and DeletionPeak
sed 's/AmplificationPeak/Amp/g' | \
sed 's/DeletionPeak/Del/g'

