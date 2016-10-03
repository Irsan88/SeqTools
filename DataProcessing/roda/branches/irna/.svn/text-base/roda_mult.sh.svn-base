set -e

# IO
DIR_INPUT=$1
DIR_OUTPUT=$2

if [ ! -d $DIR_OUTPUT ]
then
	mkdir $DIR_OUTPUT
fi

DIR_BASE=$(dirname $0) # Dir script resides in
DIR_CUR=${PWD} # Dir script is called from



echo -e "Files found to queue for processing:"
#ls $DIR_INPUT/*/*.fastq*
FILES=$( find $1 -name '*_R1*fastq*' )

#for SAMPLE in $FILES
#do
#	echo $SAMPLE
#done
COUNT=0
for SAMPLE in $FILES 
do
	#echo '' # Clear a line to keep things readable
	NAME_SAMPLE=$(basename "$SAMPLE")
	EXT_SAMPLE=${NAME_SAMPLE#*.}
	NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
	NAME_SAMPLE=$(echo $NAME_SAMPLE | cut -d _ -f 1)
	echo -e '\e[0m\tName:\t' $NAME_SAMPLE
	echo -e '\e[0;36m\tPath:\t' $SAMPLE
	COUNT=$(($COUNT+1))
done
#COUNT=$( echo $FILES )
echo -e "\e[0mTotalling\e[1;31m $COUNT \e[0msamples"
read -p "Are you sure this is correct? [Y/n] " -n 1 -r
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
	echo -e '\n\e[1;31mQueueing cancelled\e[0m'
	exit 1
fi

echo -e '\nQueueing files:' # Clear a line to keep things readable

for SAMPLE in $FILES
do
	NAME_SAMPLE=$(basename "$SAMPLE")
	EXT_SAMPLE=${NAME_SAMPLE#*.}
	NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
	NAME_SAMPLE=$(echo $NAME_SAMPLE | cut -d _ -f 1)
	#echo -e '\e[0mQueueing:' 
	#echo -e '\tPath:\t' $SAMPLE
	#echo -e '\tName:\t' $NAME_SAMPLE '\e[0;36m'
	echo -e '\e[0m\tAdding:\t' $NAME_SAMPLE '\e[0;36m'

	$DIR_BASE/roda_single.sh $SAMPLE $DIR_OUTPUT/$NAME_SAMPLE $3
	echo '' # Clear a line to keep things readable
	#sleep 1 # Hold on a second to ensure unique timestamps
done

echo -e '\e[0mFinished queueing all samples'
