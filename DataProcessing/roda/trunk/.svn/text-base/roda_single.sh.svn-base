set -e

# IO
DIR_BASE=$(dirname $0) # Dir script resides in
DIR_BASE=$(readlink -f $DIR_BASE) # Obtain the full path, otherwise nodes go crazy
DIR_CUR=${PWD} # Dir script is called from

function loadConfig ()
{
	source $DIR_BASE/etc/default.sh # Always load default variables
	LOG_VER=$(svn info $DIR_BASE | grep URL | rev | cut -d '/' -f 1 | rev)
	LOG_REV=$(svn info $DIR_BASE | grep 'Last Changed Rev:' | awk '{ print $4 }')
	echo '# RoDa '$LOG_VER' r'$LOG_REV >> $LOG_CONF
	echo '# Run date: '$(date +"%d/%m/%y") >> $LOG_CONF
	echo '' >> $LOG_CONF
	
	echo '# Default config:' >> $LOG_CONF
	cat $DIR_BASE/etc/default.sh >> $LOG_CONF

	if [ -f $1 ]; then
		FILE_CONFIG=$1
	elif [ -f  $DIR_BASE/etc/$1.sh ]; then
		FILE_CONFIG=$DIR_BASE/etc/$1.sh
	else
		echo Config file not found: $1 nor $DIR_BASE/etc/$1.sh
		exit 1;
	fi
	echo Loading $FILE_CONFIG
	
	source $FILE_CONFIG

	echo '' >> $LOG_CONF
	echo '' >> $LOG_CONF
	echo '# Config used:' >> $LOG_CONF
	cat $FILE_CONFIG >> $LOG_CONF
}

export FILE_INPUT=$(readlink -f $1)
export DIR_OUTPUT=$(readlink -f $2)

STAMP=`date +%s` # Unique identifier based on the second the script is called, may need better ideas

# Name fetching, bit pointless (and crap) perhaps, but fun!
NAME_SAMPLE=$(basename "$FILE_INPUT")
EXT_SAMPLE=${NAME_SAMPLE#*.}
NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
export NAME_SAMPLE=$(echo $NAME_SAMPLE | cut -d _ -f 1)
DIR_SAMPLE=$(dirname "$FILE_INPUT")
export FILE_OUTPUT=$DIR_OUTPUT/$NAME_SAMPLE

# Debugging info
echo 'Running Sample:'
echo -e '\tInput file:\t' $FILE_INPUT
echo -e '\tFile basename:\t' $NAME_SAMPLE
echo -e '\tIn directory:\t' $DIR_SAMPLE
echo -e '\tOutput dir:\t' $DIR_OUTPUT
# No shit sherlock...

# Create required dirs
DIR_SGELOG=$DIR_OUTPUT/sge_log
for DIR in $DIR_OUTPUT $DIR_SGELOG
do
	if [ ! -d $DIR ]
	then
		echo -e '\t  Creating:\t' $DIR
		mkdir $DIR
	fi
done

# Catch errors, remove all jobs if an error occurred
trap onexit 1 2 3 15 ERR
function onexit() {
	local exit_status=${1:-$?}
	echo Exiting $0 with $exit_status
	echo 'Removing all added jobs from queue'
	qdel ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_* # Arrr...may remove other jobs if added at the very same second using this ID
	exit $exit_status
}

LOG_CONF=$FILE_OUTPUT.conf$STAMP
loadConfig $3

# Simple pipeline function: Add list of steps to run and it will fix dependencies by itself
function runpipe()
{
	PIPELINE=$1
	for STEP in ${PIPELINE[@]}
	do
		qsub \
			-hold_jid ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${PREVSTEP} \
			-V \
			-N ${RODA_ID}_${NAME_SAMPLE}_${STAMP}_${STEP} \
			-e $DIR_SGELOG/${STEP}.e${STAMP} \
			-o $DIR_SGELOG/${STEP}.o${STAMP} \
			$DIR_BASE/bin/$STEP.sh
		PREVSTEP=$STEP
	done
}

source $DIR_BASE/pip/${RODA_PIPELINE}.sh
