DATE=$(date +"%y%m%d")

PATH_DIAG=/illumina/diagnostics
PATH_DATA=$PATH_DIAG/validation/data
PATH_OUTP=$PATH_DIAG/validation/output

SCRIPT_RODA=$PATH_DIAG/roda/roda_mult.sh

for SET in $PATH_DATA/*
do
	NAME_FOLDER=$(echo $SET | rev | cut -d '/' -f 1 | rev)
	NAME_TYPE=$(echo $NAME_FOLDER | cut -d _ -f 1)

	$SCRIPT_RODA \
		$SET \
		$PATH_OUTP/${NAME_FOLDER}_$DATE \
		$NAME_TYPE
done
