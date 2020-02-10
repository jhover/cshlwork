#!/bin/bash
#
#  Run eval (phmmer and expression) on input files. 
#  Output .csv predictions. 
#
DATE=`date "+%Y%m%d"`
TESTDIR=~/data/cafa4/TargetFiles
OUTDIR=~/play/cafa4/$DATE
PROG=~/data/git/cshl-work/fastcafa/fastcafa.py
CONF=~/data/git/cshl-work/etc/fastcafa.conf
METHODS="prior"
#DEBUG="-d "
DEBUG=""

echo "Making $OUTDIR ..."
mkdir -p $OUTDIR

echo "Running on all targets..."
for TFA in `ls $TESTDIR/*.tfa`; do
	echo "Handling $TFA..."
	echo "###############################################"
	FILENAME=`basename $TFA`
	EXTENSION="${FILENAME##*.}"
	FILEBASE="${FILENAME%.*}"
	#echo "$FILENAME $FILEBASE $EXTENSION"
	for METHOD in $METHODS; do
		# echo "method is $METHOD"
		PREDOUT=$OUTDIR/$FILEBASE.$METHOD.csv
		
		echo "running $METHOD ..."
		echo "time $PROG -C $DEBUG -c $CONF $METHOD -i $TFA -o $PREDOUT "
		echo ""
		time $PROG -C $DEBUG -c $CONF $METHOD -i $TFA -o $PREDOUT 

	done
	echo  "###############################################"
done
