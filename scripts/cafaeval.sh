#!/bin/bash
#
#
SPECLIST="HUMAN MOUSE RAT BOVIN PIG DANRE DROME CAEEL ARATH MAIZE YEAST SCHPO ECOLI PSEAI  DICDI MYCGE BACCR SALTY"
PROG=~/git/cshl-work/fastcafa/fastcafa.py
CONF=~/git/cshl-work/etc/fastcafa.conf
NUMSEQ=10

for SPECIES in $SPECLIST; do 
	TFAOUT=~/play/cafa4/sp_species.$SPECIES.test.$NUMSEQ.tfa
	PREDOUT=~/play/cafa4/sp_species.$SPECIES.phmmer.test.$NUMSEQ.prediction.csv
	EVALOUT=~/play/cafa4/sp_species.$SPECIES.test.$NUMSEQ.evaluate.csv

	echo "time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT"  
	time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT

	echo "time $PROG -C -c $CONF phmmer -i $TFAOUT -o $PREDOUT"
	time $PROG -C -c $CONF expression -i $TFAOUT -o $PREDOUT

	echo "time $PROG -C -c $CONF evaluate -p $PREDOUT -o $EVALOUT"
	time $PROG -C -c $CONF evaluate -p $PREDOUT -o $EVALOUT
done


