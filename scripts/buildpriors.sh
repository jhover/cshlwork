#!/bin/bash
#
# The CAFA target species....
SPECLIST="ARATH BACCR BOVIN CAEEL DANRE DICDI DROME ECOLI HUMAN MAIZE MOUSE MYCGE PIG PSEAI SCHPO RAT SALTY YEAST "
PROG=~/git/cshl-work/fastcafa/fastcafa.py
CONF=~/git/cshl-work/etc/fastcafa.conf
OUTDIR=~/play/cafa4

for SPECIES in $SPECLIST; do 
	echo "Building prior for $SPECIES..." 
	echo "time $PROG -c $CONF build_prior -s $SPECIES -o $OUTDIR/uniprot.prior.$SPECIES.csv" 
	time $PROG -c $CONF build_prior -s $SPECIES -o $OUTDIR/uniprot.prior.$SPECIES.csv    
done

echo "Building global prior..."
echo "time $PROG -c $CONF build_prior -o $OUTDIR/uniprot.prior.GLOBAL.csv"
time $PROG -c $CONF build_prior -o $OUTDIR/uniprot.prior.GLOBAL.csv

