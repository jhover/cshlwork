#!/bin/bash 
#   Usage: rungofocus.sh  <setup> <basefile> <outdir> 
#
#   args =  $(setup) 
#           $(basefile)
#           $(outdir)
# job name
#$ -N gofocus
#
# job indexes for array all jobs, but only run 10 at a time for disk quota 
#$ -t 1-39
#
# processes per job
#$ -pe threads 4
#
#$ -wd /grid/gillis/home/hover/play/jones
#
# Per-thread memory request. 
#$ -l m_mem_free=3G
# 
#$ -l gpu=1 
#
#

COMMON=~/git/cshl-work/scripts/jones/common.sh 
CMD=~/git/gofocus/gofocus/pytorch_goterm_pred.py 

echo "*********START*************************"
date
. $COMMON

echo "*********JOB*************************"
echo "Args are $@ Numargs is $#"
if [ $# -ne 3 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: rungofocus.sh  <setup> <basefile> <outdir>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

basefile=$2
# $SGE_TASK_ID
TASKID=gettaskid

filebase=`head -$SGE_TASK_ID $2 | tail -1 `
echo "Filebase is $filebase"

outdir=$3
echo "Outdir is $outdir"
infile="$outdir/fasta/${filebase}_hiprio.tfa "
outfile="$outdir/${filebase}_hiprio.predout"

echo "infile is $infile"
echo "outfile is $outfile"

#
#    time ~/git/gofocus/gofocus/pytorch_goterm_pred.py -d <species<_hiprio.tfa <species>_hiprio.predout
echo $CMD -v $infile $outfile
time $CMD -v $infile $outfile

RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"
date
echo "*********END***************************"
exit $RET

