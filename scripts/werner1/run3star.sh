#!/bin/bash 
#   Usage: runstar <setup> <genomedir> <infile1> <infile2> <outdir> <filebase>
#
#   args =  $(setup) 
#           $(basefile)
#           $(genomedir) 
#           $(outdir)
#						$(outdir)/$(filebase).end1.fq  
#           			$(outdir)/$(filebase).end2.fq  
#           			$(outdir) 
#           			X $(filebase)
#   request_cpus = 16
#   request_memory = 51200
#
# job name
#$ -N run3star
#
# job indexes for array all jobs, but only run 10 at a time for disk quota 
#$ -t 1-4 -tc 2
#
# processes per job
#$ -pe threads 16 
#
#$ -wd /grid/gillis/data/hover/work/werner1/
#
# Per-thread memory request. 
#$ -l m_mem_free=3G
# 
COMMON=~/git/cshl-work/scripts/werner1/common.sh 

NUMTHR=16


echo "*********START*************************"
date
. $COMMON

echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 4 ]; then
    echo "Incorrect number of arguments."
    "Usage: runstar <setup> <genomedir> <infile1> <infile2> <outdir> <filebase>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

basefile=$2
# $SGE_TASK_ID
filebase=`head -$SGE_TASK_ID $2 | tail -1 `
echo "Filebase is $filebase"

infile1="$4/$filebase.end1.fq "
infile2="$4/$filebase.end2.fq"

genomdir=$3

echo "Running job..."

mkdir $filebase
cd $filebase

echo STAR  --genomeDir $genomedir --readFilesIn $infile1 $infile2 --runThreadN $NUMTHR --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM Unsorted --quantMode GeneCounts  
time STAR  --genomeDir $genomedir --readFilesIn $infile1 $infile2 --runThreadN $NUMTHR --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM Unsorted  --quantMode GeneCounts
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"

echo "********STAGEOUT**********************"
echo "Contents of working dir:"
ls -alh
KEEPFILES="Aligned.out.bam Log.out Log.final.out Log.progress.out ReadsPerGene.out.tab SJ.out.tab"
for f in $KEEPFILES; do
    echo mv -v $f ../$basefile.$f
    mv -v $f ../$basefile.$f
done

cd ..

echo "*********END***************************"
exit $RET
