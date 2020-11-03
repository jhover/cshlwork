#!/bin/bash 
#   Usage: runstar <setup> <genomedir> <infile1> <infile2> <outdir> <filebase>
#
#   args =  $(setup) 
#           $(genomedir) 
#           $(outdir)/$(filebase).end1.fq  
#           $(outdir)/$(filebase).end2.fq  
#           $(outdir) 
#           $(filebase)
#   request_cpus = 16
#   request_memory = 51200
#

echo "*********START*************************"
date
echo "*********NODE*************************"
hostname -f
WD=`pwd`
echo "Working Directory: $WD"
cat /etc/redhat-release
CV=`condor_version | tr -d "\n"`
echo "Condor Version: $CV"
NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
echo "Processors: $NPROC "
KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
MBMEM=`expr $KMEM / 1000`
echo "Memory MB: $MBMEM"



echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 6 ]; then
    echo "Incorrect number of arguments."
    "Usage: runstar <setup> <genomedir> <infile1> <infile2> <outdir> <filebase>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

echo "Staging in..."
INTMP1=`mktemp -p ./`
INTMP2=`mktemp -p ./`
echo cp -v $3 $INTMP1
time cp -v $3 $INTMP1
echo cp -v $4 $INTMP2
time cp -v $4 $INTMP2

echo "Running job..."
echo STAR  --genomeDir $2 --readFilesIn $INTMP1 $INTMP2 --runThreadN 20 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM Unsorted --quantMode GeneCounts  
time STAR  --genomeDir $2 --readFilesIn $INTMP1 $INTMP2 --runThreadN 20 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM Unsorted  --quantMode GeneCounts
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"

echo "********STAGEOUT**********************"
echo "Contents of working dir:"
ls -alh $WD
KEEPFILES="Aligned.out.bam Log.out Log.final.out Log.progress.out ReadsPerGene.out.tab SJ.out.tab"
for f in $KEEPFILES; do
    echo mv -v $WD/$f $5/$6.$f
    mv -v $WD/$f $5/$6.$f
done

echo "*********END***************************"
exit $RET

