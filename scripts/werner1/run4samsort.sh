#!/bin/bash 
#   Usage: runsamtools <setup> <infile> <outfile>
#
#   args = $(setup) 
#          $(outdir)/$(filebase).Aligned.out.bam  
#          $(outdir)/$(filebase).Aligned.sortedByCoord.out.bam
#   request_cpus = 10
#   request_memory = 20480
#
#
#
echo "*********START*************************"
date

echo "*********NODE*************************"
hostname -f
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
if [ $# -ne 3 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: runsamtools <setup> <infile> <outfile>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

INTMP=`mktemp -p ./`

echo "Staging in..."
echo cp -v $2 $INTMP
time cp -v $2 $INTMP

echo "Running job..."
TMPFILE=`mktemp -p ./`
echo "TMPFILE is $TMPFILE"
echo samtools sort -m 2G -o $TMPFILE -O bam  -@ 20 $INTMP
time samtools sort -m 2G -o $TMPFILE -O bam  -@ 20 $INTMP
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo samtools index $TMPFILE 
time samtools index $TMPFILE 
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Job command Return code was $RET"

echo "********STAGEOUT**********************"
echo mv -v $TMPFILE $3
time mv -v $TMPFILE $3

echo mv -v $TMPFILE.bai $3.bai
time mv -v $TMPFILE.bai $3.bai 

echo "*********DONE*************************"
date
echo "*********END***************************"
exit $RET

