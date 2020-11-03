#!/bin/bash 
#   Usage: runbedtools <setup> <infile> <outfile1> <outfile2>
#
#   args =  $(setup)  $(outdir)/$(filebase).readSorted.bam $(outdir)/$(filebase).end1.fq  $(outdir)/$(filebase).end2.fq
#   request_cpus = 1
#   request_memory = 2048

#
#
#


echo "*********START*************************"
date

echo "*********NODE*************************"
hostname -f
cat /etc/redhat-release
NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
echo "Processors: $NPROC "
KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
MBMEM=`expr $KMEM / 1000`
echo "Memory MB: $MBMEM"
CV=`condor_version | tr -d "\n"`
echo "Condor Version: $CV"
NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
echo "Processors: $NPROC "
KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
MBMEM=`expr $KMEM / 1000`
echo "Memory MB: $MBMEM"


echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 4 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: runbedtools <setup> <infile> <outfile1> <outfile2>"
     
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

echo "Running job..."
TMP1=`mktemp -p ./`
TMP2=`mktemp -p ./`

INTMP=`mktemp -p ./`

echo "Staging in..."
echo cp -v $2 $INTMP
time cp -v $2 $INTMP


echo bedtools bamtofastq -i $INTMP -fq $TMP1 -fq2 $TMP2 
time bedtools bamtofastq -i $INTMP -fq $TMP1 -fq2 $TMP2
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"

echo "********STAGEOUT**********************"
echo mv -v $TMP1 $3
mv -v $TMP1 $3
echo mv -v $TMP2 $4
mv -v $TMP2 $4

echo "*********DONE*************************"
date
echo "*********END***************************"

exit $RET
