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

echo bedtools bamtofastq -i $2 -fq $3 -fq2 $4 
time bedtools bamtofastq -i $2 -fq $3 -fq2 $4
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"

echo "*********DONE*************************"
date
echo "*********END***************************"

exit $RET
