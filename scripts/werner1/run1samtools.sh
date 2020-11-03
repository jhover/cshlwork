#!/bin/bash 
#   Usage: runsamtools <setup> <infile> <outfile>
#
# args = $(setup) 
#        $(indir)/$(filebase).Aligned.sortedByCoord.out.patched.md.bam 
#        $(outdir)/$(filebase).readSorted.bam 
# request_cpus = 12
# request_memory = 16348
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

echo "Running job..."
TMP1=`mktemp -p ./`
INTMP=`mktemp -p ./`


echo samtools sort -m 2G -o $3 -O bam -n  -@ 25 $2
time samtools sort -m 2G -o $3 -O bam -n -@ 25 $2
RET=$?
echo "Job command Return code was $RET"

echo "*********DONE*************************"
date
echo "*********END***************************"
exit $RET

