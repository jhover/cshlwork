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

echo "Staging in..."
echo cp -v $2 $INTMP
time cp -v $2 $INTMP

echo samtools sort -m 2G -o $TMP1 -O bam -n  -@ 25 $INTMP
time samtools sort -m 2G -o $TMP1 -O bam -n -@ 25 $INTMP
RET=$?
echo "Job command Return code was $RET"

echo "********STAGEOUT**********************"
echo mv -v $TMP1 $3
mv -v $TMP1 $3

echo "*********DONE*************************"
date
echo "*********END***************************"
exit $RET

