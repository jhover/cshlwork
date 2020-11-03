#!/bin/bash 
#   Usage: runsamgatvigv <setup> <infile> <outfile> <outvcf> <outwig> <outcut>
# 
#   args = $(setup) 
#          $(outdir)/$(filebase).Aligned.sortedByCoord.out.bam  
#          $(outdir)/$(filebase).chrX.split.filtered.wig $(outdir)/$(filebase).chrX.filtered.vcf 
#          $(outdir)/$(filebase).chrX.intersect.wig $(outdir)/$(filebase).chrX.cut
#   request_cpus = 1
#   request_memory = 10240
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
if [ $# -ne 6 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: runsamgatvigv <setup> <infile> <outfile> <outvcf> <outwig> <outcut>"
    exit 1
fi

echo "Running job..."
echo "Running setup from $1"
. $1
IGV=`igvtools version`
echo "IGV version: $IGV"
echo "PATH=$PATH"

chr="chrX"
chr_fa="/data/jwerner/data/annot_genomes/GRCh38_Gencode25/$chr.fa"
genome_fa="/data/jwerner/data/annot_genomes/GRCh38_Gencode25/GRCh38.p7.genome.fa"
snpdb="/data/jwerner/data/snpdb/chrX.snpdb.txt"


echo "Staging in..."
INTMP=`mktemp -p ./`
echo cp -v $2 $INTMP
time cp -v $2 $INTMP


echo "samtools view -b $INTMP  chrX > $chr.bam"
samtools view -b $INTMP chrX > $chr.bam
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "samtools view -b -q 10 $chr.bam > $chr.filt.bam"
samtools view -b -q 10 $chr.bam > $chr.filt.bam
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo gatk AddOrReplaceReadGroups -I=$chr.filt.bam -O=$chr.rg.bam -SO=coordinate -RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample
gatk AddOrReplaceReadGroups -I=$chr.filt.bam -O=$chr.rg.bam -SO=coordinate -RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Marking duplicates..."
echo gatk MarkDuplicates -I=$chr.rg.bam -O=$chr.dedupped.bam -CREATE_INDEX=true -VALIDATION_STRINGENCY=SILENT -M=output.metrics
gatk MarkDuplicates -I=$chr.rg.bam -O=$chr.dedupped.bam -CREATE_INDEX=true -VALIDATION_STRINGENCY=SILENT -M=output.metrics
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Splitting and trimming..."
echo gatk SplitNCigarReads -R $genome_fa -I $chr.dedupped.bam -O $chr.split.filtered.bam
gatk SplitNCigarReads -R $genome_fa -I $chr.dedupped.bam -O $chr.split.filtered.bam
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Haplotype calling (?)..."
echo gatk HaplotypeCaller -R $genome_fa -L $chr -I $chr.split.filtered.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0  -O $chr.filtered.vcf
gatk HaplotypeCaller -R $genome_fa -L $chr -I $chr.split.filtered.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0  -O $chr.filtered.vcf
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

TMPFILE=`mktemp -p ./  tmp.XXXXXX.wig`
echo "TMPFILE is $TMPFILE"
echo "Counting SNPs..."
echo igvtools count -z 0 -w 1 --bases --strands read $chr.split.filtered.bam $TMPFILE $chr_fa
igvtools count -z 0 -w 1 --bases --strands read $chr.split.filtered.bam $TMPFILE $chr_fa
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Additional filtering/cutting..."
grep chrX chrX.filtered.vcf | cut -f2-5 > A
cut -f1 A > B
grep -Fw -f B $snpdb > D
grep -Fw -f D $TMPFILE > chrX.intersect.wig
grep -Fw -f D A > chrX.cut 


echo "**************FINALIZE************"
ls -alh

# handle TMPFILE (which is chrX.split.filtered.wig) 
echo mv -v $TMPFILE $3
time mv -v $TMPFILE $3

# handle 
echo mv -v chrX.filtered.vcf $4
mv -v chrX.filtered.vcf $4

echo mv -v chrX.intersect.wig $5
mv -v chrX.intersect.wig $5

echo mv -v chrX.cut $6
mv -v chrX.cut $6


echo "*********DONE*************************"
date
echo "*********END***************************"
sleep 30
exit $RET

