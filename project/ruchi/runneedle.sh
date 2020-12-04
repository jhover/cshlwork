#!/bin/sh
# This is a simple example of a SGE batch script

# request Bourne shell as shell for job
#$ -S "/bin/bash -l"

#  Working dir. Also where error and out are put. 
#$ -wd /grid/gillis/home/hover/play/ruchi/pairwise1

#   job name
#$ -N needle

#$ -t 6-201

# Per-processor memory request. 
#$ -l h_vmem=2G

date
whoami
hostname -f
pwd
. ~/.bash_profile
conda activate emboss
~/git/cshl-work/scripts/pair_needle.py protein_uniprot_pairwise.tsv.$SGE_TASK_ID protein_uniprot_pairwise.tsv.$SGE_TASK_ID.out
# print date and time again
date
