#!/usr/bin/env python
#
# read genome fasta(s) 
# read annotation files
# write out spliced mrnas to fasta
# and/or to tsv. 
#
#
#
import argparse
import os
import sys
import logging
import traceback

import pandas as pd

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from genome.genome import *
from cshlwork.utils import *


def handle_genome(genome_file, annot_file, outfile=None):

    logging.info(f'building transcriptome  ')
    seq_dict = load_seqrecords(genome_file)
    logging.info(f'got {len(seq_dict)} records from {genome_file}') 
    
    annot_list = load_annotrecords(annot_file)
    logging.info(f'got {len(annot_list)} annotations from {annot_file} ')
    
    transcript_list = get_genes(annot_list, srecord_dict=None)
    logging.info(f'got {len(transcript_list)} merged sequences.')
    
    if outfile is not None:
        logging.info(f'writing {len(transcript_list)} items to {outfile}')
        write_seqfeature_fasta(transcript_list, seq_dict, outfile)
        logging.info('done')
    return transcript_list
    
  
  


if __name__=='__main__':

    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARNING)

    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('-g','--genome', 
                        metavar='genome',
                        required=False,  
                        type=str,
                        default=os.path.expanduser('~/data/genomes/Mus_musculus/genome.fa'), 
                        help='Genome .fasta file') 

    parser.add_argument('-a','--annot', 
                        metavar='annot',
                        required=False,  
                        type=str,
                        default=os.path.expanduser('~/data/genomes/Mus_musculus/annotation.gtf'), 
                        help='Annotation .gtf or .gff file') 

    parser.add_argument('-o','--outfile', 
                        metavar='outfile',
                        required=False,
                        default='transcriptome.fasta', 
                        type=str, 
                        help='out file for all transcripts. ')    

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    tlist = handle_genome(args.genome, args.annot, outfile=args.outfile)



