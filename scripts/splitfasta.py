#!/usr/bin/env python
#   
#  Takes .fasta files and splits them, one entry per file, 
#  named after ID (first string after ">"
#

import argparse
import json
import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

#import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from cshlwork.utils import *

def handle_infile(infile, maxseq=1, outdir = None):
      
    filepath = os.path.abspath(infile)    
    filename = os.path.basename(filepath)
    filedir = os.path.dirname(filepath)
    (base, ext) = os.path.splitext(filename)
    ext = ext[1:]
    logging.debug(f' filepath={filepath} filedir={filedir} base={base} ext={ext} filename={filename}')
    
    if outdir is None:
        outdir = os.getcwd()
    else:
        outdir = os.path.abspath(outdir)    
    
    logging.debug(f'outdir = {outdir}')  
    records = list( SeqIO.parse(infile, 'fasta'))
    logging.debug(f'got {len(records)} sequences')
    
    splitnum = 0
    seqnum = 0
    handled = 0
    currentlist = []
    
    splitsets = {}
    for srec in records:
        currentlist.append(srec)
        if len(currentlist) >= maxseq:
            splitsets[splitnum] = currentlist
            splitnum += 1
            currentlist = []
    # get last items. 
    splitsets[splitnum] = currentlist
    
    lenlist = []
    for k,v in splitsets.items():
        lenlist.append(len(v))        
    logging.debug(f'seqs split into {len(splitsets)} sets with lengths={lenlist}')

    for k,v in splitsets.items():
        outfile = f'{filedir}/{base}.{k}.{ext}'
        logging.info(f'writing {len(v)} seqs to {outfile}')
        with open(outfile, 'w') as of:
            SeqIO.write(v, of, 'fasta')
        
    logging.debug(f'done writing fasta files. ')
               

if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('maxseq', 
                        metavar='maxseq', 
                        type=int, 
                        default=1,
                        help='max number of sequences per file. [one seq per output file]')    
    
    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='fasta input file')

    parser.add_argument('-o','--outdir', 
                        metavar='outdir',
                        required=False,  
                        type=str,
                        default=None,
                        help='output dir [<cwd>] ')
    
    args= parser.parse_args()
   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    handle_infile(args.infile, maxseq = args.maxseq, outdir=args.outdir)