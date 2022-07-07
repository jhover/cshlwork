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
    current = None
    processed = 0
    
    filepath = os.path.abspath(infile)    
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)
    logging.debug(f'base={base} ext={ext} filename={filename} filepath={filepath}')
    
    if outdir is None:
        outdir = os.getcwd()
    else:
        outdir = os.path.abspath(outdir)    
    
    logging.debug(f'outdir = {outdir}')  
    records = list( SeqIO.parse(infile, 'fasta'))
    
    logging.debug(f'got {len(records)} sequences')
    
    split=0
    seqnum = 0
    
    while seqnum < len(records):
        
        
        
        
        
        
        
        seqnum += 1
    
    #    print(record.description)
    
    
    with open(infilepath, 'r') as f:
        try:
            while True:
                line = f.readline().strip()
                if line == '':
                        break
                
                if line.startswith(">"):
                    if current is not None:
                        write_current(current)
                        processed += 1
                    current = []
                    current.append(line)

                else:
                    current.append(line)
        except Exception as e:
            traceback.print_exc(file=sys.stdout)

    logging.debug(f'done writing {processed} fasta files. ')
               

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