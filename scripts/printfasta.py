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

def handle_infile(infile, identifier, start=0, end=None, width=50):
      
    filepath = os.path.abspath(infile)    
    filename = os.path.basename(filepath)
    filedir = os.path.dirname(filepath)
    (base, ext) = os.path.splitext(filename)
    ext = ext[1:]
    logging.debug(f' filepath={filepath} filedir={filedir} base={base} ext={ext} filename={filename}')
    
    records = list( SeqIO.parse(infile, 'fasta'))
    logging.debug(f'got {len(records)} sequences')
    
    for srec in records:
        if identifier in srec.id:
            if end is None:
                end = len(srec.seq)
            print( str( srec[start:end].seq ))
             
    
               

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

    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='fasta input file')

    parser.add_argument('identifier', 
                        metavar='identifier', 
                        type=str, 
                        help='sequence id substring')

    parser.add_argument('start', 
                        metavar='start', 
                        type=int, 
                        default=0,
                        help='start')    

    parser.add_argument('end', 
                        metavar='end', 
                        type=int, 
                        default=None,
                        help='start')   
    
    args= parser.parse_args()
   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    handle_infile( args.infile, args.identifier, args.start, args.end)