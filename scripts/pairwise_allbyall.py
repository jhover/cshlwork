#!/usr/bin/env python
#
# run phmmer against fasta file all by all.  
# produce csv of pairwise match alignment. 
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
from protlib import phmmer
from cshlwork import utils

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
    
    parser.add_argument('-o','--outfile', 
                        metavar='outfile',
                        required=False,  
                        type=str,
                        default=None,
                        help='summarized data output to TSV file. ')

    parser.add_argument('fastafile', 
                        metavar='fastafile', 
                        type=str, 
                        help='')


    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    pc = phmmer.get_default_config()
    logging.info(f"Running phmmer query={args.fastafile} db={args.fastafile}")
    pdf = phmmer.get_phmmer_df(pc, args.fastafile, args.fastafile)     
    if args.outfile is not None:
        merge_write_df(pdf, args.outfile)
    else:
        print(pdf)    
