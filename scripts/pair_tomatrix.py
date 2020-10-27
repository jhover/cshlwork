#!/usr/bin/env python
#
#  comsume tsv table, 
#  create matrix from 2 columns. 
#

import argparse
import os
import sys
import logging
import traceback
import subprocess

import pandas as pd

gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

gitpath=os.path.expanduser("~/git/pyEGAD")
sys.path.append(gitpath)

from egad.egad import *


def do_matrix(infile, outfile):
    columns = ['p1','p2','len','ident','simil','gaps','score','pident','psimil']
    df = pd.read_csv(infile,sep='\t')
    df.columns=columns
    logging.debug(df)
    #genelist = make_genelist(df.values)
    #logging.debug(genelist)
    #logging.debug(f"genelist length= {len(genelist)}")

    matrix = df.pivot(index='p1', columns='p2', values='psimil')
    logging.debug(f"matrix=\n{matrix}")
    logging.debug(f"matrix shape= {matrix.shape}")
    matrix.to_csv(outfile, sep='\t', float_format='%.3f' )
    logging.debug(f"Wrote {outfile}")

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
                        help='tsv table.')

    parser.add_argument('outfile', 
                        metavar='outfile', 
                        type=str, 
                        help='matrix file')
    
    args= parser.parse_args()

   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    do_matrix(args.infile, args.outfile)