#!/usr/bin/env python
#
#  Split file by lines into N pieces.
#  splitfile.py  <file> <N>

import argparse
import os
import sys
import logging
import traceback
import subprocess


gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

def do_water(infile, outfile):
    logging.debug(f"processing {infile} to {outfile}...")
    with open(infile) as f:
        for i, l in enumerate(f):
            p1, p2 = l.split('\t')
            logging.debug(f"p1={p1} p2={p2}")
    infile.close()
    
def run_water(p1, p2):
    cmd = f'water -brief -gapopen 10.0 -gapextend 0.5 -stdout -auto -aformat score up:{p1} up:{p2}'
    logging.debug(f"command is {cmd}")
    





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
                        help='a .fasta sequence file')

    parser.add_argument('outfile', 
                        metavar='outile', 
                        type=str, 
                        help='pairwise info. <p1> <p2> <len> <dist>')
    
    args= parser.parse_args()

   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    do_water(args.infile, args.outfile)
    
    