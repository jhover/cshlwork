#!/usr/bin/env python
#
#  Split file by lines into N pieces.
#  splitfile.py  <file> <N>

import argparse
import os
import sys
import logging

from itertools import combinations

gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

def do_comb(infile, outfile=None):
    plist = []
    logging.debug(f"Opening infile {infile}")
    i = 0
    with open(infile) as f:
        for i, l in enumerate(f):
            plist.append(l.strip())
    f.close()
    logging.debug(f"Handled {i} items...")
  
    comb = combinations(plist, 2)
    comblist = list(comb)
    logging.debug(f"Generated {len(comblist)} combinations...")
    
    if outfile is None:
        outfile = f'{infile}.pairwise.tsv'
    o = open(outfile,'w')
    logging.debug(f"Opened outfile {outfile}")
    for (p1, p2 ) in comblist:
        o.write(f"{p1}\t{p2}\n")
    o.close()
  
  
  
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
                        help='input file')

    #parser.add_argument('outfile', 
    #                    metavar='outfile', 
    #                    type=str,
    #                    default=None,
    #                    required=False, 
    #                    help='pairwise info. ')
    
    args= parser.parse_args()

   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    do_comb(args.infile, None)