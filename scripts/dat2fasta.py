#!/usr/bin/env python
# takes in uniprot .dat 
# emits fasta
# 

import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from protlib.uniprot import *


def dat_to_fasta( infile, outfile):
    logging.debug(f"infile={infile}  outfile={outfile}")
    uplist = parse_uniprot_dat(infile)
    logging.debug(f"writing fasta to {outfile}")
    write_tfa_file(uplist, outfile)
    logging.debug("done.")
    



if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
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
                        help='Uniprot .dat file') 


    parser.add_argument('outfile', 
                        metavar='outfile', 
                        type=str, 
                        help='Fasta .tfa file.')    
    
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
    
    dat_to_fasta(os.path.expanduser(args.infile),
              os.path.expanduser(args.outfile))
