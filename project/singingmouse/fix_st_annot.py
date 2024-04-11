#!/usr/bin/env python
#
# Fix annotation nesting in Singing Mouse genome files. 
# Restore hierarchy, include qualifiers. 
#

import argparse
import logging
import os
import sys
import traceback


gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *
from genome.genome import *


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
                        type=str,
                        help='GTF annotation file')
    
    parser.add_argument('outfile',
                        type=str,
                        help='GTF annotation file')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    annot_list = load_annotrecords(args.infile)
    fix_annot(annot_list)
