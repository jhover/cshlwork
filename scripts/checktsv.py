#!/usr/bin/env python
#
#  checks TSVs for sanity. column number, and value number for each row. 
#
import argparse
import logging
import os
import sys
#gitpath=os.path.expanduser("~/git/cshlwork")
#sys.path.append(gitpath)

from collections import defaultdict

import pandas as pd
import numpy as np
import traceback


def check_all(infiles):
    outdf = None 
    for infile in infiles:
        basename = os.path.basename(infile)
        logging.debug(f"checking '{infile}' ...")
        try:
            df = pd.read_csv(infile, index_col=[0], header=0, sep='\t', comment="#",  on_bad_lines='error')
            cols = list(df.columns)
            logging.debug(f'got columns: {cols}')
        except:
            logging.warning(f'something went wrong with {infile}')
            logging.error(traceback.format_exc(None))
        logging.info(f'tsv reads in as DF...')
        
        lendict = defaultdict(int)
        
        with open(infile) as f:
            nlines = 0
            for line in f:
                fields = line.split('\t')
                lendict[len(fields)] += 1
                nlines += 1
        logging.info(f'handled {nlines} lines.')
        logging.info(f'column counts: {dict(lendict)}')
        


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
    
    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='TSV [ TSV TSV ...] ')

    args= parser.parse_args()
    # This is diagnostic. Always run in debug. 
    logging.getLogger().setLevel(logging.DEBUG)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    
    check_all(args.infiles)

        


