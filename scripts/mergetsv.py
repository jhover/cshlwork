#!/usr/bin/env python
#
#  Takes N .tsv file and outputs merged tsv. 
#
import argparse
import logging
import os
import sys
#gitpath=os.path.expanduser("~/git/cshlwork")
#sys.path.append(gitpath)

import pandas as pd
import numpy as np
import traceback


def merge_all(infiles):
    topdf = None 
    for infile in infiles:
        basename = os.path.basename(infile)
        # print(f"{infile}")
        try:
            df = pd.read_csv(infile, index_col=[0], sep='\t', comment="#")
            if topdf is None:
                topdf = pd.DataFrame(columns=list(df.columns))   
            topdf = topdf.append(df, ignore_index=True, sort=True)
        
        except:
            logging.warn(f'something went wrong with {infile}')
    
    topdf.drop_duplicates(inplace=True)
    topdf = topdf.reset_index(drop=True)
    return topdf


def write_tsv(df, outfile=None):
    if outfile is None:       
        outfile = sys.stdout
    logging.debug(f'writing {len(df)} lines output to {outfile}')      
    df.to_csv(outfile, sep='\t', index=False)
    


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

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='out file. stdout if not given.')  
    
    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='TSV [ TSV TSV ...] ')

    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    
    df = merge_all(args.infiles)
    logging.info(f'got df with {len(df)} rows.')
    write_tsv(df, args.outfile)
        


