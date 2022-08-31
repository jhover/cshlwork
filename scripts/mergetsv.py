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
    outdf = None 
    for infile in infiles:
        basename = os.path.basename(infile)
        # print(f"{infile}")
        try:
            df = pd.read_csv(infile, index_col=[0], sep='\t', comment="#")
            if outdf is None:
                cols = list(df.columns)
                logging.debug(f'got columns: {cols}')
                outdf = pd.DataFrame(columns=cols)
            outdf = pd.concat([outdf, df], ignore_index=True )
        
        except:
            logging.warn(f'something went wrong with {infile}')
    
    outdf.drop_duplicates(inplace=True)
    outdf = outdf.reset_index(drop=True)
    logging.debug(f'final columns: {outdf.columns}')
    return outdf


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
        


