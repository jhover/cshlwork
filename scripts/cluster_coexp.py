#!/usr/bin/env python
import logging
import os
import sys

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import pandas as pd
import numpy as np

from cshlwork.exputils import *
from cshlwork.utils import *

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

    parser.add_argument('-T', '--test', 
                        action="store_true", 
                        dest='test', 
                        help='make small subset of input for quick test.')

    parser.add_argument('-i','--infile', 
                        metavar='infile',
                        required=False,  
                        type=str,
                        default=os.path.expanduser('~/data/cococonet/yeast_AggNet.hdf5'),
                        help='hd5 input file ')

    parser.add_argument('-t','--threshold', 
                        metavar='threshold',
                        required=False,  
                        type=float,
                        default=0.95,
                        help='only consider coefficients above threshold')

    parser.add_argument('-s','--sense', 
                        metavar='sense',
                        required=False,  
                        type=str,
                        default='above',
                        help='keep above or below threshold')

    
    parser.add_argument('-o','--outfile', 
                        metavar='outfile',
                        required=False,  
                        type=str,
                        default=None,
                        help='summarized .tsv output ')

       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
        
    outdf = cluster_coexp(exphd5=args.infile, threshold=args.threshold, test=args.test, sense=args.sense)
    logging.debug(f'outdf = \n{outdf}')
    if args.outfile is not None:
        write_df(outdf, args.outfile)
    else:
        print(outdf)
