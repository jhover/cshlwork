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

    parser.add_argument('-i','--infile', 
                        metavar='infile',
                        required=False,  
                        type=str,
                        default=None,
                        help='hd5 input file ')
    
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
        
    outdf = cluster_coexp(exphd5='~/data/cococonet/yeast_AggNet.hdf5', threshold=0.95, test=True)
    logging.debug(f'outdf = \n{outdf}')
    if args.outfile is not None:
        merge_write_df(outdf, args.outfile)
    else:
        print(outdf)
