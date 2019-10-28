#!/usr/bin/env python
#
# read in starout format file. 
# create Series
# merge Series to dataframe
#

#  STAR / Output format. 
# -quantMode GeneCounts 
# -hseq-count -s yes   /  -s reverse
# columns
# geneID       unstranded_counts   1st_stranded  2nd_stranded

# N_unmapped       2674617    2674617    2674617
# N_multimapping   2828822    2828822    2828822
# N_noFeature      5598327    24570638    24512243
# N_ambiguous      3871509    813790    812892
# ENSG00000223972.5    0    0    0
# ENSG00000227232.5    302    159    143
# ENSG00000278267.1    9    7    2
# ENSG00000243485.4    1    1    0
# ENSG00000237613.2    0    0    0
#

import argparse
import logging
import os

import pandas as pd
import numpy as np


def processall(filelist):
    ## Use Series
    dslist = []
    for f in filelist:
        ds = parsefile2series(f)
        dsr = ds.rank()
        logging.info("\n%s" % dsr)
        dslist.append(dsr)
  
    df = pd.concat(dslist, axis=1)
    
    # Pairwise correlation 
    cdf = df.corr(method='spearman')
    print(cdf)

    return cdf
      
 
def parsefile2series(filename):
    logging.info("Processing file %s" % filename)

    (head, tail) = os.path.split(filename)    
    sname = tail.split('.')[0]
      
    f = open(filename)
    lines = f.readlines()
    data = {}
          
    # handle values. 
    for line in lines[4:]:
        (label, unstrand, strand1, strand2) = [ f.strip() for f in line.split('\t') ]
        data[label] = strand2

    logging.info("data length is %s" % len(data))
    ds = pd.Series( data, name=sname )
    logging.info("\n%s" % ds.head() )  
    return ds    



if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s')
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('infiles', 
                        metavar='infiles', 
                        type=str, 
                        nargs='+',
                        help='a list of single-cell expression files')
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    filelist = args.infiles 
    logging.debug("%d files to process. " % len(filelist))
    
    df = processall(filelist)
    

    
    
    