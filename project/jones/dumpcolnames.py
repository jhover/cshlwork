#!/usr/bin/env python


import argparse
import logging

import pandas as pd
import h5py


def parse_expression_hd5(filename):
    """
    Loads data in file to dataframe.
    """
    with h5py.File(filename, 'r') as f:
        logging.debug("reading matrix...")
        matrix = f['agg'][:]
        logging.debug("reading rows. converting to unicode.")
        rows = [ s.decode() for s in  f['row'][:] ]
        logging.debug("reading columns. converting to unicode")
        columns = [ s.decode() for s in  f['col'][:] ]
        logging.debug("making dataframe...")
    return columns
#        df = pd.DataFrame(matrix,  index=rows, columns = columns )
#    return df    

def write_list(itemlist, outfile):
    with open(outfile, 'w') as f:
        for item in itemlist:
            f.write(item)
            f.write('\n')
    f.close()


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

    parser.add_argument('outfile', 
                        metavar='outfile', 
                        type=str, 
                        help='outfile')
   
    args= parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    columns = parse_expression_hd5(args.infile)
    columns.sort()
    write_list(columns, args.outfile)
        