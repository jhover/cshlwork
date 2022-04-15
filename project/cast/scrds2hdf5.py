#!/usr/bin/env python
#
#  Convert .rds SingleCellExperiment to structure HDF5 file containing a sparse matrix. 
#  
#  Requires that counts matrix be exported first from within R:
#
# require(SingleCellExperiment)
# require(Matrix)
# sce = readRDS('barseq2_whitelisted.rds')
# m = counts(sce)
# writeMM(m, file='mmfile.mtx')
#

import argparse
import os
import sys
import logging

import h5py

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

import numpy as np
import pandas as pd

def readmm(mmfile):
    m = pd.read_csv(mmfile, skiprows=1, sep=' ', names=['i','j','x'] )  
    logging.debug(m)
    return (m.i, m.j, m.x)


def rdstohdf5(infile, mmfile,  outfile):
    logging.debug(f'handling {infile} -> {outfile}')
    importr('SingleCellExperiment')
    importr('Matrix')
    readRDS = robjects.r['readRDS']
    sce = readRDS(infile)
    logging.debug(sce)


    # get column info. 
    coldata = sce.slots['colData']
    colnames = list(coldata.slots['rownames'])
    logging.debug(f'got {len(colnames)} column (cells) names.')
    colnames = np.array(colnames, dtype=np.dtype('S8') )    

    # row data (genes)
    rownames = (list(sce.names))
    logging.debug(f'got {len(rownames)} row (gene) names.')    
    rownames = np.array(rownames, dtype=np.dtype('S8') )
    
    i,j,x = readmm(mmfile)
    #logging.debug(f'{i} {j} {x} ')


    logging.debug(f'opening {outfile} for HDF5...')
    with h5py.File(outfile, "w") as f:
        logging.debug(f'file {outfile} opened for HDF5...')

        dim = np.array([len(rownames), len(colnames) ], dtype='int32' )
        f.create_dataset('dim', data = dim )
        logging.debug(f"created dataset 'dim': {dim} " )
        
        dimnames = f.create_group('dimnames' )
        dimnames.create_dataset('cells', data=colnames, compression="gzip")         
        dimnames.create_dataset('genes', data=rownames, compression="gzip")
        
        matrix = f.create_group('matrix')
        matrix.create_dataset('i', data = i, dtype='int32', compression="gzip")
        matrix.create_dataset('j', data = j, dtype='int32', compression="gzip")
        matrix.create_dataset('x', data = x, dtype='f8', compression="gzip")            

    logging.debug(f'done creating {outfile}')



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
                        help='rds input file')

    parser.add_argument('mmfile', 
                        metavar='mmfile', 
                        type=str, 
                        help='writeMM output of counts matrix in .rds')

    parser.add_argument('outfile', 
                        metavar='outfile', 
                        type=str,
                        default=None,
                        help='hdf5 outfile. ')
    
    args= parser.parse_args()
   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    rdstohdf5(args.infile, args.mmfile,  args.outfile)