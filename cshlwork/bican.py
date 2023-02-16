#!/usr/bin/env python
#

import argparse
import logging
import h5py
from scipy import sparse

def read_scqc_h5ad(infile):
  
    with h5py.File(infile) as f:
        logging.debug(list(f.keys()))
        for gk in list(f.keys()):
            g = f[gk]
            for sgk in list(g.keys()):
                ds = g[sgk]
                if type(ds) == h5py._hl.dataset.Dataset:
                    print(f'{gk}/{sgk} shape={ds.shape} type={ds.dtype} ')
                
        
        #['X', 'obs', 'obsm', 'obsp', 'uns', 'var', 'varm']
        X = f['X']
        obs = f['obs']
        obsm = f['obsm']
        obsp = f['obsp']
        uns = f['uns']
        var = f['var']
        varm = f['varm']
        
        
        #dim = f['dim'][()]
        #genes = [g.decode('utf-8') for g in f['dimnames/genes'][()]]
        #cells = [c.decode('utf-8') for c in f['dimnames/cells'][()]]
        #i = f['matrix/i'][()] - 1
        #j = f['matrix/j'][()] - 1
        #x = f['matrix/x'][()]
        #M = sparse.csc_matrix((x, (i,j)) , shape =dim )
        abs = f.get('uns/abstract')
        print(abs)
        


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

   
    parser.add_argument('-o', '--outfile', 
                        action="store", 
                        dest='outfile', 
                        default='genegomatrix.csv',
                        help='Binary matrix file. ')
                    
    parser.add_argument('infiles', 
                        metavar='infiles', 
                        type=str, 
                        help='one or more infiles' ,
                        nargs='*'
                   )
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    for f in args.infiles:
        read_scqc_h5ad(f)