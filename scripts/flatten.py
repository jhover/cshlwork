#!/usr/bin/env python
#   
#      takes directory hierachy and puts all files into one outdir
#     creates a mapfile of filenames (without extension) so that other output 
#     can be restored to original hierarchy
#    takes files (ignoring extension) in indir and restores (copies) them to 
#    hierarchy as stored in mapfile, based at rootdir. 
#

import argparse
import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import flatten_tree, unflatten_tree

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

    parser.add_argument('-u', '--unflatten', 
                        action="store_true", 
                        dest='unflatten',
                        default=False, 
                        help='perform unflattening')
     
    parser.add_argument('-o','--outdir', 
                        metavar='outdir',
                        required=True,
                        type=str,
                        help='top output dir ')

    parser.add_argument('-i','--indir', 
                        metavar='indir',
                        required=True,
                        type=str,
                        help='input dir. one file per infile ')

    parser.add_argument('-m','--mapfile', 
                        metavar='mapfile',
                        required=True,
                        type=str,
                        help='file to map trees/filenames')

    parser.add_argument('-e','--extension', 
                        metavar='extension',
                        required=False,
                        type=str,
                        default='*',
                        help='restrict extension on unflatten')

    args= parser.parse_args()
    
    if args.debug:  
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
        
    if args.unflatten:
        unflatten_tree(indir=args.indir, rootdir=args.outdir, mapfile=args.mapfile, ext=args.extension)
    else:
        flatten_tree(indir=args.indir, outdir=args.outdir, mapfile=args.mapfile)  