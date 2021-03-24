#!/usr/bin/env python
#   
#  splits any file into pieces by given number of lines
#

import argparse
import logging
import os


def splitfile(filepath, nlines):
    filepath = os.path.expanduser(filepath)
    (fbase, fext) = os.path.splitext(filepath)
     

    f = open(filepath)
    lines = f.readlines()
    tlines = len(lines)
    logging.debug(f"read file with {len(lines)} lines.")
    
    fileno=0
    baseline = 0
    endline = nlines
    while baseline < tlines:
        wlines = lines[baseline:endline]
        fn = f"{fbase}_{fileno}{fext}"
        logging.info(f"Using new file {fn}")
        of = open(fn, 'w')
        of.writelines(wlines)
        of.close()
        baseline += nlines
        endline += nlines
        fileno += 1
    f.close()

if __name__=='__main__':

    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARNING)

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
                        help='')

    parser.add_argument('nlines', 
                        metavar='nlines', 
                        type=int, 
                        help='number of lines per file. ')    
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    splitfile(args.infile, args.nlines )
    