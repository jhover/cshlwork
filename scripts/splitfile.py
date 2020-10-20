#!/usr/bin/env python
#
#  Split file by lines into N pieces.
#  splitfile.py  <file> <N>

import argparse
import os
import sys
import logging
import traceback

gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

def file_len(filename):
    i = -1
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1

#
#  i=0  line
#  i=1 
#
def do_split(filename, nfiles):
    logging.debug(f"splitting {filename} into {nfiles} pieces...")
    nlines = file_len(filename)
    logging.debug(f"{filename} has {nlines} lines.")
    lpf = int(nlines / nfiles)
    logging.debug(f"will put {lpf} lines per file...")
    fnum = 0
    of = None
    with open(filename) as f:
        for i, l in enumerate(f):
            if i % lpf == 0:
                if of is not None:
                    logging.debug("closing file...")
                    of.close()
                    fnum += 1
                logging.debug(f"opening new file= {filename}.{fnum} ")
                of = open(f"{filename}.{fnum}", 'w')
            of.write(l)
    of.close()
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
                        help='a .fasta sequence file')

    parser.add_argument('nfiles', 
                        metavar='nfiles', 
                        type=int, 
                        help='Number of pieces.')
    
    args= parser.parse_args()

   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    do_split(args.infile, args.nfiles)