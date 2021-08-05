#!/usr/bin/env python
#
#  Split a uniprot/swissprot .dat file by entries into N pieces.
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

def entry_len(filename):
    i = 0
    entries = 0
    with open(filename) as f:
        for i, l in enumerate(f):
            if l.startswith('//'):
                entries += 1
    f.close()
    logging.debug(f'handled {i} lines...')
    return entries

#
def do_split(filename, nfiles):
    logging.debug(f"splitting {filename} into {nfiles} pieces...")
    nentries = entry_len(filename)
    logging.debug(f"{filename} has {nentries} entries.")
    epf = int(nentries / nfiles)
    logging.debug(f"will put {epf} entries per file...")
    
    (base, ext) = os.path.splitext(filename)
    logging.debug(f'base={base} ext={ext}')
    
    fnum = 1
    wentries = 0
    of = None
    with open(filename) as f:
        ofname = f"{base}.{fnum}{ext}"
        logging.debug(f"opening new file= {ofname} ")
        of = open(ofname, 'w')
        for i, l in enumerate(f):
            if l.startswith('//'):
                of.write(l)
                wentries += 1
                if wentries >= epf:
                    logging.debug(f'reached {wentries} written to {ofname}')
                    of.close()
                    fnum += 1
                    wentries = 0
                    ofname = f"{base}.{fnum}{ext}"
                    logging.debug(f"opening new file= {ofname} ")
                    of = open(ofname, 'w')    
            else:              
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
                        help='a .dat uniprot/swissprot file')

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