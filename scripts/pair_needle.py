#!/usr/bin/env python
#
#  Split file by lines into N pieces.
#  splitfile.py  <file> <N>

import argparse
import os
import sys
import logging
import traceback
import subprocess


gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

def do_needle(infile, outfile):
    logging.debug(f"processing {infile} to {outfile}...")
    o = open(outfile, 'w')
    with open(infile) as f:
        for i, l in enumerate(f):
            p1, p2 = l.split('\t')
            p1 = p1.strip()
            p2 = p2.strip()
            logging.debug(f"p1={p1} p2={p2}")
            run_needle(p1, p2, o)
    f.close()
    o.close()
    
    
def run_needle(p1, p2, outf):
    
    #p1 = f"{p1}_HUMAN"
    #p2 = f"{p2}_HUMAN"    
    cmd = f'needle -brief -gapopen 10.0 -gapextend 0.5 -stdout -auto -aformat srspair uph:{p1} uph:{p2}'
    cmdlist = cmd.split()
    logging.debug(f"command is {cmd}")
    try:
        p = subprocess.run(cmdlist, check=True, stdout=subprocess.PIPE, universal_newlines=True)    
        output = p.stdout
        lines = output.split('\n')
        towrite = f'{lines[0]}\n'
        outf.write(towrite)
        logging.debug(f"wrote: '{towrite}'")
    except subprocess.CalledProcessError:
        logging.warning(f"Problem with p1={p1} p2={p2}")
        



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

    parser.add_argument('outfile', 
                        metavar='outile', 
                        type=str, 
                        help='pairwise info. <p1> <p2> <len> <dist>')
    
    args= parser.parse_args()

   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    do_needle(args.infile, args.outfile)
    
    