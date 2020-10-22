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
    cmd = f'needle -brief -gapopen 10.0 -gapextend 0.5 -stdout -auto uph:{p1} uph:{p2}'
    cmdlist = cmd.split()
    logging.debug(f"command is {cmd}")
    try:
        p = subprocess.run(cmdlist, check=True, stdout=subprocess.PIPE, universal_newlines=True)    
        output = p.stdout
        lines = output.split('\n')
        towrite = parse_output(lines)
        if towrite is not None:
            outf.write(towrite)
            logging.debug(f"wrote: '{towrite}'")
        else:
            logging.warning(f"problem parsing output with  p1={p1} p2={p2} ")

    except subprocess.CalledProcessError:
        logging.warning(f"Problem with p1={p1} p2={p2}")
        towrite = f"{p1}\t{p2}\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\n"    
      
    except Exception:
        logging.warning("Some other problem.")



def parse_output(lines, p1, p2):
    '''
    assumes default srspair aformat
    
    out | head -33 | tail 17:
    #=======================================
    #
    # Aligned_sequences: 2
    # 1: GP157_HUMAN
    # 2: SPSB1_HUMAN
    # Matrix: EBLOSUM62
    # Gap_penalty: 10.0
    # Extend_penalty: 0.5
    #
    # Length: 564
    # Identity:      16/564 ( 2.8%)
    # Similarity:    21/564 ( 3.7%)
    # Gaps:         520/564 (92.2%)
    # Score: 22.0
    # 
    #
    #=======================================

    '''
    try:
        lines = lines[15:32]
        p1 = lines[3].split()[2]
        p2 = lines[4].split()[2]            
        length = int(lines[9].split()[2])
        ident = int(lines[10].split()[2].split('/')[0] )
        simil = int(lines[11].split()[2].split('/')[0] )
        gaps = int(lines[12].split()[2].split('/')[0] )
        score = float(lines[13].split()[2])
        pident = ident / length
        psimil = simil / length
        out = f"{p1}\t{p2}\t{length}\t{ident}\t{simil}\t{gaps}\t{score}\t{pident:.3f}\t{psimil:.3f}\n"
        logging.debug(out)
        return out
    
    except Exception:
        logging.debug(f"Problem parsing output for {p1} {p2}")
        return None



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
                        metavar='outfile', 
                        type=str, 
                        help='pairwise info. <p1> <p2> <len> <dist>')
    
    args= parser.parse_args()

   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    do_needle(args.infile, args.outfile)
    
