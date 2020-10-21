#!/usr/bin/env python
#
# Map Uniprot ACs to IDs
#
#

import argparse
import os
import sys
import logging
import traceback
import subprocess


gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

DATFILE=os.path.expanduser('~/data/uniprot/humandb/uniprot_all_human.dat')


def parse_uniprot_dat(filename):
    filehandle = open(filename, 'r')
    current = None
    allentries = []
    try:
        while True:
            line = filehandle.readline()
            if line == '':
                break

            if line.startswith("ID "):
                # ID   001R_FRG3G              Reviewed;         256 AA.
                #      <prot_name>_<prot_spec>
                proteinid = line[5:16].strip()
                current = defaultdict(dict)
                current['id'] = proteinid
        
            elif line.startswith("AC "):
                # AC   Q6GZX4;
                # AC   A0A023GPJ0; 
                # AC   Q91896; O57469;
                #logging.debug("Handling AC.")
                rest = line[5:]
                accession = rest.split(';')[0]
                #accession = line[5:11].strip()
                current['ac'] = accession

         
            elif line.startswith("//"):           
                allentries.append(current)
                current = None
    except:
        pass



def do_water(infile, outfile):
    pass



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
    
    do_mappings(args.infile, args.outfile)