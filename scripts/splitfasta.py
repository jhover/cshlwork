#!/usr/bin/env python
#   
#  Takes .fasta files and splits them, one entry per file, 
#  named after ID (first string after ">"
#

import argparse
import json
import logging
import os
import sys
import traceback


def write_current(linelist):
    header = linelist[0]
    tfaid = header.split()[0][1:]
    tfaid = tfaid.replace('/','-')
    st = "\n".join(linelist)
    with open(f'{tfaid}.fasta', 'w') as f:
        f.write(st)

def handle_infile(infilepath):
    current = None
    processed = 0
    
    with open(infilepath, 'r') as f:
        try:
            while True:
                line = f.readline().strip()
                if line == '':
                        break
                
                if line.startswith(">"):
                    if current is not None:
                        write_current(current)
                        processed += 1
                    current = []
                    current.append(line)

                else:
                    current.append(line)
        except Exception as e:
            traceback.print_exc(file=sys.stdout)

    logging.debug(f'done writing {processed} fasta files. ')
               

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
                        help='fasta input file')

    
    args= parser.parse_args()
   
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    handle_infile(args.infile)