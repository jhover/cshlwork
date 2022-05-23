#!/usr/bin/env python
# 
# takes input file and emits same file but with identifiers translated to target.
# infile formats can be comma, space, or tab-separated
#

import argparse
import logging
import os
import pickle
import requests
import sys
import time
import traceback

from collections import defaultdict
from configparser import ConfigParser

import pandas as pd

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from protlib.utils import *
from protlib.uniprot import *

ASPECTMAP = { 'C': 'cc',
              'F': 'mf',
              'P': 'bp'
            }

KEYMAP =  { 'OrderedLocusNames' : 'locus'  ,
            'Name' : 'gene',
            'ORFNames' : 'orfname',
            'Synonyms' : 'synonym' 
            }

VALID_IDS = ['accession', 'proteinid','locus', 'gene']


def process_file(infile, ididx, outid):
    with open(infile) as f:
        for line in f:
            fields = line.split('\t')
            newline = []
            for i, field in enumerate(fields):
                field = field.strip()
                logging.debug(f'i={i}   field={field}')
                try:
                    edict = ididx[field]
                    newid = edict[outid]
                    newline.append(newid)
                    logging.debug(f'found {newid} to replace {field}')
                except KeyError:
                    logging.debug(f'no entry found for {field}')
                    newline.append(field)
            print('\t'.join(newline))



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
    
    parser.add_argument('-u', '--uniprot', 
                        metavar='uniprot', 
                        type=str,
                        nargs='?',
                        default=os.path.expanduser("~/data/uniprot/uniprot_all_vertebrates.dat"), 
                        help='A Uniprot .dat file') 

    parser.add_argument('-i','--inid', 
                        metavar='inid', 
                        type=str,
                        required=False,
                        default='accession', 
                        help='Identifier type [ accession | proteinid | locus ] (B4FCB1 | B4FCB1_MAIZE ')    

    parser.add_argument('-o','--outid', 
                        metavar='outid', 
                        type=str,
                        required=False,
                        default='accession', 
                        help='Identifier type [ accession | proteinid | locus ] (B4FCB1 | B4FCB1_MAIZE ')  
    
    parser.add_argument('infile' ,
                        metavar='infile', 
                        type=str,
                        nargs='?',
                        default=None, 
                        help='UPID [UPID UPID ...] ')
    
        
    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    

    # Prepare for commands that require parsing DAT file. 
    logging.info(f"Parsing uniprot .dat file={args.uniprot} ...")

    config = get_default_config()
    
    entries = parse_uniprot_dat(config, datfile=args.uniprot)
    
    idx = index_by(entries, args.inid)
    
    #(entries, pididx, accidx) 
    logging.debug(f"uniprot length = {len(entries)}") 
       
    logging.debug(f'infile={args.infile}')
    process_file(args.infile, idx, args.outid)

