#! /usr/bin/env python
import argparse
import logging
import os
import sys
import pandas as pd

'''
  Take in trembl   fasta file and create Pandas-compatible tsv of accession, uniprotid, genename

'''

def parse_uniprot_fasta(filename):
    f = open(filename, 'r')
    lod = []
    nlines = 0
    
    for line in f:
        nlines += 1
        if nlines % 10000 == 0:
            logging.debug(f"processed {nlines} lines...")
        
        if line.startswith('>'):
            d = {}
            # >tr|A0A1S3RID5|A0A1S3RID5_SALSA ras and Rab interactor 2-like OS=Salmo salar OX=8030 GN=LOC106602976 PE=4 SV=1
            fields = line.split()
            (db, acc, uid) = fields[0][1:].split('|')
            d['db'] = db
            d['acc'] = acc
            d['uid'] = uid
            
            #logging.debug(f"db={db} acc={acc} uid={uid}")
            for t in fields[1:]:
                if '=' in t:
                    (key,val) = t.split('=')
                    #logging.debug(f"key={key} val={val}")
                    if key == 'GN':
                        d['gn'] = val
            lod.append(d)
    
    logging.info(f"Processed {nlines} lines. Done.")
     
    df = pd.DataFrame(lod) 
    return df
           
    


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence file')

    parser.add_argument('outfile', 
                        metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')
    
    args= parser.parse_args()
    
    logging.debug(f"infile={args.infile} outfile={args.outfile}")
    tdf = parse_uniprot_fasta(args.infile)
    logging.debug(f"tdf=\n{tdf}")
    tdf.to_csv(args.outfile, sep='\t' )
    logging.debug(f"Wrote to {args.outfile}")

    
