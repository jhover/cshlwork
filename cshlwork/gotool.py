#!/usr/bin/env python
#
# Prints one or more specified GO terms for OBO file to stdout. 
# usage gocat.py  <GOTERM>

#      

import argparse
from configparser import ConfigParser
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.ontology import GeneOntology, GOMatrix

OBOFILE=os.path.expanduser('~/data/go/go.obo') 
GAFFILE=os.path.expanduser('~/data/go/zfin.gaf')
CONFFILE=os.path.expanduser('~/git/cshl-work/etc/gotool.conf')


def gaf_to_df():
    """

db    dbobji              dbobjsym     goterm       dbref             goevidence  withfrom        goaspect  dbobjnam  dbobjtype taxonid  date  assignedby   annotext   geneprodid  
    
ZFIN  ZDB-GENE-070410-141  zgc:163098  GO:0003723  ZFIN:ZDB-PUB-170525-1  IEA  UniRule:UR000414619  F zgc:163098 protein  taxon:7955  20190914  UniProt  UniProtKB:A0A2R8QLY
    """
    pass




def gocat(gotermlist):
    filehandle = open(GOFILE, 'r')
    lines = filehandle.readlines()
    print("read in %d lines" % len(lines))
    print("got arg %s" % gotermlist)
    for gt in gotermlist:
        found = False
        for line in lines:
            if line.startswith("id: %s" % gt):
                found = True
                print('[Term]')
                #print(line.strip())
            if line.startswith("[Term]"):
                found = False
            if found and not line.startswith("[Term]"):
                print(line.strip()) 
    

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

    parser.add_argument('-b', '--obofile', 
                        action="store", 
                        dest='obofile', 
                        default=OBOFILE,
                        help='Gene ontology OBO file.')
    
    parser.add_argument('-g', '--gaffile', 
                        action="store", 
                        dest='gaffile', 
                        default=GAFFILE,
                        help='GAF file. ')

    parser.add_argument('-o', '--outfile', 
                        action="store", 
                        dest='outfile', 
                        default='genegomatrix.csv',
                        help='Binary matrix file. ')

   
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default=CONFFILE,
                        help='Config file path [~/etc/gotool.conf]') 
                    
    parser.add_argument('goterms', 
                        metavar='goterms', 
                        type=str, 
                        help='one or more space-separated goterms GO:XXXXXXXX' ,
                        nargs='*'
                   )
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)

    #gocat(args.goterms)
    
    
    
    