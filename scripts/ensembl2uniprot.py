#!/usr/bin/env python
#
# consume input file and replace ensembl IDs with uniprot ids where available 
# 
#
#
import argparse
import os
import sys
import logging
import traceback

INFOCOLS=['EntrezID','EnsemblID','GeneSymbol','LocusTag','Synonyms','dbXrefs','Chromosome','Type',
          'UniProtID','UniProtID2','OrthoID','NetworkIDs']

def process_infile(filepath, map):
    """
    map = 
      { ensembl : ( uniprot1, uniprot2),
        ensembl : ( uniprot1, uniprot2),
      }    
    """
    pairlist = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            fields = line.split(',')
            if fields[0] == fields[1]:
                logging.error(f'fields in pair file the same {fields[0]}')
            else:    
                u1 = get_uniprot(fields[0], map)
                u2 = get_uniprot(fields[1], map)
                if u1 is not None and u2 is not None:
                    print(f'{u1},{u2}')
                else:
                    logging.error(f'some field None')
                #logging.debug(fields)
            

def get_uniprot(ensemblid, map):
    answer = None
    try:
        u1, u2 = map[ensemblid]
        if u1 != 'NA':
            answer = u1
        elif u2 != 'NA':
            answer = u2
        else:
            logging.error(f'no uniprotid for {ensemblid}')
            
    except:
        logging.error(f'no ensemblid: {ensemblid}')
    return answer


def build_map(filepath):
    map = {}
    with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    fields = line.split(",")
                    ensembl = fields[1].replace('"','')
                    uniprot1 = fields[8].replace('"','')
                    uniprot2 = fields[9].replace('"','')                    
                    logging.debug(f"ensembl={ensembl} uniprot1={uniprot1} uniprot2={uniprot2}")
                    map[ensembl] = (uniprot1, uniprot2)
    return map
                    




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
                        help='Any space or comma-separated file with ensemble ids')

    parser.add_argument('infofile', 
                        metavar='infofile', 
                        default=os.path.expanduser('~/data/cococonet/geneinfo/maize_info.csv'),
                        nargs="?",
                        type=str, 
                        help='A CoCoCoNet gene info file w/ Uniprot IDs.')

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    map = build_map(args.infofile)
    process_infile(args.infile, map)


