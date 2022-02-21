#!/usr/bin/env python
#
# cococonet related functions
# convert identifiers using gene_info.csv files from cococonet
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

INFO2=   ["GeneSymbol","EntrezID","LocusTag","Synonyms","dbXrefs","Chromosome","Type","EnsemblID","NA","OrthoID","Synon","UniProtID","NetworkIDs"]

COL2IDX = {'gene' :  0,
           'entrez' : 1, 
           'locus' :  2,
           'ensembl': 7,
           'uniprot' : 8
           }

def process_infile(filepath, map):
    """
    map = 
      {  fromcol  : tocol,
         from   : tocol
      }    
    """
    with open(filepath, 'r') as f:
        outdoc = ''
        for line in f:
            linelist = []
            line = line.strip()
            if ',' in line:
                sep = ','
                fields = line.split(',')
            elif '\t' in line:
                sep = '\t'
                fields = line.split(',')
            else:
                sep = ' '
                fields = line.split()
            
            for fld in fields:
                newval = fld
                try:
                    newval = map[fld]
                except KeyError:
                    pass
                linelist.append(newval)
            outline = sep.join(linelist)
            outdoc = outdoc + f'{outline}\n'
    return outdoc            
            
            
    #logging.debug(fields)



def process_infile_old(filepath, map):
    """
    map = 
      {  fromcol  : tocol,
         from   : tocol
      }    
    """
    pairlist = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if ',' in line:
                fields = line.split(',')
            else:
                fields = line.split()
            if fields[0] == fields[1]:
                logging.error(f'fields in pair file the same {fields[0]}')
            else:    
                u1 = get_uniprot(fields[0], map)
                u2 = get_uniprot(fields[1], map)
                if u1 is not None and u2 is not None:
                    print(f'{u1},{u2}')
                    pairlist.append(u1,u2)
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


def build_map(filepath, fromcol, tocol):
    map = {}
    logging.debug(f'opening {filepath}')
    f = open(filepath, 'r')
    linezero = next(f)
    logging.debug(f'{linezero}')
    for line in f:
        line = line.strip()
        fields = line.split(",")
        fromidx = COL2IDX[fromcol]
        toidx = COL2IDX[tocol]
        fromvalue = fields[fromidx].replace('"','')
        tovalue = fields[toidx].replace('"','')
        #gene = fields[0].replace('"','')
        #ensembl = fields[1].replace('"','')
        #locus = fields[2].replace('"','')
        #uniprot1 = fields[8].replace('"','')
        #uniprot2 = fields[9].replace('"','')                    
        #logging.debug(f"ensembl={ensembl} locus={locus} uniprot1={uniprot1} uniprot2={uniprot2}")
        #map[ensembl] = (uniprot1, uniprot2)
        map[fromvalue] = tovalue
    f.close()
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

    parser.add_argument('-f', '--fromid', 
                        metavar='fromid',
                        default='ensembl',
                        type=str,
                        help='identifier from [gene, entrez, locus, ensembl, uniprot]')

    parser.add_argument('-t', '--toid', 
                        metavar='toid',
                        default='uniprot',
                        type=str,
                        help='identifier to [gene, entrez, locus, ensembl, uniprot]')

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

    map = build_map(args.infofile, args.fromid, args.toid)
    #for k in map.keys():
    #    val = map[k]
    #    print(f'{k} -> {val}')
    outtext = process_infile(args.infile, map)
    print(outtext)
