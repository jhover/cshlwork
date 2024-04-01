#!/usr/bin/env python
#
#   get mappings between transcripts, genes, and gene symbol/name. 
#
import os
import sys

gitpath=os.path.expanduser("~/git/cellreader")
sys.path.append(gitpath)

import logging
from pybiomart import Server
from cellreadr.caching import *
 

DATASETS = { 'Mus_musculus': 'mmusculus_gene_ensembl',
             'Homo_sapiens': 'hsapiens_gene_ensembl',
             'Rattus_norvegicus': 'rnorvegicus_gene_ensembl',
             'Macaca_mulatta' : 'mmulatta_gene_ensembl',
             'Macaca_fascicularis' : 'mfascicularis_gene_ensembl',
             'Callithrix_jacchus' : 'cjacchus_gene_ensembl',
             'Taeniopygia_guttata': 'tguttata_gene_ensembl',
    }


def get_mapping_df(species='Mus_musculus', cache=True):
    
    atts = ['external_gene_name',
            #'external_gene_source',
            'ensembl_gene_id',
            'ensembl_transcript_id',
            #'ensembl_peptide_id'
            ]

    
    server = Server(host="http://www.ensembl.org")
    dataset_id = DATASETS[species]
    ege = server.marts['ENSEMBL_MART_ENSEMBL']
    dataset = ege.datasets[dataset_id]
   
    df = dataset.query(attributes = atts)
    df.columns=atts
    return df

def get_mapping_by_gene(species='Mus_musculus', cache=True):
    key = f'tid_geneid_map.{species}'
    if cache_object_exists(key):
        tid_geneid_map = get_cache_object(key)
    else:
        tid_geneid_map = {}
        df = get_mapping_df(species=species)
        for index, row in df.iterrows():
            tid_geneid_map[ row['ensembl_transcript_id'] ] = row['ensembl_gene_id']
        logging.debug(f'mapped {len(tid_geneid_map)} transcripts from {df.shape[0]} rows')
    
    store_cache_object(tid_geneid_map, key)
    return tid_geneid_map
    
    

if __name__ == '__main__':
    df = get_mapping_df()
    print(df)