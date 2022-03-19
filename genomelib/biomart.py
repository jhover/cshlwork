#!/usr/bin/env python
#
#   get mappings between transcripts, genes, and gene symbol/name. 
#
import logging
from pybiomart import Server


DATASETS = { 'Mus_musculus': 'mmusculus_gene_ensembl',
             'Homo_sapiens': 'hsapiens_gene_ensembl',
    }

def get_mapping_df(species='Mus_musculus'):
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

def get_mapping_by_gene(species='Mus_musculus'):
    tid_geneid_map = {}
    
    df = get_mapping_df(species=species)
    for index, row in df.iterrows():
        tid_geneid_map[ row['ensembl_transcript_id'] ] = row['ensembl_gene_id']
    logging.debug(f'mapped {len(tid_geneid_map)} transcripts from {df.shape[0]} rows')
    return tid_geneid_map
    
        



if __name__ == '__main__':
    df = get_mapping_df()
    print(df)