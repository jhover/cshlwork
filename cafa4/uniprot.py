#!/usr/bin/env python
#
# Encapsulates bioservices UniProt API for use with CAFA
#
import configparser
import logging
import os

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

from bioservices.uniprot import UniProt

class UniProtGOPlugin(object):
    '''
    Aux info plugin.
    Takes dataframe, extracts entry_ids, adds info from uniprot.  
    Returns modified dataframe. 
    '''

    def __init__(self, config):
        self.log = logging.getLogger()
        self.uniprot = UniProt()
        self.outdir = os.path.expanduser( config.get('global','outdir') )


    def get_df(self, dataframe):
        entries = dataframe['proteinid'].unique().tolist()
        self.log.debug("Querying uniprot for %d unique entries" % len(entries))
        udf = self.uniprot.get_df(entries)
        self.log.debug(udf)
        udf.to_csv("%s/uniprot.csv" % self.outdir)
        udfslim = udf[['Entry','Gene names','Gene ontology IDs']] 
        # df.tacc corresponds to udf.Entry  ...
        
        newrowdict = {} 
        ix = 0 
        for row in udfslim.itertuples():
            (entry, gene_names, golist) = row[1:] 
            #print("gene_names is %s" % gene_names) 
            try: 
                namestr = gene_names[0] 
                gene = namestr.split()[0] 
            except: 
                gene = '' 
            for goterm in golist: 
                #print("creating new row: %s : %s %s %s" % (ix, entry, gene, goterm)) 
                newrow = [ entry, gene, goterm ] 
                newrowdict[ix] = newrow 
                ix += 1       
        
        godf = pd.DataFrame.from_dict(newrowdict, orient='index', columns=['entry','gene','goterm']) 

        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            self.log.debug("inbound row = %s" % str(row))
            #(query, evalue, score, bias, db, tacc, protein, species) = row[1:]
            (cafaid, evalue, score, bias, db, proteinid, protein, species, cafaprot, cafaspec) = row[1:]
            gomatch = godf[ godf.entry == proteinid ]
            for gr in gomatch.itertuples():
                (entry, gene, goterm) = gr[1:]
                newrow = [cafaid, evalue, score, bias, db, proteinid , protein, species, cafaprot, cafaspec, gene, goterm ]
                newdfdict[ix] = newrow
                ix += 1
        newdf = pd.DataFrame.from_dict(newdfdict, orient='index', columns = ['cafaid', 'evalue', 
                                                                             'score', 'bias', 
                                                                             'db', 'proteinid' , 
                                                                             'protein', 'species', 
                                                                             'cafaprot', 'cafaspec',
                                                                             'gene', 'goterm'])
        self.log.debug("\n%s" % str(newdf))        
        return newdf


if __name__ == '__main__':
    
    c = configparser.ConfigParser()
    upg = UniProtGOPlugin(c)
    #entrylist = ['Q9CQV8', 'P35213', 'A4K2U9', 'P31946', 'Q4R572', 'P68250']
    #out = upg._query_entries(entrylist)
    #print(out)    

    
    