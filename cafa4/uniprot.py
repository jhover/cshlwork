#!/usr/bin/env python
#
# Encapsulates bioservices UniProt API and local files for use with CAFA
#

import argparse
from configparser import ConfigParser
import logging
import os

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

from bioservices.uniprot import UniProt   # to query online REST interface
from Bio import SeqIO   # to parse uniprot.dat

#
#  Possibly compare  https://pypi.org/project/PyUniProt/
#

class UniProtGOPlugin(object):
    '''
    Aux info plugin.
    Takes dataframe, extracts entry_ids, adds info from uniprot.  
    Returns modified dataframe. 
    '''

    def __init__(self, config):
        self.log = logging.getLogger('UniProtGOPlugin')
        self.uniprot = UniProt()
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.sprotdatfile = os.path.expanduser( config.get('goplugin','sprotdatfile') )


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
        newdf = pd.DataFrame.from_dict(newdfdict, orient='index', columns = ['cafaid', 
                                                                             'evalue', 
                                                                             'score', 
                                                                             'bias', 
                                                                             'db', 
                                                                             'proteinid' , 
                                                                             'protein', 
                                                                             'species', 
                                                                             'cafaprot', 
                                                                             'cafaspec',
                                                                             'gene', 
                                                                             'goterm'])
        self.log.debug("\n%s" % str(newdf))        
        return newdf

    def _open_swissprot_file(self):
        '''
         Read uniprot_sprot.dat and get dataframe of relevant fields.

        '''
        self.log.debug("Getting swissprot DF")
        filehandle = None
        try:
            self.log.debug("opening file %s" % self.sprotdatfile )
            filehandle = open(self.sprotdatfile, 'r')
            self._parsefile(filehandle)
            filehandle.close()
                
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)        
        
        finally:
            if filehandle is not None:
                filehandle.close()
        out = self._parsefile(filehandle)
        self.log.debug("Parsed data file.")


    def get_swissprot_df(self):
        from Bio import SeqIO
        rgen = SeqIO.parse(self.sprotdatfile,"swiss")
        for record in rgen:
            pass


    def _parsefile(self, filehandle):
        '''
    
        '''
        current = None
        try:
            for line in filehandle:
                if line.startswith("ID "):
                    # ID   001R_FRG3G              Reviewed;         256 AA.
                    val = line[5:]
                    self.log.debug("Beginning of entry.")                
                elif line.startswith("AC "):
                    #AC   Q6GZX4;
                    val = line[5:]
                elif line.startswith("DR "):
                    # DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
                    val = line[5:]              
                elif line.startswith("// "):
                    self.log.debug("End of entry.")
                
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        self.log.debug("Parsed file with %d terms" % keyid )
  
        

def test_uniprot(config):
    upg = UniProtGOPlugin(config)
    entrylist = ['Q9CQV8', 'P35213', 'A4K2U9', 'P31946', 'Q4R572', 'P68250']
    out = upg._query_entries(entrylist)
    print(out)     

def test_datparse(config):
    upg = UniProtGOPlugin(config)
    upg.get_swissprot_df()
    


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
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default='~/etc/cafa4.conf',
                        help='Config file path [~/etc/cafa4.conf]')
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)

    # test_uniprot(config)
    test_datparse(cp)
  

    
    