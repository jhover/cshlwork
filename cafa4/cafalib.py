#!/usr/bin/env python
# 
__author__ = "John Hover"
__copyright__ = "2019 John Hover"
__credits__ = []
__license__ = "Apache 2.0"
__version__ = "0.99"
__maintainer__ = "John Hover"
__email__ = "hover@cshl.edu"
__status__ = "Testing"

import argparse
from configparser import ConfigParser
import logging
import os
import sys
import tempfile
from bioservices import uniprot

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

class CAFA4Run(object):
    
    def __init__(self, config, targetlist):
        '''
        Embodies all the processing for a single run against all targets.
        Overall input is a set of Target sequence files. 
        Overall output is a properly-formatted CAFA4 prediction file.   
        
        '''
        self.config = config
        self.targetlist = targetlist
        self.log = logging.getLogger()
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.author = config.get('global','author')
        self.modelnumber = config.get('global','modelnumber')
        self.keywords = config.get('global','keywords')
        self.pipeline = [ x.strip() for x in config.get( 'global','pipeline').split(',')]
    
    def __repr__(self):
        s = "CAFA4Run:"
        for atr in ['outdir','targetlist','pipeline']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s


    def execute(self):
        self.log.info("Begin run...")
        
        phm = Phmmer(self.config, self.targetlist)
        self.log.debug(phm)
        df = phm.execute()
        df.drop(columns=['target','bias','prot_spec'], inplace=True)
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/phmmer.csv" % self.outdir)
        
        ortho = Orthologs(self.config)
        self.log.debug(ortho)
        df = ortho.execute(df)
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/ortho.csv" % self.outdir)
        self.log.info("Ending run...")


class Phmmer(object):
    '''
    Pipeline object. Takes list of Fasta files to run on, returns pandas dataframe. 
    Input:   List of FASTA files 
    Output:  Pandas DataFrame
    
    '''

    def __init__(self, config, targetlist):
        self.log = logging.getLogger()
        self.config = config
        self.targetlist = targetlist
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.database = os.path.expanduser( config.get('phmmer','database') )
        self.score_threshold = config.get('phmmer','score_threshold')
        self.cpus = config.get('phmmer','cpus')
    
    def __repr__(self):
        s = "Phmmer:"
        for atr in ['targetlist','database', 'score_threshold','cpus']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s        
    
    def execute(self):
        outlist = self.run_phmmer_files(self.targetlist, self.database)
        outfile = outlist[0]
        df = self.read_phmmer_table(outfile)
        return df
            
    
    def run_phmmer_files(self, targetlist, database="/data/hover/data/uniprot/uniprot_sprot.fasta"):
    #
    #  time phmmer --tblout 7955.phmmer.2.txt 
    #              --cpu 16 
    #              --noali 
    #              ~/data/cafa4/TargetFiles/sp_species.7955.tfa 
    #              ~/data/uniprot/uniprot_sprot.fasta 
    #              > 7955.phmmer.console.out 2>&1
        
        dbase = database
        outfiles = []
        for file in targetlist:
            (cmd, outfile) = self._make_phmmer_cmdline(file)
            self.log.debug("Running cmd='%s' outfile=%s " % (cmd, outfile))
            cp = subprocess.run(cmd, 
                                shell=True, 
                                universal_newlines=True, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            
            outfiles.append(outfile)
            self.log.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd,outfile, cp.returncode))
        return outfiles
            
    def _make_phmmer_cmdline(self, filename):
        outpath = os.path.dirname(filename)
        filebase = os.path.splitext(os.path.basename(filename))[0]
        outfile = "%s/%s.phmmer.tbl.txt" % (self.outdir, filebase)
        #self.log.debug("outfile=%s" % outfile)
        cmdlist = ['time', 'phmmer']
        cmdlist.append( '--tblout  %s ' % outfile )
        cmdlist.append('--noali' )
        cmdlist.append('--cpu %s ' % self.cpus)
        if self.score_threshold is not None:
            cmdlist.append("-T %s " % self.score_threshold)
        cmdlist.append(' %s ' % filename )
        cmdlist.append(' %s ' % self.database )
        cmd = ' '.join(cmdlist).strip()
        #self.log.debug("command is '%s'" % cmd)
        return (cmd, outfile)

    def read_phmmer_table(self, filename):
        df = pd.read_table(filename, 
                         names=['target','t-acc','query','q-acc',
                                'evalue', 'score', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                                'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                         skip_blank_lines=True,
                         comment='#',
                         index_col=False,
                         skiprows=3,
                         engine='python', 
                         sep='\s+')
        df = df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 
                 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 
                 'description'], axis=1)
        df['db'] = df.apply(lambda row: row.target.split('|')[0], axis=1)
        df['tacc'] = df.apply(lambda row: row.target.split('|')[1], axis=1)
        df['prot_spec'] = df.apply(lambda row: row.target.split('|')[2], axis=1)
        df['protein'] =   df.apply(lambda row: row.prot_spec.split('_')[0], axis=1)
        df['species'] =   df.apply(lambda row: row.prot_spec.split('_')[1], axis=1)
        return df
    

class Orthologs(object):
    '''
    Pipeline object. Takes Pandas dataframe, looks up orthologs and GO values by 
    Input:  Pandas DataFrame
    Output: Pandas DataFrame
    '''

    def __init__(self, config):
        '''
        
        '''
        self.log = logging.getLogger()
        self.config = config
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.uniprot = uniprot.UniProt()

    def __repr__(self):
        s = "Orthologs: uniprot"
        for atr in ['outdir']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s            

        
    def execute(self, dataframe):
        '''
        for each row of dataframe, look up ortholog in uniprot and for each GO code
        add a new row with gene, goterm, gocategory
        
        iterate input df fully, putting new info in new df. 
        merge old + new df, return resulting dataframe
        
        '''
        self.log.info("Looking up each ortholog")
        entries = dataframe['tacc'].unique().tolist()
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
            #print("entry is %s" % entry) 
            #print("gene_names are %s" % gene_names) 
            #print("gene is %s" % gene) 
            #print("golist is %s" % golist)     
            #print("newrowdict is %s" % newrowdict) 
        
        godf = pd.DataFrame.from_dict(newrowdict, orient='index', columns=['entry','gene','goterm']) 
        #print(newdf)
        
        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            (query, evalue, score, db, tacc, protein, species) = row[1:]
            gomatch = godf[ godf.entry == tacc ]
            for gr in gomatch.itertuples():
                (entry, gene, goterm) = gr[1:]
                newrow = [query, evalue, score, db, tacc , protein, species, gene, goterm ]
                newdfdict[ix] = newrow
                ix += 1
        newdf = pd.DataFrame.from_dict(newdfdict, orient='index', columns = ['query', 'evalue', 'score', 'db', 'tacc' , 
                                                                             'protein', 'species', 'gene', 'goterm'])
        self.log.debug("\n%s" % str(newdf))
        return newdf
            
    
    def _query_uniprot(self, tacc):
        '''
        query http://www.uniprot.org/uniprot/<tacc>.txt and parse for gene, goterm, gocategory. 
        
        u.search("id:P35213")
        
        get_df(self, entries, nChunk=100, organism=None):
            entries=['id:P35213'] 
            df = u.get_df(entries)
        
        columns:
            'Gene ontology IDs'
            'Gene ontology (GO)'
            'Gene ontology (biological process)'
            'Gene ontology (molecular function)'
            'Gene ontology (cellular component)'
        
        
        In [45]: df['Gene ontology (biological process)'][0]                                                                                            
        Out[45]: 'cytoplasmic sequestering of protein [GO:0051220]; 
                 negative regulation of G protein-coupled receptor signaling pathway [GO:0045744]; 
                 negative regulation of protein dephosphorylation [GO:0035308]; 
                 negative regulation of transcription, DNA-templated [GO:0045892]; 
                 positive regulation of catalytic activity [GO:0043085]; 
                 protein heterooligomerization [GO:0051291]; protein targeting [GO:0006605]'    
        
        In [47]: df['Gene names  (primary )'][0]                                                                                                        
        Out[47]: 'Ywhab'

        In [48]: df['Gene names'][0]                                                                                                                    
        Out[48]: ['Ywhab']
        
        retrieve(self, uniprot_id, frmt="xml", database="uniprot"):   'txt' good too
        
         
        '''    
        pass



def get_plugin(klassname):
    return getattr(sys.modules[__name__], klassname)
    

       

   



   

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s')
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('infiles', 
                        metavar='infiles', 
                        type=str, 
                        nargs='+',
                        help='a list of .fasta sequence files')
    
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
           
    c4run = CAFA4Run(cp, args.infiles)
    logging.debug(c4run)
    c4run.execute()
    
    
    
    
 

    
    
    