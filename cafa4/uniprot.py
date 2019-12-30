#!/usr/bin/env python
#
# Encapsulates bioservices UniProt API and local files for use with CAFA
#

import argparse
from configparser import ConfigParser
import logging
import os
import pprint

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess
import sys
import traceback

from bioservices.uniprot import UniProt   # to query online REST interface
from Bio import SeqIO   # to parse uniprot.dat


#
#  Possibly compare  https://pypi.org/project/PyUniProt/
#

class UniProtRecord(object):
    
    def __init__(self, record):
        '''
        
        ID: Q6GZW9
        Name: 006R_FRG3G
        Description: RecName: Full=Uncharacterized protein 006R;
        Database cross-references: EMBL:AY548484, RefSeq:YP_031584.1, SMR:Q6GZW9, GeneID:2947778, KEGG:vg:2947778, Proteomes:UP000008770
        Number of features: 1
        /accessions=['Q6GZW9']
        
        
        record:
            id:
            name:
            description: 
            annotations: -> dict
            dbxrefs: -> list     
            features: -> SeqFeature

        'letter_annotations', 
        'lower', 
        'name', 'reverse_complement', 'seq', 'translate', 'upper'
        
        E.g 
        annotations -> dict
            'accessions': ['Q6GZW9'], 
            'protein_existence': 4, 
            'date': '28-JUN-2011', 
            'sequence_version': 1, 
            'date_last_sequence_update': '19-JUL-2004', 
            'date_last_annotation_update': '05-JUN-2019', 
            'entry_version': 28, 
            'gene_name': 'ORFNames=FV3-006R;', 
            'organism': 'Frog virus 3 (isolate Goorha) (FV-3)', 
            'taxonomy': ['Viruses', 'Iridoviridae', 'Alphairidovirinae', 'Ranavirus'], 
            'ncbi_taxid': ['654924'], 
            'organism_host': ['Ambystoma (mole salamanders)', 
                               'Dryophytes versicolor (chameleon treefrog)', 
                               'Lithobates pipiens (Northern leopard frog) (Rana pipiens)', 
                               'Notophthalmus viridescens (Eastern newt) (Triturus viridescens)',
                               'Rana sylvatica (Wood frog)'], 
                               'host_ncbi_taxid': ['8295', 
                                                    '30343', 
                                                    '8404', 
                                                    '8316', 
                                                    '45438'], 
                                'references': --> Reference 
        
        dbxrefs -> list
            ['EMBL:AY548484', 'RefSeq:YP_031579.1', 'SwissPalm:Q6GZX4', 'GeneID:2947773', 
            'KEGG:vg:2947773', 'Proteomes:UP000008770', 
            'GO:GO:0046782', 
            'InterPro:IPR007031', 'Pfam:PF04947']
    
        [Reference(title='Comparative genomic analyses of frog virus 3, type species of the genus Ranavirus (family Iridoviridae).',
                     ...)], 
                  'keywords': ['Complete proteome', 
                              'Reference proteome']
        '''
        self.proteinid = record.id
        self.protein = record.name
        self.goterms = []
        for xf in record.dbxrefs:
            if xf.startswith("GO:"):
                gt = xf[3:]
                self.goterms.append(gt)
        self.accessions = record.annotations['accessions']
        self.taxonid = record.annotations['ncbi_taxid'][0]
        #self.species = 
        #if 'accessions' in record.annotations.keys():
        #    self.accessions = record.annotations['accessions']
            
            
            
    def __repr__(self):
        REPR_ATTRS=['proteinid','protein','accessions', 'taxonid','goterms']
        s="UniProtRecord: "
        for atr in REPR_ATTRS:
            s+="%s=%s " % (atr, 
                           self.__getattribute__(atr))
        return s


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
        self.speciesmap = os.path.expanduser( config.get('global','species_mapfile'))


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
        '''
    
        Dataframe:
        proteinid   name   

cf... 
             cafaid         evalue  score  bias  db proteinid protein species cafaprot cafaspec       gene      goterm  probest goaspect
0     T100900000001  1.100000e-156  523.6   8.5  sp    Q9CQV8   1433B   MOUSE    1433B    MOUSE      Ywhab  GO:0005634      1.0       cc
    
        '''
        rlist = self._dat2upr()
        
        for r in rlist:
            for gt in r.goterms:
                pass
                
        #self.id = record.id
        #self.name = record.name
        

            
        
    def _dat2upr(self):
        self.log.debug("opening swissprot dat file %s" % self.sprotdatfile)
        rgen = SeqIO.parse(self.sprotdatfile,"swiss")
        i = 0
        uprlist = []
        self.log.debug("Completed SeqIO.parse(). Handling records...")
        for record in rgen:
            upr = UniProtRecord(record)
            uprlist.append(upr)
            #print(record)
            i += 1
            if i % 10000 == 0:
                self.log.debug("Handled %d records..." % i)
            #    break
        self.log.debug("parsed dat file of %d records" % len(uprlist))
        return uprlist

    def get_annotation_df(self):
        self.log.debug("opening swissprot dat file %s" % self.sprotdatfile)
        rgen = SeqIO.parse(self.sprotdatfile,"swiss")
        self.log.debug("rgen type is %s" % type(rgen))
        #self.log.debug("Created generator with %d records" % len(rgen))
        i = 0
        alltuples = []
        for record in rgen:
            #print(record)
            i += 1
            if i % 10000 == 0:
                self.log.debug("Handled %d records..." % i)
            goterms = []
            for xf in record.dbxrefs:
                if xf.startswith("GO:"):
                    gt = xf[3:]
                    goterms.append(gt)
            if len(goterms) > 0:
                proteinid = record.id
                protein = record.name
                taxonid = record.annotations['ncbi_taxid'][0]                
                for gt in goterms:
                    t = (taxonid, proteinid, protein, gt)
                    alltuples.append(t)
                # fan out over goterms
            else:
                # ignore un-annotated entries. 
                pass
                        
            if i >= 1000:
                break
        self.log.debug("generated %d tuples" % len(alltuples)) 
        df = pd.DataFrame(alltuples, columns=['taxonid','proteinid','protein','goterm'])
        
        return df


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

    def _make_species_map(self):
        '''
        Parses uniport speclist.txt    https://www.uniprot.org/docs/speclist.txt
        to local .CSV
        
        taxonid   species   lineanname       commonname
        72259      ABANI    Abaeis nicippe   Sleepy orange butterfly
                                             
        OXYMO E  475340: N=Oxytenis modestia
                         C=Costa Rica leaf moth
                         S=Dead-leaf moth
        
        
        '''
        listfile = self.speciesmap
        self.log.debug("Opening species map file %s" % listfile)
        try:
            fh = open(listfile, 'r')
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)                
       
        species = None
        kingdom = None
        taxonid = None
        lineanname = None
        commonname = None
        
        columnnames = ['species','kingdom','taxonid','lineanname','commonname']
        datalist = []
        # list of tuples
              
        try:
            for line in fh:
                #self.log.debug("handling line %s" % line)
                if 'N=' in line and not line.startswith('Code')  :
                    #self.log.debug("handling N= line. taxonid is %s" % taxonid)
                    if species is not None:
                        tup = (species, kingdom, taxonid, lineanname, commonname)
                        #self.log.debug("Adding tuple: %s" % str(tup))
                        datalist.append( tup )
                        # reset all varaiables
                        species = kingdom = taxonid = lineanname = commonname = None
                    species = line[:5]
                    kingdom = line[6]
                    taxonid = line[7:15].strip()
                    lineanname = line[19:].strip()
                    #self.log.debug("handling N= line. taxonid is %s" % taxonid)         
                elif 'C=' in line :
                    commonname = line[19:].strip()
                elif 'S=' in line :
                     pass
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        finally:
            fh.close()
        
        self.log.debug("Parsed file with %d terms" % len(datalist) )
    
        df = pd.DataFrame(datalist, columns = columnnames)
        outfile = "%s/speclist.csv" % self.outdir
        self.log.debug("Writing dataframe to %s" % outfile )
        df.to_csv(outfile )
        print(str(df))
        return df

        

def test_uniprot(config):
    upg = UniProtGOPlugin(config)
    entrylist = ['Q9CQV8', 'P35213', 'A4K2U9', 'P31946', 'Q4R572', 'P68250']
    out = upg._query_entries(entrylist)
    print(out)     

def test_datparse(config):
    upg = UniProtGOPlugin(config)
    df = upg.get_swissprot_df()
    return df

def test_speciesmap(config):
    upg = UniProtGOPlugin(config)
    upg._make_species_map()

def test_testset(config):
    upg = UniProtGOPlugin(config)
    df = upg.get_annotation_df()
    return df
    


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

    #test_uniprot(cp)
    df = test_testset(cp)
    df.to_csv('testset.csv')
    print(str(df))
    #test_speciesmap(cp)
