#!/usr/bin/env python
#
#  Consume .obo file and .gaf
#  output dataframe
#  save/load CSV
#          
#            F1   F2   f3   f4  ...
#   gene 1   0     1    0    1
#   gene 2   1     0    0    0
#   gene 3   0     1    1    0 
#    ...  
#

import argparse
import collections
from configparser import ConfigParser
import itertools
import logging
import os
import sys
import traceback

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


class GOMatrix(object):
    """
    matrix of genes/proteins x goterms
    optionally extended by all is_a relationships..



    
    """
    def __init__(self, config, gaffile):
        #self.kname = self.__class__.__name__
        self.log = logging.getLogger(self.__class__.__name__)
        self.config = config
        self.cachedir = os.path.expanduser(config.get('ontology','cachedir'))
        self.cachefile = 'gomatrix.csv'
        self.cachepath = "%s/%s" % (self.cachedir, self.cachefile)
        self.df = None
        self.log.debug("GOMatrix initialized.")


    def get_df(self, cache=True):
        self.df = self._df_from_cache()
        if self.df is None:
            self.execute()
        return self.df
    
    
    def _df_from_cache(self):
        self.log.debug("Trying to load from cache: %s" % self.cachepath )  
        self.df = pd.read_csv(self.cachepath)
        return self.df


class GOTerm(object):
    #  Cache of showing all (non-redundant) is_a parents for given GO Term. 
    ISA_LISTCACHE = {}
    
    def __init__(self ):
        self.log = logging.getLogger(self.__class__.__name__)
        self.goterm = None
        self.name = None
        self.goaspect = None
        self.definition = None
        self.synonym = []
        self.is_a = []
    
    def get_isalist(self):
        '''
        Return non-redundant list of all goterm strings that are is_a parents of this term
        '''
        # Root terms
        nl = None
        
        if self.goterm in GOTerm.ISA_LISTCACHE.keys() :
                self.log.debug("%s Cache hit, using..." % self.goterm)
                nl = GOTerm.ISA_LISTCACHE[self.goterm]
        
        elif len(self.is_a) == 0:
            self.log.debug("%s Cache miss" % self.goterm)
            self.log.debug("%s ROOT TERM!" % self.goterm)
            nl = []
            self.log.debug("%s Storing cache entry= %s" % (self.goterm, nl))
            GOTerm.ISA_LISTCACHE[self.goterm] = nl 

        # Only has one parent, simply construct list
        elif len(self.is_a) == 1:
            # no entry for this isa_list. Must construct...
            self.log.debug("%s Cache miss" % self.goterm)
            item = self.is_a[0]
            itemlist = [item.goterm] + item.get_isalist()
            self.log.debug("%s Got itemlist %s for %s" % (self.goterm, itemlist, item.goterm))
            nl =  itemlist
            self.log.debug("%s new list is %s" % (self.goterm, nl))
            # add item term to item's isalist. 
            #nl = list(collections.OrderedDict.fromkeys(gl))
            #nl = list(self.goterm).extend(itemlist)
            self.log.debug("%s Storing cache entry= %s" % (self.goterm, nl))
            GOTerm.ISA_LISTCACHE[self.goterm] = nl
        
        # Has multiple paths to root, must use sets to de-duplicate. 
        else:
            self.log.debug("%s Cache miss" % self.goterm)   
            gl = []              
            for item in self.is_a:
                gl.extend( [item.goterm] + item.get_isalist())
            
            nl = list(collections.OrderedDict.fromkeys(gl))
            # store de-duped list to cache, preserving order
            self.log.debug("%s Storing cache entry= %s" % (self.goterm, nl))
            GOTerm.ISA_LISTCACHE[self.goterm] = nl

        self.log.debug("%s returning %s" % (self.goterm, nl))
        return nl    
    
    def get_isastr(self):
        '''
        Return string representation of all is_a parents of this term. 
        
        '''
        il = self.get_isalist()
        s = ' '.join(il)
        return s   
    
    def __repr__(self):
        s = "GOTerm:"
        for atr in ['goterm','name', 'goaspect']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s
        

class GeneOntology(object):
              
    NSMAP= { 'biological_process' : 'bp',
             'cellular_component' : 'cc',
             'molecular_function' : 'mf',
             'external'           : 'ex'
            }
        
    def __init__(self, config, obofile=None):
        self.config = config
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.cachedir = os.path.expanduser(config.get(self.lkname ,'cachedir'))
        self.cachefile = "%s.csv" % self.lkname
        self.cachepath = "%s/%s" % (self.cachedir, self.cachefile)
                
        if obofile is None:
            self.obofile = os.path.expanduser(self.config.get('ontology','obofile'))
        self.log = logging.getLogger(self.__class__.__name__)
        self.goidx = {}   #  { 'GO:XXXXXX' : GoTermObject, }      
        self.isalists = {}
        self.df = None
        self.log.debug("Object initialized.")
        
        
    def get_df(self):
        '''
        Create dataframe from local file for further usage...
        Used cached version on disk if available. 
        '''
        self._df_from_cache()
        if self.df is not None:
            self.log.debug("Cache hit. Using DataFrame from cache...")
        else:
            self.log.debug("Cache miss. Regenerating DF...")
            data = self.get_dict()
            #self.log.debug("godict = %s" % data)
            df = pd.DataFrame.from_dict(data, orient='index', columns=['goterm','name','goaspect']) 
            df.set_index('goterm')
            df.to_csv(self.cachepath)
            self.df = df
        self.log.debug(str(self.df))
        return self.df
    
    def get_dict(self):
        dict = self._handle_obo(self.obofile)
        return dict

    def get_tree_termindex(self):
        '''
        Builds full object tree, 
        provides dictionary with GOTerms as keys, object in tree as value. 
        object is_a/part of relations can be used to traverse tree.
        
        '''
        filehandle = None
        try:
            self.log.debug("opening file %s" % self.obofile)
            filehandle = open(self.obofile, 'r')
            self._parse2tree(filehandle)
            filehandle.close()
                
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)        
        
        finally:
            if filehandle is not None:
                filehandle.close()
        self._add_references()
        
        return self.goidx

    def _df_from_cache(self):
        if os.path.exists(self.cachepath):
            self.log.debug("Trying read from cachepath %s" % self.cachepath)
            self.df = pd.read_csv(self.cachepath, index_col=0)
            self.log.debug("Loaded DF from %s" % self.cachepath)
        else:
            self.log.debug("No file at cachepath: %s" % self.cachepath)
        
        
    def _handle_obo(self, filename):
        '''
        Create dictionary from obo:
        
        { 0  : [ 'GO:0000002', 'mitochondrial genome maintenance','biological_process'],
          1  : [ 'GO:0000186', 'activation of MAPKK activity','biological_process']
        }  
        
        '''       
        godict = {}
        try:
            self.log.debug("opening file %s" % filename)
            filehandle = open(filename, 'r')
            godict = self._parsefile(filehandle)
            filehandle.close()
                
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)                
        return godict    

    def _handle_file(self, filename):
        try:
            self.log.debug("opening file %s" % filename)
            filehandle = open(filename, 'r')
            self._parsefile(filehandle)
            filehandle.close()
                
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)                  

        self._add_references()

    def _parse2tree(self, filehandle):
        '''
        Create in-memory GO tree with GOTerm index. 
        Two passes:
            1. Do not resolve foreign references, just add them as strings...
            2. Resolve foreign refs, replacing strings w/ object refs. 
                
        goidx = { 'GO:0000001' : <GOTerm Object>,
                  'GO:0000002' : <<GOTerm Object>,
        }
        '''      
        self.goidx = {}
        current = None
        self.log.info("Parsing file %s" % self.obofile)
        try:
            for line in filehandle:
                #self.log.debug("line is %s" % line)
                if line.startswith("[Term]"):     
                    #self.log.debug("found term")
                    if current is not None:
                        self.goidx[current.goterm]= current
                        current = GOTerm()
                    else:
                        #self.log.debug("Creating new GOTerm object...")
                        current = GOTerm()
                    
                elif line.startswith("id: "):
                    current.goterm = line[4:].strip()
                    
                elif line.startswith("name: "):
                    current.name = line[6:].strip()
                
                elif line.startswith("namespace: "):
                    current.namespace = line[11:].strip()
                
                elif line.startswith("def: "):
                    current.definition = line[5:].strip()

                elif line.startswith("synonym: "):
                    current.synonym.append(line[9:].strip())

                elif line.startswith("is_a: "):
                    current.is_a.append(line[6:16].strip())
                  
        
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        self.log.info("Parsed file with %d terms" % len(self.goidx) )

    def _add_references(self):
        '''
        Go through goidx and add foreign refs...
        '''
        for gt in self.goidx.keys():
            gto = self.goidx[gt]
            isalist = gto.is_a
            newisa = []
            for gt in isalist:
                newisa.append(self.goidx[gt])
            gto.is_a = newisa        
            
            
    def _parsefile(self, filehandle):
        od = {}
        keyid = 0
        current = None
        try:
            for line in filehandle:
                if line.startswith("[Term]"):
                    #self.log.debug("found term")
                    if current is not None:
                        od[keyid] = current
                        current = []
                        keyid += 1
                    else:
                        current = []
                    
                elif line.startswith("id: "):
                    #self.log.debug("found id:")
                    (key, val ) = line.split(":",1) # only parse to first occurence of ':'
                    val = val.strip()
                    current.append(val)
                
                elif line.startswith("name: "):
                    #self.log.debug("found name")
                    (key, val ) = line.split(":",1)
                    val = val.strip()
                    current.append(val)
                
                elif line.startswith("namespace: "):
                    #self.log.debug("found namespace")
                    (key, val ) = line.split(":",1)
                    val = val.strip()
                    val = GeneOntology.NSMAP[val]
                    current.append(val)

                elif line.strip().startswith("#"):
                    pass          
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        self.log.debug("Parsed file with %d terms" % keyid )
        return od  

def test(config):
    go = GeneOntology(config)
    df = go.get_df()
    print(str(df))
    df = go.get_df()
    
    goidx = go.get_tree_termindex()
    obj = goidx['GO:0009987']
    print(obj)
    
    testterms = os.path.expanduser('~/play/cafa4/gocodes.txt')
    #testterms = os.path.expanduser('~/play/cafa4/smallgocodes.txt')

    golists = {}
    fh = open(testterms, 'r')
    for line in fh:
        gt = line.strip()
        print("looking up term: %s" % gt)
        gtobj = goidx[gt]
        print("goterm: %s -> is_a: %s" % (gt,  gtobj.get_isastr()) )
    
    #print("ISA_LISTCACHE=%s" % GOTerm.ISA_LISTCACHE)    
    #df = go.get_df()
    #print(str(df))



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

    parser.add_argument('-t', '--test', 
                        action="store_true", 
                        dest='test', 
                        help='run tests')

    parser.add_argument('-g', '--gaffile', 
                        action="store", 
                        dest='gaffile', 
                        default=None,
                        help='Gene annotation file')

    parser.add_argument('goterms', 
                        metavar='goterms', 
                        type=str, 
                        help='one or more space-separated goterms GO:XXXXXXXX' ,
                        nargs='*')
                   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)    

    config = ConfigParser() 
    config.read(args.conffile)
       
    if args.test:
        test(config)
    
    if len(args.goterms) > 0:
        go = GeneOntology(config)
        goidx = go.get_tree_termindex()
        for gt in args.goterms:
            gobj = goidx[gt]
            print("%s -> %s " % (gt, gobj.get_isalist()))
    
    if args.gaffile is not None:
        gm = GOMatrix(config, args.gaffile)
    