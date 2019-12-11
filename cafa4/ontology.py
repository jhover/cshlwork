#!/usr/bin/env python
#
#
#
import argparse
from configparser import ConfigParser
import logging
import os
import sys
import traceback

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

GOFILE='~/data/go/go.obo' 

class GeneOntologyGOInfoPlugin(object):
    
            
    NSMAP= { 'biological_process' : 'bp',
             'cellular_component' : 'cc',
             'molecular_function' : 'mf',
             'external'           : 'ex'
            }
    
    
    def __init__(self, config):
        self.gopath = os.path.expanduser(GOFILE)
        self.log = logging.getLogger()
        self.config = config
        
    def get_df(self):
        '''
        Create dataframe from local file for further usage...
        Used cached version on disk if available. 
        '''
        data = self.get_dict()
        df = pd.DataFrame.from_dict(data, orient='index', columns=['goterm','name','namespace']) 
        df.set_index('goterm')
        df.to_csv('go.csv')
        self.log.debug(str(df))
        return df
    
    
    def get_dict(self):
        dict = self._handle_obo(self.gopath)
        return dict
    
        
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
                   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)    
    
    go = GeneOntology()
    df = go.get_df()
    print(str(df))