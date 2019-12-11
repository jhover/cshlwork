#!/usr/bin/env python
#
#  Look up GO annoatations from UniProtKB:entrys. 
#  E.g.  P35213  (1422B_RAT  gene: Ywhab
# 

import configparser
import logging
import requests
import sys

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

class QuickGOPlugin(object):
    
    def __init__(self, config):
        self.log = logging.getLogger()
        self.requestbase = "https://www.ebi.ac.uk/QuickGO/services/annotation"
        self.config = config

    def get_df(self,dataframe):
        '''
        Takes 
        Returns Pandas DataFrame
        
        '''
        entries = dataframe['tacc'].unique().tolist()
        txt = self._query_entries(entries)

    def _query_entries(self, entrylist):
        self.log.debug("querying entry list: %s" % entrylist)
        entrystr=','.join(entrylist)
        self.log.debug("querying entry string: %s" % entrystr)
        requestURL = "%s/downloadSearch?geneProductId=%s" % (self.requestbase, entrystr )
        self.log.debug("RequestURL=%s"% requestURL )        
        r = requests.get(requestURL, headers={ "Accept" : "text/tsv"})
        if not r.ok:
            r.raise_for_status()
            #sys.exit()
        response = r.text
        self.log.debug("response=%s" % response)
        return response
        

if __name__ == '__main__':
    config = configparser.ConfigParser()
    
    qg = QuickGOPlugin(config)
    entrylist = ['Q9CQV8', 'P35213', 'A4K2U9', 'P31946', 'Q4R572', 'P68250']
    out = qg._query_entries(entrylist)
    print(out)    

    #  curl -X GET --header 'Accept:text/tsv' 
    # 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=P0A9Z8'
    #
    # GENE PRODUCT DB    GENE PRODUCT ID    SYMBOL    QUALIFIER    GO TERM    GO ASPECT    ECO ID    GO EVIDENCE CODE    REFERENCE    WITH/FROM    TAXON ID    ASSIGNED BY    ANNOTATION EXTENSION    DATE
    #UniProtKB    P0A9Z8    bla    enables    GO:0008800    F    ECO:0000256    IEA    GO_REF:0000002    InterPro:IPR000871    573    InterPro    20191123
    #UniProtKB    P0A9Z8    bla    involved_in    GO:0030655    P    ECO:0000256    IEA    GO_REF:0000002    InterPro:IPR000871    573    InterPro        20191123
    #UniProtKB    P0A9Z8    bla    involved_in    GO:0046677    P    ECO:0000256    IEA    GO_REF:0000002    InterPro:IPR000871    573    InterPro        20191123
    #UniProtKB    P0A9Z8    bla    enables    GO:0008800    F    ECO:0000501    IEA    GO_REF:0000003    EC:3.5.2.6    573    UniProt        20191123
    #UniProtKB    P0A9Z8    bla    involved_in    GO:0046677    P    ECO:0000322    IEA    GO_REF:0000037    UniProtKB-KW:KW-0046    573    UniProt        20191123
    #UniProtKB    P0A9Z8    bla    enables    GO:0016787    F    ECO:0000322    IEA    GO_REF:0000037    UniProtKB-KW:KW-0378    573    UniProt        20191123
    #
    #






