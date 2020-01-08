#!/usr/bin/env python
#
#
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

import pandas as pd

from cafa4.cafalib import CAFAPlugin
from cafa4.ontology import GeneOntology



class GOPlugin(CAFAPlugin):
    """
    Pipeline object. 
    Takes Pandas dataframe 
        1) Walks up GO tree, fanning out each line with different GOterm, *increasing* cafaprob estimate.  

    Input:  Pandas DataFrame
    Output: Pandas DataFrame
        Added rows
        Added columns: cafaprob

    """
    REPR_ATTRS=['outdir']

    def __init__(self, config):
        """
        
        """
        super(GOPlugin, self).__init__(config)
        #self.configname = self.__class__.__name__.lower()
        self.go = GeneOntology(self.config)
        self.log.debug("GOPlugin initialized.")

    def execute(self, dataframe):
        # add column for probability
        
        dataframe['cafaprob'] = 0.01

        # fan out for each go term up to root...
        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            self.log.debug("inbound row = %s" % str(row))
            (cafaid, evalue, score, bias, db, proteinacc, protein, species, 
            cafaprot, cafaspec, goterm, goaspect, goevidence, cafaprob) = row[1:]
            self.log.debug(f"Searching for parents of '{goterm}'...")        
            gt = self.go.get_term(goterm)
            isa_list = gt.get_isalist()
            self.log.debug(f"Got list of parents of {goterm} : {isa_list}")
            # Add row for existing base goterm:
            newrow = [cafaid, evalue, score, bias, db, 
                         proteinacc, protein, species, cafaprot, 
                         cafaspec, goterm, goaspect, goevidence, cafaprob ]
            newdfdict[ix] = newrow
            ix += 1
            
            # Add row for each parent:
            for parentterm in isa_list:
                cafaprob = cafaprob + .01 
                newrow = [cafaid, evalue, score, bias, db, 
                          proteinacc, protein, species, cafaprot, 
                          cafaspec, parentterm, goaspect, goevidence, cafaprob ]
                newdfdict[ix] = newrow
                ix += 1
                
        self.log.debug(f"Processed {ix} rows for new DF.")
        newdf = pd.DataFrame.from_dict(newdfdict, 
                                       orient='index', 
                                       columns = ['cafaid', 
                                                  'evalue', 
                                                  'score', 
                                                  'bias', 
                                                  'db', 
                                                  'proteinacc', 
                                                  'protein', 
                                                  'species', 
                                                  'cafaprot', 
                                                  'cafaspec',
                                                  'goterm',
                                                  'goaspect',
                                                  'goevidence',
                                                  'cafaprob',
                                                  ] )    

        self.log.debug(f"New dataframe with {len(newdf.index)} rows.")
        return newdf
    

    
    