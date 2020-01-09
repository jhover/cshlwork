#!/usr/bin/env python
#
#
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

import pandas as pd
import pronto as pt

from cafa4.cafalib import CAFAPlugin
from cafa4.ontology import GeneOntology



class OntologyPlugin(CAFAPlugin):
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
        super(OntologyPlugin, self).__init__(config)
        self.log.debug("Creating ontology with pronto...")
        self.pgo = pt.Ontology(os.path.expanduser(self.config.get(self.lkname,'obofile')))
        self.initprob = float(self.config.get(self.lkname,'initial_probability'))
        self.probstep = float(self.config.get(self.lkname,'probability_step'))
        self.log.debug("OntologyPlugin initialized.")


    def execute(self, dataframe):
        # add column for probability
        dataframe['cafaprob'] = self.initprob

        # fan out for each go term up to root...
        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            self.log.debug("inbound row = %s" % str(row))
            (cafaid, evalue, score, bias, db, proteinacc, protein, species, 
            cafaprot, cafaspec, goterm, goaspect, goevidence, cafaprob) = row[1:]
            self.log.debug(f"Searching for parents of '{goterm}'...")        
            gt = self.pgo.get_term(goterm)
            #try:
            #    getattr(gt, 'cafaprob')
            #except AttributeError:
            #    gt.cafaprob = self.initprob
            
            superclasses = gt.superclasses()
            self.log.debug(f"Got generator of superclasses for {goterm}")
            
            # Add row for existing base goterm:
            newrow = [cafaid, evalue, score, bias, db, 
                         proteinacc, protein, species, cafaprot, 
                         cafaspec, goterm, goaspect, goevidence, cafaprob ]
            newdfdict[ix] = newrow
            ix += 1
            
            # Add row for each parent:
            for sclass in superclasses:
                #try:
                #    getattr(gt, 'cafaprob')
                #except AttributeError:
                #    gt.cafaprob = self.initprob
                #cafaprob = sclass.cafaprob + self.probstep 
                cafaprob = cafaprob + self.probstep
                if cafaprob >= 1.0:
                    cafaprob = 0.99
                newrow = [cafaid, evalue, score, bias, db, 
                          proteinacc, protein, species, cafaprot, 
                          cafaspec, sclass.id , goaspect, goevidence, cafaprob ]
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

 
    
    