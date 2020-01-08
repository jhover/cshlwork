#!/usr/bin/env python
#
#
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

from cafa4.cafalib import CAFAPlugin
from cafa4.ontology import GeneOntology

class GOPlugin(CAFAPlugin):
    """
    Pipeline object. 
    Takes Pandas dataframe 
        1) Walks up GO tree, fanning out each line with different GOterm, *increasing* likelihood/prob estimate.  

    Input:  Pandas DataFrame
    Output: Pandas DataFrame
        Added rows
        Added columns: goaspect  goevidence cafapest

    """
    REPR_ATTRS=['outdir']

    def __init__(self, config):
        """
        
        """
        super(GOPlugin, self).__init__(config)
        #self.configname = self.__class__.__name__.lower()
        self.go = GeneOntology(self.config)
        


    def execute(self, dataframe):
        gdf = self.go.get_df()
        self.log.debug("\n%a" % str(gdf))        
        dataframe['probest'] = 1.00
        
        #igdf = gdf.set_index('goterm')
        #gdict = igdf.to_dict('index')
        #dataframe['goaspect'] = dataframe.apply(
        #    lambda row: gdict[row.goterm]['goaspect'], 
        #    axis=1)
        self.log.debug(str(dataframe))
        return dataframe
    
    
    
    
    
    
    
    