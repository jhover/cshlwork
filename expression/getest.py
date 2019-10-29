#!/usr/bin/env python
#
# read in starout format file. 
# create Series, rank,
# merge Series to DataFrame
# 

#  STAR / Output format. 
# -quantMode GeneCounts 
# -hseq-count -s yes   /  -s reverse
# columns
# geneID       unstranded_counts   1st_stranded  2nd_stranded

# N_unmapped       2674617    2674617    2674617
# N_multimapping   2828822    2828822    2828822
# N_noFeature      5598327    24570638    24512243
# N_ambiguous      3871509    813790    812892
# ENSG00000223972.5    0    0    0
# ENSG00000227232.5    302    159    143
# ENSG00000278267.1    9    7    2
# ENSG00000243485.4    1    1    0
# ENSG00000237613.2    0    0    0
#

import argparse
import logging
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


class ExpressionDataset(object):
    
    def __init__(self, filelist = None, fileformat='starcounts'):
        self.log = logging.getLogger()
        self.dslist = []
        if filelist != None:
            self.handlefiles(filelist, fileformat)
            
    def handlefiles(self, filelist, fileformat='starcounts'):
        ## Use Series
        for f in filelist:
            ds = ExpressionDataset.parsefile2series(f)
            dsr = ds.rank()
            self.log.info("\n%s" % dsr)
            self.dslist.append(dsr)
      
        self.dataframe = pd.concat(self.dslist, axis=1)
        
        # (Re-)do pairwise correlation 
        self.corrdataframe = self.dataframe.corr(method='spearman')


    @classmethod
    def parsefile2series(cls, filename):
        logging.info("Processing file %s" % filename)
        (head, tail) = os.path.split(filename)    
        sname = tail.split('.')[0]
          
        f = open(filename)
        lines = f.readlines()
        data = {}          
        # handle values. 
        for line in lines[4:]:
            (label, unstrand, strand1, strand2) = [ f.strip() for f in line.split('\t') ]
            data[label] = strand2
    
        logging.info("data length is %s" % len(data))
        ds = pd.Series( data, name=sname )
        logging.info("\n%s" % ds.head() )  
        return ds    
   

    def plot(self):
        mask = np.zeros_like(self.corrdataframe)
        mask[np.triu_indices_from(mask)] =True
        
        ax = sns.heatmap( self.corrdataframe,
                          mask=mask,
                          vmin=.5, vmax=1, 
                          #cmap="YlGnBu",
                          center=0,
                          #cmap=sns.diverging_palette(20,220,n=200),
                          square=True )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            horizontalalignment='right')           
        plt.show()

   


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
                        help='a list of single-cell expression files')
    
    parser.add_argument('-p','--plot',
                        action="store_true",
                        default=False,
                        dest='plot', 
                        help='create plot')
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    filelist = args.infiles 
    logging.debug("%d files to process. " % len(filelist))
    
    eds = ExpressionDataset(filelist)
    
    if args.plot:
        eds.plot()
    

    
    
    