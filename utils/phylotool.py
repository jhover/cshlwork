#!/usr/bin/env python
#
#
#  https://www.biostars.org/p/224572/
#  preferred breakdown:    animals,   fungi,  plants,   protists,   bacteria,  archaea
#  ncbi taxon category codes:  chordata     7711   vertebrata   7742  Insecta  50557 
#                                # in DB
#   1 => 'Eukaryota',
#        2 => 'Animals',         1653
#          Vertebrata            877
#            7 => 'Mammals', 
#            3 => 'Birds',
#            9 => 'Reptiles',
#            2 => 'Amphibians', 
#            4 => 'Fishes', 
#          Invertebrates  
#            5 => 'Flatworms',    38
#            10 => 'Roundworms'   119
#            6 => 'Insects',      412
#            8 => 'Other Animals', },
#        3 => 'Fungi', 
#        4 => 'Other', 
#        5 => 'Plants',         459
#        6 => 'Protists'},DB
#   2 => 'Bacteria' -> all, 
#   3 => 'Archaea',
#   4 => 'Viroids', 
#   5 => 'Viruses');
#   
#   areas?:   eukaryote->animals
#                          amph, birds, fishes, mammals, reptiles   |  flatworms , roundworms 
#            eukaryote->plants               
#            eukaryotes->fungi
#            eukaryotes->protists
#   
#  Convert newick tree to distance matrix in R:  https://www.biostars.org/p/312148/
#    treeText <- readLines(tree.phy)
#        treeText <- paste0(treeText, collapse="")
#        library(treeio)
#        tree <- read.tree(text = treeText) ## load tree 
#    distMat <- cophenetic(tree) ## generate dist matrix
#
#     See: https://biopython.org/wiki/Phylo_cookbook
#     for distance matrix calculation using biopython
#
#
#

import argparse
import itertools
import logging
import os

from collections import OrderedDict

from Bio import Phylo
import numpy as np
import pandas as pd

#pd.set_option("display.max_rows", 15)
#pd.set_option("display.max_columns", 15)


class Phylogeny(object):
    
    def __init__(self):
        self.log = logging.getLogger(self.__class__.__name__)
        self.tree = None
        self.distmatrix = None
        self.df = None
        self.filepath = None
        
    def __repr__(self):
        pass
        
    def parsefile(self, filepath):
        """
         Reads NHX format file,...
        """
        self.log.debug("Reading file %s" % filepath)
        self.filepath = filepath
        self.tree = Phylo.read(filepath, 'newick')
        self.log.debug("tree is %s" % self.tree )
        #print(tree)
       
        
    def get_distance_matrix(self):
        """
        
        We want a list of row/column names, and a corresponding 2D matrix
        
        names = ['75743','137246','7950']
        
        data = 
        
        Note: Takes about ~2 seconds with 119 terminals on a 2019 Macbook Pro
              Takes about ~38 seconds with 459 terminals on a 2019 Macbook Pro    
              Takes about ~7 minutes with 876 terminals on a 2019 Macbook Pro
        
        """
        self.distmatrix = OrderedDict()
        terminals = self.tree.get_terminals()
        terminals.sort(key=lambda x: x.name, reverse=True)
        mdim = len(terminals)
        self.log.debug("%d  terminals.." % mdim)
        i = 0
        for x,y in itertools.combinations_with_replacement(terminals, 2):
            if i % 1000 == 0:
                self.log.debug("Handling combination %d" % i)
            v = self.tree.distance(x, y)
            self.distmatrix[x.name] = self.distmatrix.get(x.name, OrderedDict())
            self.distmatrix[x.name][y.name] = v
            self.distmatrix[y.name] = self.distmatrix.get(y.name, OrderedDict())
            self.distmatrix[y.name][x.name] = v
            i += 1
        
        self.log.debug("Done computing distances. Filling diagonal...")
        for x in terminals:
            self.distmatrix[x.name][x.name] = 0
        
        self.log.debug(self.distmatrix)

        colnames = list(self.distmatrix.keys())
        return ( colnames, self.distmatrix )        


    def to_df(self):
        #(allclades, distmatrix) = self.to_distance_matrix()
        #df = pd.DataFrame(data = distmatrix,
        #                  index = allclades,
        #                  columns = allclades)
        csvpath = "%s.csv" % self.filepath
        if os.path.exists(csvpath):
            self.df = pd.read_csv(csvpath)
        
        else:
            if self.distmatrix is not None:
                self.log.debug("Found completed distmatrix. Converting...")
            else:
                self.log.debug("No distmatrix found. Computing...")
                self.get_distance_matrix()
            self.df = pd.DataFrame(self.distmatrix, columns = self.distmatrix.keys())    
        
        return self.df

    def to_csv(self):
        if self.df is not None:
            self.log.debug("DataFrame found. Using...")
        else:
            self.log.debug("No DataFrame found. Computing...")
            self.to_df()
        self.df.to_csv("%s.csv" % self.filepath)
    

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

    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='a phylegeny file NHX')
    
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default='~/etc/phylo.conf',
                        help='Config file path [~/etc/phylo.conf]')

                    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    #cp = ConfigParser()
    #cp.read(args.conffile)
         
    p = Phylogeny()
    p.parsefile(args.infile)
    #(terminals, matrix) = p.get_distance_matrix()
    #print(terminals)
    df = p.to_df()
    print(df)
    p.to_csv()
         