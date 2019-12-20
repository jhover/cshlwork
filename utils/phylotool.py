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
import logging

from Bio import Phylo
import numpy as np


class Phylogeny(object):
    
    def __init__(self):
        self.kname = self.__class__.__name__
        self.log = logging.getLogger(self.kname)
        
    def parsefile(self, filepath):
        """
         Reads NHX format file,...
        """
        self.log.debug("Reading file %s" % filepath)
        tree = Phylo.read(filepath, 'newick')
        print(tree)
        
    def to_distance_matrix(tree):
        """Create a distance matrix (NumPy array) from clades/branches in tree.
    
        A cell (i,j) in the array is the length of the branch between allclades[i]
        and allclades[j], if a branch exists, otherwise infinity.
    
        Returns a tuple of (allclades, distance_matrix) where allclades is a list of
        clades and distance_matrix is a NumPy 2D array.
        """
        allclades = list(tree.find_clades(order='level'))
        lookup = {}
        for i, elem in enumerate(allclades):
            lookup[elem] = i
        distmat = numpy.repeat(numpy.inf, len(allclades)**2)
        distmat.shape = (len(allclades), len(allclades))
        for parent in tree.find_clades(terminal=False, order='level'):
            for child in parent.clades:
                if child.branch_length:
                    distmat[lookup[parent], lookup[child]] = child.branch_length
        if not tree.rooted:
            distmat += distmat.transpose
        return (allclades, numpy.matrix(distmat))



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
         
         