#!/usr/bin/env python
#

import argparse
import itertools
import logging
import os
import random

from collections import OrderedDict

import numpy as np
import pandas as pdd

import dendropy as dp


class Phylogeny(object):

    def __init__(self, filepath=None):
        self.log = logging.getLogger(self.__class__.__name__)
        self.tree = None
        self.distmatrix = None
        self.df = None
        self.filepath = filepath
        
    def __repr__(self):
        pass
        
    def parsefile(self, filepath=None):
        """
         Reads NHX format file,...
        """
        self.log.debug("Reading file %s" % filepath)
        self.filepath = filepath
        self.tree = dp.Tree.get_from_path(args.infile, 
                                 schema='newick',
                                 #suppress_internal_node_taxa=True, 
                                 #suppress_leaf_node_taxa=True,
                                 )
        
        self.tree = Phylo.read(filepath, 'newick')
        self.log.debug("tree is %s" % self.tree )
        #print(tree)
  


def process_node(node, start=1.0):
    if node.parent_node is None:
        node.value = start
    else:
        node.value = random.gauss(node.parent_node.value, node.edge.length)
    for child in node.child_nodes():
        process_node(child)
    if node.taxon is not None:
        print("%s : %s" % (node.taxon, node.value))


def test_dendropy(infile):
        
    p = Phylogeny()
    p.parsefile(infile)
        
    
    #print(f"description = {tree.description()}")
    print(f"seed_node = {tree.seed_node}")
    print(f"seed_node.taxon = {tree.seed_node.taxon}")
    nt = tree.find_node_with_taxon_label("Gluconobacter japonicus")
    print(f"nt = {nt}")
    print(nt.adjacent_nodes())
    print(nt.edge)
    print(nt.edge.length)
    print(nt.parent_node)
    print(nt.edge_length)
    print(nt.distance_from_root() )
    #process_node(nt)
    
    #process_node(tree.seed_node)





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
                        help='a phylogeny file NHX/PHYLIP')
    
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
        
    logging.info("dendrotool...")
    
    test_dendropy(args.infile)

    
    
        