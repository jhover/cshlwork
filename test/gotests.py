#!/usr/bin/env python
#
#
#
import argparse 
import logging

import numpy as np
from scipy import sparse


def testgo():
    
    # propagated go matrix. representing simple tree. 
    #
    #       1
    #      / \
    #     2   3
    #        /  \
    #        4  5

    gm = [(1,0,0,0,0), 
          (1,1,0,0,0), 
          (1,0,1,0,0), 
          (1,0,1,1,0), 
          (1,0,1,0,1), 
        ]
    gomatrix = np.asarray(gm, dtype=bool)

    #  gene golist
    #  
    #   G1:  GO:3
    #   G2:  GO:4, Go:5
    #   G3:  GO:1
    #  
    gg = [(0,0,1,0,0), 
          (0,0,0,1,1), 
          (0,1,0,0,0), 
         ]
    
    genebygo = np.asarray(gg, dtype=bool)         
    
    #print(genebygo)
    print(genebygo @ gomatrix)




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
          
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    testgo()