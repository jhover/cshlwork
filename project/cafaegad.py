#!/usr/bin/env python
#
# consume CAFA prediction file. 
# 
#
#
#
#
  
import os
import sys
import logging
import traceback

gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

from expression.egad import *


def parse_predout(filepath):
    
#
# G803000000001    GO:0005667    0.10
# G803000000001    GO:0043966    0.10
# G803000000001    GO:0045893    0.10
# G803000000001    GO:0001932    0.11
# G803000000001    GO:0048589    0.11
# G803000000001    GO:0000993    0.12
#
# return list of lists 
#
    try:
        logging.debug(f" attempting to open '{filepath}'")
        f = open(filepath, 'r')
        allentries = []
       
        for line in f.readlines(): 
            fields = line.split()
            geneid = fields[0]
            goterm = fields[1]
            prob = float(fields[2])
            allentries.append( [geneid, goterm, prob]   )

    except FileNotFoundError:
        logging.error(f"No such file {filepath}")   

    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    if f is not None:
        f.close()      
    
    logging.info(f"Parsed file with {len(allentries)} entries" )
    logging.debug(f"Some entries:  {allentries[1000:1005]}")
    return allentries









    
if __name__ == "__main__":

    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    #logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().setLevel(logging.DEBUG)

    logging.debug("CAFA to EGAD..")
    
    predfile = os.path.expanduser("~/play/jones/gillis_seqs.predout")
    
    lol = parse_predout(predfile)
    
    
    
    
