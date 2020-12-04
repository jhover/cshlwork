#!/usr/bin/env python
import os
import sys
import logging
import traceback

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

from fastcafa.fastcafa import *

phmmerdf = os.path.expanduser("~/play/hamsini/dupe_targets_phmmer_20200926.csv")
dupepairs = os.path.expanduser("~/play/hamsini/mouse_dup_pairs_uniprot.txt")
evaltable = os.path.expanduser("~/play/hamsini/mouse_dupe_scores.csv")


def parse_dupepairs():
    f = open(dupepairs)
    lines = f.readlines()
    dupelist = []
    for line in lines:
        (p1, p2) = line.split(',')
        p1 = p1.strip()
        p2 = p2.strip()
        dupelist.append( (p1, p2) )
    logging.debug(f"dupelist[1] = {dupelist[1]}")
    return dupelist

def get_match(query, target, df):
    logging.debug(f"query={query} target={target}")
    qdf = df[df.cid == query]
    row = qdf[qdf.pacc == target]
    if len(row) > 1 :
        logging.warning(f'multiple matches for query={query} target={target} ')
        return None
    elif len(row) == 1:
        r = row.iloc[0]
        eval = r['eval']
        pscore =r['pscore']
        bias = r['bias']
        return (eval, pscore, bias)
    else:
        logging.warning(f'no matches for query={query} target={target} ')    
        return None


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    #logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().setLevel(logging.DEBUG)
    logging.debug("make evaltable...")
    config = get_default_config()

    pdf  = pd.read_csv(phmmerdf, index_col=0)
    pdf.drop_duplicates(inplace=True,ignore_index=True)
    
    dupelist = parse_dupepairs()
    
    lod = []
        
    for tup in dupelist:
        (p1, p2) = tup
        logging.debug(f"looking for {p1} -> {p2}")
        rv = get_match(p1, p2, pdf)
        if rv is not None:
            (eval, pscore, bias ) = rv
            lod.append( { 'query' : p1,
                          'target' : p2,
                          'eval' : eval,
                          'pscore' : pscore,
                          'bias' : bias,  
                            }
            )
    logging.debug(f"dupelist length: {len(dupelist)}")
    logging.debug(f"matchlist length: {len(lod)}")
    edf = pd.DataFrame(lod)
    edf.to_csv(evaltable)
    logging.debug(f"wrote match df to {evaltable}")

    # E.g. Q9QY40,Q3UH93
    



    
    
