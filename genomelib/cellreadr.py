import logging
import operator
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq 
# https://github.com/chapmanb/bcbb/tree/master/gff

from pprint import pformat

def most_shared_feature(gene_sf, min_length=200, feature_type='exon',top_n=3):
    '''
    Take SeqFeature (gene)
            SeqFeature.sub_features -> List of SeqFeatures  (transcripts/splice-variants)
                ->SeqFeature.sub_features -> List of SeqFeatures (exon|UTR|stop_codon|CDS|transcript)
    
    Extract exons. 
    Find longest exon common to all transcripts.             
    Order exons by commonality. ?
    
    @arg SeqFeature (gene w/ transcript SeqFeature sub_features.   
    
    @return ordered list of "best" exon
      
    Create dict of exon start position   exon_dict -> { '100934649' : [SeqFeature(exon), SeqFeature(exon) ],
                                                        '100900312' : [SeqFeature(exon)],
                                                    }
    
    Preference:  most-shared feature between variants that meets minimum length. 
    Note:  Assumes all versions of a feature at a given location start are the same between transcripts.   
    '''
    logging.debug(gene_sf)
    exon_dict = defaultdict(list)
    len_list = []
    
    #
    #    
    #
    for tsf in gene_sf.sub_features:
        if tsf.type=='inferred_parent':
            for sf in tsf.sub_features:
                if sf.type == feature_type:
                    sf_loc = sf.location.start.position
                    sf_len = len(sf)
                    if sf_len >= min_length:
                        exon_dict[sf_loc].append(sf)
                        len_list.append( (sf_loc, sf_len) ) 
    
    logging.debug(pformat(dict(exon_dict)))
    logging.debug(len_list)
    len_list = list(set(len_list))
    len_list.sort(key=operator.itemgetter(1), reverse=True)
    logging.debug(len_list)                       
    len_list = len_list[:top_n]
    out_sf_list = []
    for (sf_loc, sf_len) in len_list:
        print(f'looking up loc {sf_loc} with len {sf_len}')            
        out_sf_list.append(exon_dict[sf_loc][0])                
    logging.debug(f'found list of {len(out_sf_list)} exons > {min_length} bp. Top {top_n}. ')
    return out_sf_list   


def calc_cg_percent(feature_sf):
    pass


  
    

    
    
    
