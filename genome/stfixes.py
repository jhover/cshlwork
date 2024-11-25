import logging
import os
import sys
import traceback


gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *


from BCBio import GFF
from Bio import SeqIO


def fix_annot(annot_list):
    '''
    restore hierarchy of UTexas GTF 
    features=['gene','inferred_parent'],  # top level object 
              sub_features=['mRNA'],  # mid-level to handle.        
              seq_feature='exon', #  lowest-level sequences to be extracted.
    '''
    # MAP = { 'transcript': ['exon', 'CDS']   }
  
    logging.debug(f'handling annot_list of {len(annot_list)} SeqRecord(s)...')   
    for sr in annot_list:
        logging.debug(f'handling SeqRecord(id = {sr.id})')
        # no type
        # e.g. SeqRecord(id ='scf7180000321923.62891').features
        for f in sr.features:
            logging.debug(f'handling SeqFeature(type= {f.type} id = {f.id})')
            f.type = 'gene'
            # types:  biological_region, enhancer, pseudogene, gene,   
            # e.g. SeqFeature(type='gene', id='evm.TU.scf7180000321923.62891.1').sub_features
            mrna = None
            sflist = []
            for sf in f.sub_features:
                # e.g.
                #   [SeqFeature(SimpleLocation(ExactPosition(11726), ExactPosition(14477), strand=-1), type='transcript', qualifiers=...),
                #    SeqFeature(SimpleLocation(ExactPosition(11726), ExactPosition(12644), strand=-1), type='exon', qualifiers=...),
                #    SeqFeature(SimpleLocation(ExactPosition(12787), ExactPosition(12842), strand=-1), type='exon', qualifiers=...),
                logging.debug(f'handling sub-SeqFeature(type= {sf.type} id = {sf.id})')
                if sf.type == 'transcript':
                    logging.debug(f'found transcript for {f.id}')
                    mrna = sf
                elif sf.type == 'exon':
                    sflist.append(sf)
                else:
                    logging.warning(f'sf.type = {sf.type}')
            logging.debug(f'sflist is {sflist}')        
            if ( mrna is not None) and ( len(sflist) > 0):
                logging.debug('found sane transcript.')
                for sf in sflist:
                    sf.qualifiers['Parent'] = mrna.qualifiers['transcript_name']
                mrna.sub_features = sflist
                #mrna.qualifiers['Parent'] = f.qualifiers['ID']
                mrna.qualifiers['ID'] = mrna.qualifiers['transcript_name']
                f.sub_features = [ mrna ]        
    return annot_list



def save_annotation(annot_list, outbase='new_annotation'):
    with open(f'{outbase}.gff3', 'w') as fh:
        GFF.write(annot_list, fh )
        logging.debug(f'done writing {outbase}.gff3')
    
     

