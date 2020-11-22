#! /usr/bin/env python
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

from fastcafa.fastcafa import *

gitpath=os.path.expanduser("~/git/pyEGAD")
sys.path.append(gitpath)

from egad.egad import *

#SPECIES_NET=os.path.expanduser('~/data/cococonet/atlanticspecies_prioAggNet.hdf5')
SPECIES_NET=os.path.expanduser('~/data/cococonet/networks/afrog_prioAggNet.hdf5')


#PREDOUT=os.path.expanduser('~/play/jones/gillis_seqs.predout')
PREDOUT=os.path.expanduser('~/play/jones/afrog_prio.predout')
#  G803000000001    GO:0005667    0.10
#  G803000000001    GO:0043966    0.10
#  G803000000001    GO:0045893    0.10

SEQ_IDMAP=os.path.expanduser('~/play/jones/species_hiprio_seqmap.tsv')
#  G803000000001 A0A1S3SK04_SALSA
#  G803000000002 A0A1S3RA14_SALSA
#  G803000000003 A0A1S3RDQ3_SALSA

UID_GN_MAP=os.path.expanduser('~/play/jones/uniprot-trembl-species.8030.map.tsv')
#  db  acc          uid                 gn
#  tr  A0A1S3RID5   A0A1S3RID5_SALSA    LOC106602976
#  tr  B5XFF4       B5XFF4_SALSA        WRB
#

OUTFILE=f"{PREDOUT}.results.tsv"


def read_network_hdf5(filename):
    """
    Loads data in file to dataframe.
    """
    with h5py.File(filename, 'r') as f:
        logging.debug("reading matrix...")
        matrix = f['agg'][:]
        logging.debug("reading rows. converting to unicode.")
        rows = [ s.decode() for s in  f['row'][:] ]
        logging.debug("reading columns. converting to unicode")
        columns = [ s.decode() for s in  f['col'][:] ]
        logging.debug("making dataframe...")
        df = pd.DataFrame(matrix,  index=rows, columns = columns )
        logging.debug(f"network shape: {df.shape}")
        logging.debug(f"network:\n {df}")        
    return df        
    

def read_predout(predout):
    columns = ['seqid','goterm','prob'] 
    df = pd.read_csv(predout, sep='\t', header=None, names=columns)
    logging.debug(f"predout\n{df.head()}")
    logging.debug(f"predout shape: {df.shape}")
    return df


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)
    
    logging.info(f"Reading network: {SPECIES_NET}")
    nw = read_network_hdf5(SPECIES_NET)
    logging.info(f"SPECIES_NETwork:\n{nw}")
    #nw = fix_rowcol_names(nw, UID_GN_MAP)
    #logging.info(f"fixed SPECIES_NETwork:\n{nw}")        
    
    logging.info(f"Reading predictions: {PREDOUT}")
    po = read_predout(PREDOUT)   
    #po.to_csv(f"{PREDOUT}.csv", sep="\t")
    logging.info(f"\n{po}")    
    
    amdf = build_annotation_matrix(po, 'seqid','goterm')
    logging.info(f"\n{amdf}")    
    
    logging.info(f"input to run_egad: genesXgo:\n{po}\ngenesXgenes:\n{nw}")    
    outdf = run_egad(amdf, nw )
    logging.info(f"\n{outdf}")
    outdf.to_csv(f"{OUTFILE}", sep='\t')    
    logging.info(f"Wrote to {OUTFILE}")