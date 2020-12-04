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

#SALMON_NET=os.path.expanduser('~/data/cococonet/atlanticsalmon_prioAggNet.hdf5')
#SALMON_NET=os.path.expanduser('~/data/cococonet/atlanticsalmon_metaAggNet.Rdata')
HUMAN_NET=os.path.expanduser('~/data/cococonet/human_prioAggNet.hdf5')

HUMAN_GOA=os.path.expanduser('~/data/goa/goa_human_gomatrix.csv')


PREDOUT=os.path.expanduser('~/play/jones/gillis_seqs.predout')
#  G803000000001    GO:0005667    0.10
#  G803000000001    GO:0043966    0.10
#  G803000000001    GO:0045893    0.10

SEQ_IDMAP=os.path.expanduser('~/play/jones/salmon_hiprio_seqmap.tsv')
#  G803000000001 A0A1S3SK04_SALSA
#  G803000000002 A0A1S3RA14_SALSA
#  G803000000003 A0A1S3RDQ3_SALSA

# UID_GN_MAP=os.path.expanduser('~/play/jones/uniprot-trembl-salmon.8030.map.tsv')
#  db  acc          uid                 gn
#  tr  A0A1S3RID5   A0A1S3RID5_SALSA    LOC106602976
#  tr  B5XFF4       B5XFF4_SALSA        WRB
#

UID_GN_MAP=os.path.expanduser('~/data/cococonet/human_uid_map.tsv')


OUTFILE=os.path.expanduser('~/play/jones/human_goa_results.tsv') 


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
    

def read_predout(predout, seqidmap):
    columns = ['seqid','goterm','prob'] 
    df = pd.read_csv(predout, sep='\t', header=None, names=columns)
    logging.debug(f"predout shape: {df.shape}")
    columns = ['seqid','uid'] 
    smdf = pd.read_csv(seqidmap, sep='\t', header=None, names=columns) 
    logging.debug(f"seqmap shape: {smdf.shape}")
    logging.debug(f"seqmap:\n{smdf}")

    fixedpredout = pd.merge(df, smdf,  how='left', on=['seqid'])
    fixedpredout.drop(['seqid'], inplace=True, axis=1)
    logging.debug(f"fixed pred out is \n{fixedpredout}")

        
    return fixedpredout

def fix_rowcol_names(network, mapfile):
    ugm = pd.read_csv(mapfile, sep='\t', header=0, index_col=0)
    logging.debug(f"uid_gn_map:\n{ugm}")
    #mapdict = pd.Series(ugm.uid.values, index=ugm.gn).to_dict()
    mapdict = pd.Series(ugm.gn.values, index=ugm.uid).to_dict()
    logging.debug(f"mapdict={mapdict}")
    gncolumns = list(network.columns)
    logging.debug(f"columnlist={gncolumns}")
    newcols = []
    for g in gncolumns:
        try:
            n = mapdict[g]
            logging.debug(f"got mapping {g} ->{n}")
            if pd.isna(n):
                newcols.append(g)
            else:
                newcols.append(n)
        except KeyError:
            logging.debug(f"mapping error with {g}")
            newcols.append(g)
    
    logging.debug(f"newcols={newcols[:10]} length={len(newcols)}")
    logging.debug(f"network shape={network.shape} assigning columns..")
    network.columns=newcols
    logging.debug("assigning row index..")
    network.index = newcols
    logging.debug("done.")
    return network
    


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)
    
    logging.info(f"Reading network: {HUMAN_NET}")
    nw = read_network_hdf5(HUMAN_NET)
    logging.info(f"network:\n{nw}")
    nw = fix_rowcol_names(nw, UID_GN_MAP)
    logging.info(f"fixed network:\n{nw}")        
    
    #logging.info(f"Reading predictions: {PREDOUT}")
    #po = read_predout(PREDOUT, SEQ_IDMAP)   
    #po.to_csv(f"{PREDOUT}.csv", sep="\t")
    #ogging.info(f"\n{po}")    
    
    #amdf = build_annotation_matrix(po, 'uid','goterm')
    #logging.info(f"\n{amdf}")    
    logging.debug(f"Reading in {HUMAN_GOA} ...")
    adf = pd.read_csv(HUMAN_GOA, sep=',', index_col=0)
       
    
    logging.info(f"input to run_egad: genesXgo:\n{adf}\ngenesXgenes:\n{nw}")    
    outdf = run_egad(adf, nw )
    logging.info(f"\n{outdf}")
    outdf.to_csv(f"{OUTFILE}", sep='\t')    
    #logging.info(f"Wrote to {OUTFILE}")