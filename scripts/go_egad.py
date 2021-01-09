#! /usr/bin/env python
import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)
from utils.cococonet import *

gitpath=os.path.expanduser("~/git/pyEGAD")
sys.path.append(gitpath)
from egad.egad import *

SPECIES_NET=os.path.expanduser('~/data/cococonet/networks/afrog_prioAggNet.hdf5')
PREDOUT=os.path.expanduser('~/work/jones/afrog_prio.predout')
OUTFILE=os.path.expanduser('~/work/jones/afrog_prio.egad.tsv')


def read_predout(predout):
    columns = ['seqid','goterm','prob'] 
    df = pd.read_csv(predout, sep='\t', header=None, names=columns)
    logging.debug(f"predout\n{df.head()}")
    logging.debug(f"predout shape: {df.shape}")
    return df

def do_egad(species_net=SPECIES_NET, go_pred=PREDOUT ):
    logging.info(f"Reading network: {species_net}")
    nw = read_network_hdf5(species_net)
    
    logging.info(f"Reading predictions: {go_pred}")
    po = read_predout(go_pred)   
    logging.info(f"\n{po}")    
    
    amdf = build_annotation_matrix(po, 'seqid','goterm')
    logging.info(f"\n{amdf}")    
    
    logging.info(f"input to run_egad: genesXgo:\n{po}\ngenesXgenes:\n{nw}")    
    outdf = run_egad(amdf, nw )
    logging.info(f"\n{outdf}")
    return outdf

def write_outfile(df, outfile=OUTFILE): 
    outdf.to_csv(f"{outfile}", sep='\t')    
    logging.info(f"Wrote to {outfile}")


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('network', 
                        metavar='network', 
                        type=str, 
                        help='A cococonet hd5 network')

    parser.add_argument('goterms', 
                        metavar='goterms', 
                        type=str, 
                        help='a GO prediction file.')

    parser.add_argument('outfile', 
                        metavar='outfile', 
                        type=str, 
                        help='EGAD output file.')    
    
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
        
    outdf = do_egad(args.network, args.goterms)
    write_outfile(outdf, args.outfile)
    
