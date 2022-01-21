#!/usr/bin/env python

import argparse
import logging
import os
import subprocess
import sys

from configparser import ConfigParser

import pandas as pd

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)
from utils import *

def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cshlwork/etc/phmmer.conf"))
    return cp


def execute_phmmer(config, queryfile, database=None):
    """
    cpus = 8
    eval_threshold = 1.0e-120
    database=~/data/uniprot/uniprot_sprot.fasta
    remove_self_hits = True

    First item in header (after '>') in queryfile must be the same identifier 
    after the '>' in the database, in order to exclude self-hits. 


    """
    logging.debug(f"executing with filename={queryfile} database={database} ") 
    outdir = os.path.expanduser(config.get('phmmer','cachedir'))
    filename =os.path.expanduser(queryfile)
    outpath = os.path.dirname(queryfile)
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outfile = "%s/%s.phmmer.tbl.txt" % (outdir, filebase)
    logging.debug(f"outfile={outfile}")
    cpus = config.get('phmmer','cpus')
    eval_threshold = config.get('phmmer','eval_threshold')
    score_threshold = config.get('phmmer','score_threshold')
    if database is None:
        database = os.path.expanduser(config.get('phmmer', 'database'))
        logging.debug(f"Using config file database: {database}")
    
    cmd = ["/usr/bin/which","phmmer"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readlines()
    logging.debug(f'which phmmer gave {res}')

    cmd = [ 'phmmer',
           '--tblout', outfile , 
           '--noali',
           '--cpu', cpus,
           #'-E', eval_threshold,
           '-T', score_threshold,
           queryfile,
           database 
           ]
    
    logging.debug(f"Running: {cmd}")
    logging.info(f"Command line: {' '.join(cmd)} ")
    cp = subprocess.run(cmd)
   
    logging.debug(f"Ran cmd='{cmd}' outfile={outfile} returncode={cp.returncode} " )
    logging.debug(f"returning outfile={outfile} ")
    return outfile

#def parse_phmmer(config, filename, excludelist, cidcgidmap):
def parse_phmmer(config, filename):
    """
    Read phmmer tbl out. Return dict  
    
    """
    logging.info(f"Reading {filename}")
    df = pd.read_table(filename, 
                     names=['target','t-acc','query','q-acc',
                            'eval', 'score', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                            'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                     skip_blank_lines=True,
                     comment='#',
                     index_col=False,
                     skiprows=3,
                     engine='python', 
                     sep='\s+')
    
    
    logging.debug("Dropping unneeded columns..")
    df = df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 
             'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 
             'description'] , axis=1)
    logging.debug(f"sorting by query and score...")
    df = df.sort_values(by=['query','score'], ascending=[True, False])
    logging.debug(f"got phmmer df:\n{df}")
           
    topx = config.getint('phmmer','topx_threshold')
    logging.debug(f"topx ={topx}")
        
    if topx is not None:
        df = df.groupby('query').head(topx).reset_index(drop=True) 
    
    dict = df.to_dict(orient='index')
    logging.debug(f"dict is {dict}")
    idxtodel = []
    
    
    if logging.root.level >= logging.DEBUG: 
        s ="{"
        for k in dict.keys():
            s+=f"{k} : {dict[k]}\n"
        s += "}"
        logging.debug(f"{s}")
    logging.debug(f"returning OK. dict={dict}")
    return dict



'''    for idx in dict.keys():      
        (db,pacc,pid) = dict[idx]['target'].split('|')
        logging.debug(f"split |-separated target field...") 
        cid = dict[idx]['cid']
        # only exclude hits for the exact target protein...
        if pid == cidcgidmap[cid]:
            logging.debug(f"Found the pid {pid} excluded to be excluded for this cid {cid}")
            idxtodel.append(idx)
        else:
            (protein, species) = pid.split('_')
            dict[idx]['pacc'] = pacc
            dict[idx]['pid'] = pid
            dict[idx]['cgid'] = cidcgidmap[dict[idx]['cid']]
            del dict[idx]['target']
    for idx in idxtodel:
        del dict[idx]
'''


def get_phmmer_df(config, queryfile, database=None ):
    """
    orders by target, evalue ascending (best fit first).  
    cache output for later usage....

    """
    logging.debug(f"Handling input file {queryfile}")
    infile = os.path.expanduser(queryfile)
    pcachedir = config.get('phmmer','cachedir') 

    filename = os.path.basename(os.path.expanduser(queryfile))
    (filebase, e) = os.path.splitext(filename)
    pcachefile =f"{pcachedir}/{filebase}.phmmerdf.csv"
    
    df = None
    try:
        df = pd.read_csv(pcachefile, index_col=0)
        logging.info("Cached phmmer run found. ")
        
    except FileNotFoundError:
        logging.info("No cached phmmer data. Running...")      
        outfile = execute_phmmer(config, queryfile, database)
        logging.debug(f"outfile is {outfile} parsing...")
        phd = parse_phmmer(config, outfile)
        
        if len(phd) > 0:
            df = pd.DataFrame.from_dict(phd, orient='index')
            df = df.sort_values(by=['query','score'], ascending=[True, False])
        if df is not None:
            logging.debug(f"Caching phmmer output to {pcachefile}")
            df.to_csv(pcachefile)
    return df


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
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
                        help='Query .tfa file. ') 

    parser.add_argument('-D','--dbasefile', 
                        metavar='dbasefile', 
                        type=str,
                        required=False,
                        default=None, 
                        help='Fasta .tfa search database.')  
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
    
    c = get_default_config()
    if args.dbasefile is not None:
        dbfile = os.path.expanduser(args.dbasefile)
    else:
        dbfile = None
        
    queryfile =  os.path.expanduser(args.infile)
    logging.debug(f"queryfile={queryfile} dbfile={dbfile}")
    df = get_phmmer_df(c, queryfile, dbfile  )
    logging.info(f'phmmer outfile={df}')
