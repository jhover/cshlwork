#!/usr/bin/env python
#   ************************CANONICAL COLUMNS*********************************
#   
# COLUMN    DESCRIPTION               MAPPINGS                         EXAMPLES
# cid       cafa4 target identifier   N/A                              T2870000001
# cafaprot  cafa4 target protein id                                    1433B
# cafaspec  cafa4 target protien species                               MOUSE
# id        UniProtKB: entry/ID                                        1433B_RAT
# pacc      UniProtKB: accession  ?quickgo:gene product_id             P63103
# protein   all caps name                                              1433B
# gene      Free-text gene name.                                       Lrrk2  Ywahb
# geneid    Gene name+species.                                         LRRK2_MOUSE     
# taxonid   NCBI taxon id                                              9606                 
# species   all caps code                                              MOUSE   PONAB
# goterm    Gene Ontology Annotation ID                                GO:0005634
# goasp     biological process|molecular function|cellular component   bp       mf   cc
# goev      evidence codes for GO annotation.                          IEA 
# eval      BLAST/HMMER/PHMMER expect statistic                        1.000000e-126
# bias      Adjustement to score for char prevalence                   3.5
# score     BLAST/HMMER/PHMMER bit-score                               400.3
# db        database against which orthology query is done             sp (swissprot)
# probest   Probability estimate for prediction.                       0.68  [.01-1.0]  
#
#
#
__author__ = "John Hover"
__copyright__ = "2019 John Hover"
__credits__ = []
__license__ = "Apache 2.0"
__version__ = "0.99"
__maintainer__ = "John Hover"
__email__ = "hover@cshl.edu"
__status__ = "Testing"

import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

import argparse
from configparser import ConfigParser
import logging
import subprocess
import tempfile
import traceback


import pandas as pd
import numpy as np


GONSMAP= { 'biological_process' : 'bp',
             'cellular_component' : 'cc',
             'molecular_function' : 'mf',
             'external'           : 'ex'
            }


def dorun(config, filename, runname):
    logging.info("starting...")
    outfile = runphmmer(config, filename)
    phdict = parsephmmer(config, outfile )
    logging.debug(phdict)
    ont = build_ontology(config)

def runphmmer(config, filename):
    """
    cpus = 8
    eval_threshold = 1.0e-120
    database=~/data/uniprot/uniprot_sprot.fasta
    """
    outdir = config.get('global','outdir')
    outpath = os.path.dirname(filename)
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outfile = "%s/%s.phmmer.tbl.txt" % (outdir, filebase)
    cpus = config.get('phmmer','cpus')
    eval_threshold = config.get('phmmer','eval_threshold')
    database = config.get('phmmer','database')
    cmd = f"time phmmer --tblout {outfile} --noali --cpu {cpus} -E {eval_threshold} {filename} {database}"
    cp = subprocess.run(cmd, 
                        shell=True, 
                        universal_newlines=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    logging.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd, outfile, cp.returncode))
    return outfile

def parsephmmer(config, filename):
    """
    Read phmmer tbl out. Return dict  
    """
    logging.info("Reading %s" % filename)
    df = pd.read_table(filename, 
                     names=['target','t-acc','cid','q-acc',
                            'eval', 'score', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                            'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                     skip_blank_lines=True,
                     comment='#',
                     index_col=False,
                     skiprows=3,
                     engine='python', 
                     sep='\s+')
    logging.debug(str(df))
    logging.debug("Dropping unneeded columns..")
    df = df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 
             'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 
             'description'] , axis=1)
    dict = df.to_dict(orient='index')
    for idx in dict.keys():
         (db,pacc,pid) = dict[idx]['target'].split('|')
         (protein, species) = pid.split('_')
         dict[idx]['pacc'] = pacc
         dict[idx]['pid'] = pid
         del dict[idx]['target']
    if logging.root.level >= logging.DEBUG: 
        s ="{"
        for k in dict.keys():
            s+=f"{k} : {dict[k]}\n"
        s += "}"
        logging.debug(f"{s}")
    """
    { 0 : {'cid': 'T100900000001', 'eval': 1.0999999999999997e-156, 'score': 523.6, 'bias': 8.5, 'pacc': 'Q9CQV8', 'pid': '1433B_MOUSE'}
      1 : {'cid': 'T100900000001', 'eval': 4.099999999999999e-155, 'score': 518.4, 'bias': 7.7, 'pacc': 'P35213', 'pid': '1433B_RAT'}
    }
    
    """
    return dict


def build_ontology(config):
    """
    obofile=~/data/go/go.obo
    
    """
    dict = parse_obo(config)
    #for k in dict.keys():
    #    po = dict[k]['part_of']
    #    print(po)
    logging.debug(f"got dict: {dict}")
    
    

def calc_prior(config,species=None):
    """
    
    
    """
    

def parse_obo(config):
    """
    creates dict of dicts. key is goterm, contents is dict of 
       goterm  ""
       goname  ""
       goasp   ""
       godef   ""
       goisa   [] of goterms
       gohasa  [] 
    
    """
    obofile = os.path.expanduser(config.get('ontology','obofile'))
    filehandle = open(obofile)
    goidx = {}
    current = None
    logging.info(f"Parsing file {obofile}")
    try:
        for line in filehandle:
            if line.startswith("[Term]"):     
                if current is not None:
                    goidx[current['goterm']]= current
                # create new item...
                current = {}
                current['is_a'] = []
                current['part_of'] = []
                
            elif line.startswith("id: "):
                current['goterm'] = line[4:].strip()
                
            elif line.startswith("name: "):
                current['goname'] = line[6:].strip()
            
            elif line.startswith("namespace: "):
                asp = line[11:].strip()
                current['goasp'] = GONSMAP[asp]
            
            elif line.startswith("def: "):
                current['godef'] = line[5:].strip()

            #elif line.startswith("synonym: "):
            #    current.synonym.append(line[9:].strip())

            elif line.startswith("is_a: "):
                current['is_a'].append(line[6:16].strip())
            
            elif line.startswith("relationship"):
                if "part_of" in line:
                    current['part_of'].append(line[22:32])
                                 
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    logging.info(f"Parsed file with {len(goidx)} terms")    
    return goidx




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
                        help='a .fasta sequence files')
    
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default='~/etc/cafa4.conf',
                        help='Config file path [~/etc/cafa4.conf]')

    parser.add_argument('-n', '--name', 
                        action="store", 
                        dest='runname', 
                        default='default',
                        help='Run-specific identifier to use in file output.')    
                    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
           
    dorun(cp, args.infile, args.runname)

