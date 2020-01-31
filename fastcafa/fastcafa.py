#!/usr/bin/env python
#   ************************CANONICAL COLUMNS*********************************
#   
# COLUMN    DESCRIPTION               MAPPINGS                         EXAMPLES
# cid       cafa4 target identifier   N/A                              T2870000001
# cafaprot  cafa4 target protein id                                    1433B
# cafaspec  cafa4 target protien species                               MOUSE
# cgid      cafa4 target gene id                                       LRRK2_MOUSE 
# pid/gid   UniProtKB: entry/ID                                        LRRK2_MOUSE 
# pacc      UniProtKB: accession  ?quickgo:gene product_id             P63103
# protein   all caps name                                              1433B
# gene      Free-text gene name.                                       Lrrk2  Ywahb   
# taxonid   NCBI taxon id                                              9606                 
# species   all caps code                                              MOUSE   PONAB
# goterms
# goterm    Gene Ontology Annotation ID                                GO:0005634
# goasp     biological process|molecular function|cellular component   bp       mf   cc
# goev      evidence code for GO annotation.                          IEA 
# eval      BLAST/HMMER/PHMMER expect statistic                        1.000000e-126
# bias      Adjustement to score for char prevalence                   3.5
# pscore    PHMMER bit-score                                           400.3
# db        database against which orthology query is done             sp (swissprot)
# score     General numeric score of a prediction. [ any scalar ]
# pest      Probability estimate for prediction.   [.01-1.0]           0.68    
# seq
# seqlen

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
from collections import defaultdict
from configparser import ConfigParser
import logging
import pickle
import pprint as pp
import random
import subprocess
import tempfile
import traceback

import pandas as pd
import numpy as np
np.set_printoptions(threshold=400)
from scipy import sparse
from sklearn import metrics

import h5py



GOASPECTMAP= { 'biological_process' : 'bp',
               'cellular_component' : 'cc',
               'molecular_function' : 'mf',
               'external'           : 'ex'
            }

UPASPECTMAP = { 'C': 'cc',
                'F': 'mf',
                'P': 'bp'
              }

# filled in by build_specmaps()
# SPECMAPS = None


class Ontology(dict):
    """
    gomatrix:  goterm x goterm np.ndarray fully propagated relationships. 
    gotermidx: { <str> goterm : <int> index in gomatrix, ... } 
    gotermlist:
    
    NOTE: dict key order is now assumed stable as of Python 3.7. If this is run on 
    earlier version, unstable key order will cause chaos. s 
    
    """
    
    instance = None
    
    def __init__(self, gomatrix, gotermidx, altidx):
        # a list of np.arrays
        self.data = gomatrix
        # a dictionary of goterms -> array_index
        self.gotermidx = gotermidx
        # a dictionary of (alternate) goterms -> (real) goterms
        self.altidx = altidx
        # list of go terms
        # **depends on key sort order stablility**
        self.gotermlist = list(gotermidx)
        
        
    def __getitem__(self, key):
        #
        try:
            val = self.data[self.gotermidx[key]]
        except KeyError:
            realgt = self.altidx[key]
            realidx = self.gotermidx[realgt]
            val = self.data[realidx]
        return val

    def __repr__(self):
        #return repr(self.__dict__)
        s = "{"
        for k in self.gotermidx.keys(): 
            s += f"'{k}' : {self.data[self.gotermidx[k]]} "
        s += "}"
        return s 

    def __len__(self):
        return len(self.data)

    def keys(self):
        return self.gotermidx.keys()
    
    def indexof(self, key):
        return self.gotermidx[key]

    def keyof(self, indexid ):
        return self.gotermlist[indexid]



class UniprotByGene(dict):
    """
    dictionary from geneid -> propagated boolean go vector for that geneid. 
    
    Assumes given Ontology object above, provided on construction.
    
    Provides easy contains to check if a given geneid is annotated 
    with a given goterm.  
    
    """
    
    instance = None
    
    def __init__(self, dict,  ontology):
        super(UniprotByGene, self).__init__(dict)
        self.ontobj = ontology
        
    def contains(self, gid, goterm):
        """
        Return boolean   True if goterm is annotated for that geneid  
        """
        gv = self[gid]
        #logging.debug(f"govector fo {gid} is {gv}")
        gtidx = self.ontobj.indexof(goterm)
        #logging.debug(f"goterm index for {goterm} is {gtidx}")
        return gv[gtidx]



def get_uniprot_bygene_object(config, usecache=True):
    goobj = get_ontology_object(config, usecache=True)
    ubgdict = get_uniprot_bygene(config, usecache=True)
    uobj = UniprotByGene(ubgdict, goobj )
    UniprotByGene.instance = uobj
    return uobj
  
        
def get_ontology_object(config, usecache=True):
    if Ontology.instance is None:
        build_ontology(config, usecache)
    return Ontology.instance

def do_build_ontology(config, usecache=False):
    """
    Rebuild all cached information.
    """
    logging.info("building ontology...")
    build_ontology(config, usecache)
    logging.info(f"made Ontology object: {Ontology.instance} ")
  


def do_build_uniprot(config, usecache=False):

    logging.info("getting uniprot/sprot info..")
    ubt = get_uniprot_byterm(config, usecache)
    logging.info(f"got ubt list:\n{ubt[0:20]}")


def get_uniprot_bygene(config, usecache=True):
    """
    loads from cache, or triggers build...
    
    
    """
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprotbygene.pickle"    
    ubgdict = None
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            ubgdict = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:
        
        ubgdict = build_uniprot_bygene(config)
        logging.debug(f"saving dict: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(ubgdict, cf )
            logging.debug(f"saved dict: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()        
    return ubgdict


def build_uniprot_bygene_pandas(config):
    """
    Since this is only used by do_evaluate we will only use experimentally
    validated annotations. 
    
    builds dictionary geneid -> propagated govector. 
    
    Returns dict for later use..
    """
    ontobj = get_ontology_object(config, usecache=True)
    ubt = get_uniprot_byterm_exponly_df(config, usecache=True)
    pids = list(ubt.pid.unique())
    logging.debug(f"pid list, length={len(pids)} e.g. : {pids[:50]}")
    bygenedict = {}
    sumreport = 1
    suminterval = 10000
    repthresh = sumreport * suminterval
    gtlength = len(ontobj.gotermidx)
    logging.debug(f"gv length is {gtlength}")
    
    for pid in pids:
        piddf = ubt[ubt.pid == pid]        
        current = np.zeros(gtlength, dtype=bool)
        # create vector of all correct goterms
        for (i, row) in piddf.iterrows():
            #logging.debug(f"goterm is {row.goterm}")
            gv = ontobj[row.goterm]
            current = current + gv
        #logging.debug(f"Constructed goterm vector for {pid} sum={current.sum()}")
        bygenedict[pid] = current
        if len(bygenedict) >= repthresh:
            logging.info(f"Processed {len(bygenedict)} entries... ")
            sumreport +=1
            repthresh = sumreport * suminterval
    
    logging.debug(f"Constructed dict by geneid {bygenedict[bygenedict.keys()[0]]}")
    return bygenedict

def build_uniprot_bygene(config):
    """
    Since this is only used by do_evaluate we will only use experimentally
    validated annotations. 
    
    builds dictionary geneid -> propagated govector. 
    
    Returns dict for later use..
    """
    ontobj = get_ontology_object(config, usecache=True)
    ubt = get_uniprot_byterm_exponly_df(config, usecache=True)
    pids = list(ubt.pid.unique())
    logging.debug(f"pid list, length={len(pids)} e.g. : {pids[:50]}")
    ubtd = ubt.to_dict(orient = 'index')
    logging.debug(f"converted DF to dict: e.g. \n{ [ubtd[i] for i in range(0,3)] } ")
   
    bygenedict = {}
    sumreport = 1
    suminterval = 10000
    repthresh = sumreport * suminterval
    gtlength = len(ontobj.gotermidx)
    logging.debug(f"gv length is {gtlength}")
    
    i = 0
    currentpid = None
    currentgv = np.zeros(gtlength, dtype=bool)
    while i < len(ubtd):
        row = ubtd[i]
        pid = row['pid']
        #logging.debug(f"row {i} is pid {pid}")
        if currentpid is None:
            currentpid = pid
            currentgv = currentgv + ontobj[row['goterm']]

        elif currentpid == pid:
            # same pid, continue...
            currentgv = currentgv + ontobj[row['goterm']] 
        else:
            # new pid
            bygenedict[currentpid] = currentgv 
            currentpid = pid
            currentgv = np.zeros(gtlength, dtype=bool)

        if len(bygenedict) >= repthresh:
            logging.info(f"Processed {len(bygenedict)} entries... ")
            sumreport +=1
            repthresh = sumreport * suminterval    
        
        i += 1
        
    samplekeys = list(bygenedict.keys())[:3]
    logging.debug(f"Made dict by geneid: {[ bygenedict[k] for k in samplekeys]} ")
    return bygenedict


def do_build_prior(config, usecache=True ):    

    logging.info("calculating prior..")
    priormatrix = calc_prior(config, usecache, species=None)
    logging.info(f"priormatrix = {matrix_info(priormatrix)}")
    priordf = get_prior_df(config, usecache)
    logging.info(f"priordf:\n{priordf}")    


def do_expression(config, infile, outfile, usecache=True):
    """
    perform phmmer query on each target 
    for best scoring ortholog hits, follow gene expression correlation
          annotate target with co-expressed genes annotations
          score is score of hit. 
    
    """

    logging.info("running expression")
    df = get_target_df(config, infile)
    logging.info(f"got target df:\n{df}")

    logging.info("making expression prediction...")
    df = calc_expression_prediction(config, df, usecache)
    logging.debug(f"prediction=\n{df}")#

    logging.info(f"writing to outfile {outfile}")
    df.to_csv(outfile)
    logging.info("done.")
    print(df)



def do_phmmer(config, infile, outfile, usecache=True):
    """
    Perform phmmer on infile sequences. 
    Output prediction to outfile.
    
    """
    
    logging.info("running phmmer")
    pdf = get_phmmer_df(config, infile)
    logging.info(f"got phmmer df:\n{pdf}")

    logging.info("making phmmer prediction...")
    df = calc_phmmer_prediction(config, pdf, usecache)
    logging.debug(f"prediction=\n{df}")

    logging.info(f"writing to outfile {outfile}")
    df.to_csv(outfile)
    logging.info("done.")
    print(df)


def parse_tfa_file( infile):
    """
    Reads .tfa file, determines species, target ids, geneids. 
    
    returns dataframe:
    cid  cgid
    
    """
    listoflists = []
       
    try:
        f = open(infile, 'r')
    except FileNotFoundError:
        logging.error(f"file not readable {filename} ")
    for line in f:
        # >T100900000004 1433G_MOUSE
        if line.startswith(">"):
            fields = line[1:].split()
            cid = fields[0].strip()
            cgid = fields[1].strip()
            listoflists.append( [cid, cgid] )
    logging.debug(f"got {len(listoflists)} cids with geneids.") 
    df = pd.DataFrame(listoflists, columns=['cid','cgid']) 
    return df    




def do_prior(config, infile, outfile, usecache=True):
    """
    Apply prior likelihood to all infile sequences. 
    Output prediction to outfile for later eval. 
    
    """
    logging.info("making prior prediction...")
    df = make_prior_prediction(config, infile)
    logging.debug(f"prediction=\n{df}")

    logging.info(f"writing to outfile {outfile}")
    df.to_csv(outfile)
    logging.info("done.")
    print(df)
  


def run_evaluate_pr(config, predictfile, outfile, goaspect=None):
    """
    Consume a prediction.csv file, and score based on accuracy. 
    X.prediction.csv
   
    """
    df = pd.read_csv(os.path.expanduser(predictfile), index_col=0)
    logging.debug(f"got predictdf types:\n{df.dtypes}\n{df}")
    edf = do_evaluate_pr(config, df, goaspect)
    logging.debug(f"got evaluation df:\n{edf}")
    edf.to_csv(outfile)
    
    max_goterms = config.get('global','max_goterms')
    eval_threshold = config.get('phmmer','eval_threshold')
    topx_threshold = config.get('phmmer','topx_threshold')
    score_method = config.get('phmmer','score_method')
    
    logging.info(f"hyperparams:\nmax_goterms={max_goterms}\neval_threshold={eval_threshold}\ntopx_threshold={topx_threshold}\nscore_method={score_method}  ")
    print(edf)



def is_correct_apply(row):
    """
    Function to check prediction for use with .apply and Pandas DF. 

    """
    return UniprotByGene.instance.contains(row.cgid, row.goterm)


def do_evaluate_pr(config, predictdf, goaspect):
    """
    Calculate precision recall number. 
    
    i    cid           goterm       score    cgid
    0    G960600000001 GO:0086041   53.0   Q9Y3Q4_HUMAN
    1    G960600000001 GO:0086089   49.0   Q9Y3Q4_HUMAN
    2    G960600000001 GO:0086090   49.0   Q9Y3Q4_HUMAN
    
    Return:
    max_goterms=1499
    eval_threshold=1.0e-120
    topx_threshold=200
    score_method=phmmer_score_weighted  
              cgid             cid  correct      goterm      pauc      pest     score     pr
        CHIA_MOUSE  G1009000000001     True  GO:0008150  0.759799  0.990000  0.046959  0.230651
        CHIA_MOUSE  G1009000000001    False  GO:0005575  0.759799  0.442188  0.020743  0.230651
        CHIA_MOUSE  G1009000000001    False  GO:0110165  0.759799  0.423416  0.019845  0.230651

    
    """
    logging.debug(f"got predictdf:\n{predictdf}")
    ubgo = get_uniprot_bygene_object(config, usecache=True)
    ontobj = get_ontology_object(config, usecache=True)
    logging.debug(f"got known uniprot and ontology object.")  
    
    outdf = pd.DataFrame(columns = ['cid','goterm','score','cgid','correct','pest','pr'])

    cidlist = list(predictdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")
    # Normalize all estimates from score to .01 - .99
    cmax = predictdf.score.max()
    cmin = predictdf.score.min()
    predictdf['pest'] = np.interp(predictdf['score'], (cmin, cmax ), (.01,.99))
    #logging.debug(f"cdf after score normalization -> pest is:\n{cdf.dtypes}\n{cdf}")
    
    for cid in cidlist:
        cdf = predictdf[predictdf.cid == cid].copy()
        # get gene id. 
        cgid = cdf.cgid.unique()[0]
        logging.debug(f"cgid is {cgid}")
        #logging.debug(f"geneid for this target is is {cgid}")
        cdf['correct'] = cdf.apply(is_correct_apply, axis=1)
        #cdf.reset_index(drop=True, inplace=True) 
        logging.debug(f"cdf after assessment:\n{cdf.dtypes}\n{cdf}")
        #logging.debug(f"cdf is:\n{cdf}")
        outdf = outdf.append(cdf, ignore_index=True)
    
    outdf['correct'] = outdf['correct'].astype(np.bool)
    logging.debug(f"outdf before pr is:\n{outdf}")
    outdf.sort_values(by='pest')
    pos = outdf[outdf.correct == True]
    logging.debug(f"positives is {pos}")
    poslist = list(pos.index.values)
    numtotal = len(predictdf)
    logging.debug(f"numtotal is {numtotal}")    
    logging.debug(f"len(poslist) is {len(poslist)}")

    pr = calc_precision_recall(poslist, numtotal)
    outdf['pr'] = pr

    logging.debug(f"outdf after pr is:\n{outdf}")
    return outdf

def calc_precision_recall(posidxlist, totalnum):
    """
        Calculates precision recall number...
        
    """
    n = len(posidxlist)
    sum = 0
    i = 1
    for posi in posidxlist:
        #val = ( i / posi + 1 )
        #logging.debug(f" {i} / {posi} = {val}")
        sum = sum + ( i / ( posi + 1 ) )
        logging.debug(f"sum is {sum}")
        i += 1

    pr = (1 / n) * sum 
    return pr


def do_evaluate_auroc(config, predictdf, goaspect):
    """
    i    cid           goterm       score    cgid
    0    G960600000001 GO:0086041   53.0   Q9Y3Q4_HUMAN
    1    G960600000001 GO:0086089   49.0   Q9Y3Q4_HUMAN
    2    G960600000001 GO:0086090   49.0   Q9Y3Q4_HUMAN
    
    Return:
    max_goterms=1499
    eval_threshold=1.0e-120
    topx_threshold=200
    score_method=phmmer_score_weighted  
              cgid             cid  correct      goterm      pauc      pest     score     auroc
        CHIA_MOUSE  G1009000000001     True  GO:0008150  0.759799  0.990000  0.046959  0.830652
        CHIA_MOUSE  G1009000000001    False  GO:0005575  0.759799  0.442188  0.020743  0.830652
        CHIA_MOUSE  G1009000000001    False  GO:0110165  0.759799  0.423416  0.019845  0.830652

    
    from sklearn import metrics
    metrics.auc(fpr, tpr)
    
    """
    
    logging.debug(f"got predictdf:\n{predictdf}")
    ubgo = get_uniprot_bygene_object(config, usecache=True)
    ontobj = get_ontology_object(config, usecache=True)
    logging.debug(f"got known uniprot and ontology object.")  
    
    outdf = pd.DataFrame(columns = ['cid','goterm','score','cgid','correct','pest','pauc'])

    cidlist = list(predictdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")
    for cid in cidlist:
        cdf = predictdf[predictdf.cid == cid].copy()
        # get gene id. 
        cgid = cdf.cgid.unique()[0]
        logging.debug(f"cgid is {cgid}")
        #logging.debug(f"geneid for this target is is {cgid}")
        cdf['correct'] = cdf.apply(is_correct_apply, axis=1)
        cdf.reset_index(drop=True, inplace=True) 
        logging.debug(f"cdf after assessment:\n{cdf.dtypes}\n{cdf}")
        # Normalize all estimates from score to .01 - .99
        cmax = cdf.score.max()
        cmin = cdf.score.min()
        cdf['pest'] = np.interp(cdf['score'], (cmin, cmax ), (.01,.99))
        logging.debug(f"cdf after score normalization -> pest is:\n{cdf.dtypes}\n{cdf}")
        try:                 
            pauc = metrics.roc_auc_score(cdf['correct'], cdf['pest'])
            logging.debug(f"pauc is {pauc}")
        except ValueError:
            pauc = .50
        
        cdf['pauc'] = pauc
        logging.debug(f"cdf is:\n{cdf}")
        outdf = outdf.append(cdf, ignore_index=True)
    outdf['correct'] = outdf['correct'].astype(np.bool)
    
    logging.debug(f"outdf before auc is:\n{outdf}")
    auroc = metrics.roc_auc_score(outdf['correct'], outdf['pest'])
    outdf['auroc'] = auroc
    #f1scr = metrics.f1_score(outdf['correct'], outdf['pest']         )
    #outdf['f1score'] = f1scr
    logging.debug(f"outdf after auroc is:\n{outdf}")
    return outdf


def get_evaluate_df(config, predictdf, goaspect=None,  threshold=None ):
    """
    lol:
         cid            numpredict   numcorrect  numannotated
    [['G960600000001',  4317,        539,        1637], 
     ['G960600000002',  10782,       0,          2], 
     ['G960600000003',  5634,        371,        555], 
    ]
    
    """
    
    df = do_evaluate_auroc(config, predictdf, goaspect)
    return df


def run_phmmer(config, filename):
    (outfile, exclude_list, cidcgidmap) = execute_phmmer(config, filename)
    logging.debug(f"cidcgidmap is {cidcgidmap}")
    phdict = parse_phmmer(config, outfile, exclude_list, cidcgidmap )
    return phdict

def get_phmmer_dict(config, filepath):
    logging.info("running phmmer")
    phdict = run_phmmer(config, filepath)
    logging.debug(f"got phmmer dict length: {len(phdict)}")
    return phdict

def get_phmmer_df(config, filepath):
    """
    orders by target,  evalue ascending (best fit first). 
    
    """
    phd = get_phmmer_dict(config, filepath)
    df = pd.DataFrame.from_dict(phd, orient='index')
    df = df.sort_values(by=['cid','eval'], ascending=[True, True]) 
    return df

def execute_phmmer(config, filename):
    """
    cpus = 8
    eval_threshold = 1.0e-120
    database=~/data/uniprot/uniprot_sprot.fasta
    remove_self_hits = True
    
    *Excludes geneids ( <protein>_<species> ) of sample from phmmer hit results*
    So the inbound sequence files *must* contain correct geneids (or some other
    method must be used to exclude self-hits). 

    
    """
    exclude_list = [] 
    cidgidmap = get_tfa_geneids(filename)
    for c in cidgidmap.keys():
        exclude_list.append(c)     
    outdir = config.get('global','outdir')
    outpath = os.path.dirname(filename)
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outfile = "%s/%s.phmmer.tbl.txt" % (outdir, filebase)
    cpus = config.get('phmmer','cpus')
    eval_threshold = config.get('phmmer','eval_threshold')
    database = config.get('phmmer','database')
    cmd = f"time phmmer --tblout {outfile} --noali --cpu {cpus} -E {eval_threshold} {filename} {database}"
    logging.debug(f"Running: {cmd}")
    cp = subprocess.run(cmd, 
                        shell=True, 
                        universal_newlines=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    logging.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd, outfile, cp.returncode))
    return (outfile, exclude_list, cidgidmap)


def get_tfa_geneids(filename):
    """
    Get the geneids of the TargetFiles in order to exclude these results from queries, e.g. 
    phmmer runs. This wouldn't necessarily cause a problem with the real files, because they 
    will be un-annotated. But for validation testing, the targets *will* be annotated (and
    need to be excluded). 
    
    Also need the mappings between cid and geneids for validation, where we need
    to look up the correct, annotated geneid to check predicitons. 
        
    """
    cidgids = {}
    
    try:
        f = open(filename, 'r')
    except FileNotFoundError:
        logging.error(f"file not readable {filename} ")
    for line in f:
        # >T100900000004 1433G_MOUSE
        if line.startswith(">"):
            fields = line[1:].split()
            cid = fields[0].strip()
            geneid = fields[1].strip()
            cidgids[cid] = geneid
    logging.debug(f"got {len(cidgids)} cids with geneids to exclude from results.")
    return cidgids
            

def parse_phmmer(config, filename, excludelist, cidcgidmap):
    """
    Read phmmer tbl out. Return dict  
    
    excludelist = [ '<geneids of inbound target proteins.>'] e.g. 2A5D_ARATH
     
    remove_self_hits = True|False
        shouldn't be necessary in standard case. 
    
    topx_threshold=<int>
        only include top X hits 
    
    """
    logging.info(f"Reading {filename}")
    df = pd.read_table(filename, 
                     names=['target','t-acc','cid','q-acc',
                            'eval', 'pscore', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
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
           
    topx = config.getint('phmmer','topx_threshold')
    
    if topx is not None:
        df = df.groupby('cid').head(topx).reset_index(drop=True) 
    
    dict = df.to_dict(orient='index')
    
    idxtodel = []
    for idx in dict.keys():
         (db,pacc,pid) = dict[idx]['target'].split('|')
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


def parse_hd5_expression(filename):
    """
    Loads data in file to dataframe.
    """
    with h5py.File(filename, 'r') as f:
        logging.debug("reading matrix...")
        matrix = f['agg'][:]
        logging.debug("reading rows...")
        rows = f['row'][:]
        logging.debug("reading columns...")
        columns = f['col'][:]
        logging.debug("making dataframe...")
        df = pd.DataFrame(matrix,  index=rows, columns = columns )
    return df    


def get_expression_dataset(config, species='YEAST'):
    """
    
    
    """
    
    

def calc_expression_prediction(config, dataframe, usecache):
    """
    
    """
    logging.debug("getting uniprot_byterm_df..")
    ubtdf = get_uniprot_byterm_df(config, usecache)
    ontobj = get_ontology_object(config, usecache)
    gtlength = len(ontobj.gotermidx)
    max_goterms = config.getint('global','max_goterms')


    return dataframe


def calc_phmmer_prediction(config, dataframe, usecache):
    """
    Takes phmmer df PDF: 
                cid           eval  pscore  bias    pacc          pid   cgid
    1     T100900000001  4.100000e-155  518.4   7.7  P35213    1433B_RAT   1A1L1_MOUSE
    2     T100900000001  5.400000e-155  518.0   7.2  A4K2U9  1433B_PONAB
    
    Gets uniprot_byterm_df UBTDF:
         pacc     species      goterm   goev    pid
    0        Q6GZX4   FRG3G     GO:0046782   IEA  001R_FRG3G
    1        Q6GZX3   FRG3G     GO:0033644   IEA  002L_FRG3G
    
    Algorithm:
        Determines predictive score by building goterm vector (term + parents) of all orthologs
        Score assigned is the score phmmer gave to hit. 
    
    Outputs prediction:
   
                   cid      goterm         score         cgid
    0     G1009000000001  GO:0008150  7.164505e+63   CHIA_MOUSE
    1     G1009000000001  GO:0003674  2.687680e+63   CHIA_MOUSE
    2     G1009000000001  GO:0008152  2.687680e+63   CHIA_MOUSE
   
    """
    logging.debug("getting uniprot_byterm_df..")
    ubtdf = get_uniprot_byterm_df(config, usecache)
    ontobj = get_ontology_object(config, usecache)
    gtlength = len(ontobj.gotermidx)
    max_goterms = config.getint('global','max_goterms')
    score_method = config.get('phmmer','score_method')

    pdf = dataframe
    cidlist = list(pdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")    
    gtarray = np.array(list(ontobj.gotermidx))
    
    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cid','goterm','score','cgid'])
    
    # Handle each cafa target in input list. 
    for cid in cidlist:
        #gv = np.zeros(gtlength)       
        
        cdf = pdf[pdf.cid == cid]
        cgid = cdf.reset_index().iloc[0].cgid       
        logging.debug(f"cgid for cid {cid} is {cgid}")
    
        #logging.debug(f"one cid df:\n{cdf}")
        
        # handle each ortholog protein, calculate score for inferred goterm set...
        for (i, row) in cdf.iterrows():
            gv = np.zeros(gtlength)       
            udf = ubtdf[ubtdf.pacc == row.pacc]
            pscore = row.pscore
            #           pacc   species   goterm     goev      pid
            # 2500682  O55201   MOUSE    GO:0032044  ISS  SPT5H_MOUSE
            # 2500683  O55201   MOUSE    GO:0005654  ISO  SPT5H_MOUSE
            # handle each goterm for this ortholog
            for (j, prow) in udf.iterrows():
                gv = gv + ontobj[prow.goterm].astype(np.int64)
                
            #logging.debug(f"sum is {gv.sum()} ")
            # we now have a govector with all goterms indicated by this ortholog.
            # each entry is sum of number of times that goterm was indicated (as annotated or
            # parent of an annotation term).
            # gv = array([123,   0,   3,   7, 345, 0])
        if score_method == 'phmmer_score':
            ones = np.ones(gtlength)
            gv = np.minimum(gv, ones)
            gv = gv * pscore
            #logging.debug(f"cdf is {cdf}")               
            #  cdf:
            #     cid           eval             pscore  bias    pacc          pid        cgid
            #74   G1009000000009  3.500000e-256  854.2  13.3  P33534   OPRK_MOUSE  OPRK_MOUSE
            #75   G1009000000009  5.900000e-254  846.9  13.3  P34975     OPRK_RAT  OPRK_MOUSE

        if score_method == 'phmmer_score_weighted':
            logging.debug(f"gv.dtype={gv.dtype} max={gv.max()} pscore={pscore}")
            gv = gv * pscore
            logging.debug(f"after gv.dtype={gv.dtype}") 
            

        #logging.debug(f"gv: {matrix_info(gv)}")
        gvnz = gv > 0
        #logging.debug(f"gvnz: {matrix_info(gvnz)}")
        #logging.debug(f"gtarray:length: {len(gtarray)} type:{type(gtarray)}")
        gotermar = gtarray[gvnz]
        govalar = gv[gvnz]
    
        cidar = np.full( len(govalar), fill_value=cid)        
        cgidar = np.full(len(govalar), fill_value=cgid)
        df = pd.DataFrame({ 'cid': cidar, 
                           'goterm': gotermar, 
                           'score' : govalar, 
                           'cgid' : cgidar })
        
        # This pulls out values, sorted by whatever 'score' is...
        df = df.nlargest(max_goterms, 'score')
        #df.sort_values(by='pest', ascending=False)
        #logging.debug(f"made dataframe for cid {cid}:\n{df}")
        topdf = topdf.append(df, ignore_index=True)
        
    logging.debug(f"made dataframe for all:\n{topdf}")
    return topdf


def make_prior_prediction(config, infile, species=None):
    """
     Same as calc_phmmer_prediction, but assigns prior likelihoods as score
    
    """
    cdf = parse_tfa_file(infile)
    logging.debug(f"Got cid/cgid frame:\n{cdf}") 
    pdf = get_prior_df(config, True, species)
    logging.debug(f"Got prior frame:\n{pdf}")     
    max_goterms = config.getint('global','max_goterms')
    logging.debug(f"max_goterms={max_goterms}")
    
    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cid','goterm','score','cgid'])
    # This pulls out values, sorted by whatever 'score' is...
    pdf = pdf.nlargest(max_goterms, 'pest')

    for (i, row) in cdf.iterrows():  
        #logging.debug(f"Row is:\n{row}")
        cid = row.cid
        cgid = row.cgid
        logging.debug(f"cid is {cid} cgid is {cgid}")
        df = pd.DataFrame(pdf,copy=True ) 
        #df['goterm'] = pdf['goterm']
        #df['score'] = pdf['pest']
        df['cid'] = cid
        df['cgid'] = cgid
        df.rename(columns = {'pest':'score'}, inplace = True)   
        logging.debug(f"made dataframe for cid {cid}:\n{df}")
        topdf = topdf.append(df, ignore_index=True)
    logging.debug(f"got out df types:\n{topdf.dtypes}\n{topdf}")
    return topdf

# outdf['correct'] = outdf['correct'].astype(np.bool)


def get_goprior(config, usecache, species=None):
    gp = calc_prior(config, usecache,species)
    return gp

def get_altiddict(config, usecache):
    if ALTIDDICT is not None:
        return ALTIDDICT
    else:
        build_ontology(config, usecache)    
        return ALTIDDICT


def build_ontology(config, usecache):
    """
    obofile=~/data/go/go.obo
    cachedir = ~/play/cafa4      
    
    from parse_obo():
    { 'GO:2001315': 
         {'is_a': ['GO:0009226', 'GO:0046349', 'GO:2001313'], 
         'part_of': [], 
         'goterm': 'GO:2001315', 
         'goname': 'UDP-4-deoxy-4-formamido-beta-L-arabinopyranose biosynthetic process', 
         'goasp': 'bp', 
       ...  
    }
 
    result:  Numpy boolean matrix of all relations in ontology, ordered by sorted goterm name.  
       
    """
    #global GOMATRIX
    #global GOTERMIDX
    #global ALTIDDICT
    #global GOTERMLIST
    # def __init__(self, gomatrix, gotermidx, altidx):
    
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('ontology','cachedir'))
    ontologyfile = f"{cachedir}/ontology.npy"
    termidxfile = f"{cachedir}/gotermindex.pickle"
    altiddictfile = f"{cachedir}/altiddict.pickle"
    include_partof = config.getboolean('ontology','include_partof')
    
    gomatrix = None
    
    if os.path.exists(ontologyfile) and usecache:
        logging.debug("Cache hit. Using existing matrix...")
        gomatrix = np.load(ontologyfile)
        logging.debug(f"loaded matrix: {matrix_info(gomatrix)} from {ontologyfile}")
        
        f = open(termidxfile, 'rb')
        gotermidx = pickle.load(f)
        f.close()

        f = open(altiddictfile, 'rb')
        altiddict = pickle.load(f)
        f.close()
                
        logging.debug(f"goterm index, e.g. : \n{list(gotermidx)[0]} :  {gotermidx[list(gotermidx)[0]]} ")
    
    else:
        (godict, altiddict) = parse_obo(config)
        
        # get keys from dict
        gotermlist = list(godict)
        logging.debug(f"parsed obo with {len(gotermlist)} entries. ")
        logging.debug(f"example entry:\n{gotermlist[0]}")
        logging.debug("sorting goterms")
        gotermlist.sort()
        logging.debug(f"sorted: e.g. {gotermlist[0:5]} ")
        logging.debug("creating goterm index dict.")
        #
        i = 0
        gotermidx = {}
        for gt in gotermlist:
            gotermidx[gt] = i
            i = i + 1
              
        logging.debug(f"creating zero matrix of dimension {len(gotermlist)}")
        shape = (len(gotermlist), len(gotermlist))
        gomatrix = np.zeros( shape, dtype=bool )
        logging.debug(f"filling in parent matrix for all goterms...")
        for gt in godict.keys():
            for parent in godict[gt]['is_a']:
                    gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        if include_partof:
            logging.debug("Including part_of relationships as is_a")
            for gt in godict.keys():
                for parent in godict[gt]['part_of']:
                        gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        
        #logging.debug(f"initial matrix:\n{print_square(gomatrix, GOTERMLIST)}")
        logging.debug("Calculating sparsity...")
        logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")
        logging.debug("converting to sparse matrix.")
        gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
        logging.debug(f"converging matrix: {matrix_info(gomatrix)}")
        gomatrix = converge_sparse(gomatrix)
        logging.info(f"got converged matrix:\n{matrix_info(gomatrix)} ")
        logging.debug(f"converged matrix sum={gomatrix.sum()}")
        #logging.debug("Calculating sparsity...")
        #sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
        #logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")        
        gomatrix = gomatrix.todense()
        gomatrix = np.asarray(gomatrix)
            
        logging.debug(f"Caching all values/indexes...")
        logging.debug(f"Saving matrix: {matrix_info(gomatrix)} to {ontologyfile}")
        np.save(ontologyfile, gomatrix)

        logging.debug(f"Saving gotermidx {len(gotermidx)} items to {termidxfile}")
        f = open(termidxfile, 'wb')   
        pickle.dump(gotermidx, f )
        f.close()
        
        logging.debug(f"Saving altiddict {len(altiddict)} items to {altiddictfile}.")
        f = open(altiddictfile, 'wb')
        pickle.dump(altiddict, f)
        f.close()
        
        logging.debug("Done constructing input for Ontology().")

    ontobj = Ontology(gomatrix, gotermidx, altiddict)
    # set global instance
    Ontology.instance = ontobj  
    logging.debug("Done creating Ontology object.")



def print_square(matrix, labels):
    """
    pretty prints square matrix with labels for debugging.
    """
    df = pd.DataFrame(matrix, columns=labels, index=labels)
    logging.debug(f"\n{df}")
    


def converge_sparse(matrix):
    logging.debug(f"starting matrix: \n{matrix_info(matrix)}")
    #logging.debug(f"{print_square(matrix.todense(), GOTERMLIST)}")
    oldval = 0
    logging.debug("Summing inbound matrix...")
    newval = matrix.sum()
    logging.debug("Beginning convergence loop.")
    icount = 0
    while oldval != newval:
        #logging.debug(f"Inbound matrix:\n{matrix_info(matrix)}")
        #logging.debug(f"oldval={oldval} newval={newval}")
        oldval = newval
        if not isinstance(matrix,  sparse.lil.lil_matrix): 
            logging.debug(f"{type(matrix)} is not scipy.sparse.lil.lil_matrix, converting... ")
            matrix = sparse.lil_matrix(matrix, dtype=np.bool)
        else:
            pass
            #logging.debug("matrix already lil_matrix...")
        #logging.debug("Multiplying...")
        mult = matrix @ matrix
        #logging.debug("Adding back original...")
        matrix = mult + matrix
        #logging.debug("Getting new sum...")
        newval = matrix.sum()
        #logging.debug(f"New matrix {icount}:\n{matrix.todense()}")
        #logging.debug(f"{print_square(matrix.todense(), GOTERMLIST)}")
        icount += 1
    logging.debug(f"Convergence complete. {icount} iterations. matrix:\n{matrix_info(matrix)}")
    return matrix

def matrix_info(matrix):
    """
    Returns string of inexpensive info on large matrix for logging. 
    """
    return f"type: {type(matrix)} shape: {matrix.shape} dtype: {matrix.dtype}"

def matrix_debug(matrix):
    """
    Returns string of somewhat more expensive info on large matrix for logging.
    
    """
    return f"type: {type(matrix)} shape: {matrix.shape} dtype: {matrix.dtype} sum: {matrix.sum()} min: {matrix.min()} max {matrix.max()}\n{matrix}"


def get_prior_df(config, usecache=True, species = None):

    ontobj = get_ontology_object(config, usecache=True)
    priormatrix = calc_prior(config, usecache, species)    
    priormatrix = priormatrix.transpose()
    logging.debug(f"priormatrix shape: {priormatrix.shape}")
    df = pd.DataFrame(priormatrix, index=ontobj.gotermlist, columns=['pest'])
    df.reset_index(inplace=True)
    map = {'index' : 'goterm'}
    df.rename(axis='columns', mapper=map, inplace=True)
    df = df.sort_values(by='pest', ascending=False)
    df.reset_index(drop=True, inplace=True)
    return df
   
    
def calc_prior(config, usecache, species=None):
    """
    take ontology matrix. goterms x goterms
    make total_termvector[ 17k ]
    
       0 proteinacc   1 species  2 goterm       3 goevidence
    0  '3AHDP'       'RUMGV'    'GO:0016491'    'IEA'
    1  '3AHDP'       'RUMGV'    'GO:0006694'    'TAS'

    for each protein:
        for each protein.goterm:
            gtvector = get_vector(goterm)
            total_termvector += gtvector
    
    if species is specified (as species code, e.g CAEEL), calculate prior for species only, 
    otherwise global
    
    ->
    
    vector of prob (.0001 - .99999)  [ 0.001, .0020, ... ] indexed by sorted gotermlist. 

    """ 
    #logging.debug(f"building sprot. species={species}")
    #sprot = get_uniprot_byterm(config, usecache=True)
    
    #df = pd.DataFrame(sprot,columns=['pacc','species','goterm','goev', 'pid'])
    #logging.debug(f"Built dataframe:\n{df}")    
    #exp_goev=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]
    #df = df[df['goev'].isin(exp_goev)] 
    #df.reset_index(drop=True, inplace=True)
    #sprot = df.values.tolist()
    #logging.debug(f"sprot, e.g.:\n{pp.pformat(sprot[0:5])} ... ")
   
    freqarray = None
    
    if species is None:
        filespecies = 'GLOBAL'
    else:
        filespecies = species
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprot.goprior.{filespecies}.npy"       
    
    if usecache and os.path.exists(cachefile):
        freqarray = np.load(cachefile)
        logging.debug(f"Loaded prior freqarray from file: {cachefile}")
    
    else:
        ontobj = get_ontology_object(config, usecache=True)
        gtlength = len(ontobj.gotermlist)
        logging.debug(f"building sprot. species={species}")
        sprot = get_uniprot_byterm(config, usecache=True)
        sdf = pd.DataFrame(sprot,columns=['pacc','species','goterm','goev', 'pid'])
        logging.debug(f"Built dataframe:\n{sdf}")    
        exp_goev=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]
        sdf = sdf[sdf['goev'].isin(exp_goev)] 
        sdf.reset_index(drop=True, inplace=True)
        sprot = sdf.values.tolist()
        
        if species is not None: 
            logging.debug(f"species {species} specified. converted to df with {df.shape[0]} rows: {df}")
            df = df[df.species == species ]
            logging.debug(f"removed other species. {df.shape[0]} rows left.")
                             
        sumarray = np.zeros(gtlength, dtype=np.int)
        logging.debug(f"made array {sumarray}")
        i = 0
        e = 0 
        elist = []
        logging.debug("summing over all terms and parents...")
        for r in sprot:
            gt = r[2]
            sumarray = sumarray + ontobj[gt]
            elist.append(gt)
            e += 1
        eset = set(elist)
        logging.debug(f"missing keys: {list(eset)}")
        logging.debug(f"added {i} gomatrix lines, with {len(eset)} distinct missing keys. ")
        logging.debug(f"sumarray: {sumarray} max={sumarray.max()} min={sumarray.min()} dtype={sumarray.dtype}")    
        divisor = sumarray.sum()
        logging.debug(f"divisor={divisor}")
        freqarray = sumarray / divisor
        np.save(cachefile, freqarray)
        logging.debug(f"Saved freqarray to {cachefile}")

    logging.debug(f"freqarray: {freqarray.dtype} {freqarray} max={freqarray.max()} min={freqarray.min()}")    
    return freqarray


def get_uniprot_df(config, usecache=True):
    lod = build_uniprot(config, usecache)
    df = pd.DataFrame(lod )
    logging.debug(f"Built dataframe:\n{df}")    
    return df    


def get_uniprot_byterm_df(config, usecache=True):
    lol = get_uniprot_byterm(config, usecache)
    df = pd.DataFrame(lol,columns=['pacc','species','goterm','goev', 'pid'])
    logging.debug(f"Built dataframe:\n{df}")    
    return df 

def get_uniprot_byterm_exponly_df(config, usecache=True):
    lol = get_uniprot_byterm(config, usecache)
    df = pd.DataFrame(lol,columns=['pacc','species','goterm','goev', 'pid'])
    logging.debug(f"Built dataframe:\n{df}")    
    exp_goev=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]
    df = df[df['goev'].isin(exp_goev)] 
    df.reset_index(drop=True, inplace=True)
    return df




def get_uniprot_byterm(config, usecache):
    """
    [ {'proteinid': '001R_FRG3G', 
       'protein': '001R', 
       'species': 'FRG3G', 
       'proteinacc': 'Q6GZX4', 
       'taxonid': '654924', 
       'goterms': {'GO:0046782': 'IEA'}, 
       'seqlength': 256, 
       'sequence': 'MAFSAEDVL......LYDDSFRKIYTDLGWKFTPL'},
       .
       .
       .
    ]
    
    
    Generate non-redundant uniprot with a row for every goterm:
       pacc         species    goterm       goev    pid
    0  '3AHDP'      'RUMGV'   'GO:0016491'  'IEA'
    1  '3AHDP'      'RUMGV'   'GO:0006694'  'TAS'
    
    
    """    
    lod = build_uniprot(config, usecache=True)
    ubt = []
    for p in lod:
        for gt in p['goterms'].keys():
            item = [ p['proteinacc'],
                     p['species'],
                     gt, 
                     p['goterms'][gt],
                     p['proteinid'], 
                    ]
            ubt.append(item)
    logging.debug(f"created uniprot_byterm with {len(ubt)} entries.")          
    return ubt          

def get_uniprot_testset(config, usecache, species, evidence=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]):
    """
    species  = 'MOUSE' 'HUMAN' -> taxonid matched. 
    evidence = list of OK codes | None means any
        Experimental: [ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ] q
    
    Generate DataFrame with data for export as CAFA test file:
    G<taxonid>_<datetime>_<5-digit-number>    <proteinid>   <sequence>
      
    """
    lodt = build_uniprot_test(config, usecache)
    #print(lodt)

    tdf = pd.DataFrame(lodt, columns=['pacc','protein','species','goterm','goev','seqlen','seq'])
    if evidence is not None:
        tdf = tdf[tdf['goev'].isin(evidence)] 
    tdf = tdf[tdf['species'] == species]
    tdf.reset_index(inplace=True)
    un = tdf.protein.unique()
    evcodes = tdf.goev.unique()
    logging.debug(f"generated {len(tdf)} item DF with {len(un)} proteins. evcode={evcodes}")
    return tdf



def do_testset(config, numseq, species, outfile):
    """
    Creates random target set from *annotated* swissprot proteins for given species. 
    Creates FASTA files exactly like CAFA TargetFiles
        
    :param    type      name:   desc
    :param    str       species          Species to generate: 
    :return   testfile  Path to test file generated in fasta format.
    :rtype
    :raise:
    
    """
    outfile = os.path.expanduser(outfile)
    numseq = int(numseq)
    logging.debug(f"numseq={numseq} species={species} outfile={outfile}")
    tdf = get_uniprot_testset(config, usecache=True, species=species)

    logging.info(f"got testset dataframe {len(tdf)} entries. ")
    
    #  [  { taxid     : speccode, ...  },
    #     { speccode  : taxonid, ...},
    #     { linnean   : taxonid, ...}   ]
    specmaps = get_specmaps(config)
    taxonid = specmaps[1][species.strip()]
    logging.debug(f"taxonid is {taxonid} for {species}")    
    up = tdf.protein.unique()
    upl = up.tolist()

    # throws error if 
    spl = random.sample(upl, numseq)
    snum = 1
    x = 60
    s = ""
    for p in spl:
        r = tdf[tdf.protein == p].reset_index().iloc[0]
        s += f">G{taxonid}{snum:08} {r.protein}_{r.species}\n"
        chunklist = [ r.seq[y-x:y] for y in range(x, len(r.seq)+x, x) ] 
        #logging.debug(f"chunklist {chunklist}")
        for c in chunklist:
            s += f"{c}\n"

        snum += 1
    #logging.debug(s)
    
    try:
        f = open(outfile, 'w')
        f.write(s)
        f.close()
        logging.debug(f"Wrote test data to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 

    return outfile



def build_uniprot_test(config, usecache):
    """
   
    [ {'proteinid': '001R_FRG3G', 
       'protein': '001R', 
       'species': 'FRG3G', 
       'proteinacc': 'Q6GZX4', 
       'taxonid': '654924', 
       'goterms': {'GO:0047043': 'IEA', 'GO:0006694': 'IEA'}, 
       'seqlength': 256, 
       'sequence': 'MAFSAEDVL......LYDDSFRKIYTDLGWKFTPL'},
       .
       .
    ]
   
    Create redundant dataframe for later slimming. Include sequence. Cache. 
   
    """    
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprottest.pickle"    
    lodt = None
    
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            lodt = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:    
        lod = build_uniprot(config, usecache=True)
       
    
        lodt = []
        for p in lod:
            newgts = {}
            for gt in p['goterms'].keys():
                evcode = p['goterms'][gt]
                item = [ p['proteinacc'],
                         p['protein'],
                         p['species'],
                         gt, 
                         evcode,
                         p['seqlength'],
                         p['sequence'] 
                      ]
                lodt.append(item)
                
        logging.debug(f"saving listofdicts: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(lodt, cf )
            logging.debug(f"saved listofdicts: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close() 
        logging.debug(f"created uniprot test source with {len(lodt)} entries.")
    return lodt


def build_prior(config, usecache, species=None):
    logging.debug(f"running calc_prior with species {species}")
    out = calc_prior(config, usecache, species=species)
    print(out)

def build_uniprot(config, usecache):
    """
    Builds list of dictionaries, each element is item in uniprot/sprot
    
    """
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprot.pickle"    
    listofdicts = None
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            listofdicts = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:
        
        listofdicts = parse_uniprot_dat(config)
        logging.debug(f"saving listofdicts: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(listofdicts, cf )
            logging.debug(f"saved listofdicts: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()        
    
    #df = pd.DataFrame(listofdicts)
    #logging.debug(f"Built dataframe {df}")
    return listofdicts


def build_specmaps(config, usecache):
    """
    builds three maps in form of list of dicts:
       [
         { taxid     : speccode, ...  },
         { speccode  : taxonid, ...},
         { linnean   : taxonid, ...}
       ]
    caches as pickle object. 
    
    """
    global SPECMAPS
        
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/specfile.pickle"
    specfile = os.path.expanduser(config.get('uniprot','specfile'))
    #     tax2spec  spec2tax  lin2tax
    map = [ {},     {},       {} ] 
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")
        try:
            cf = open(cachefile, 'rb')    
            map = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()
        
    else:
        entries = parse_speclist(config, specfile )
        # [ ('CALCT', 'E', '227173', 'Calidris canutus', 'Red knot'), ... ]
        for entry in entries:
            map[0][entry[2]] = entry[0]
            map[1][entry[0]] = entry[2]
            map[2][entry[3]] = entry[2]
                    
        logging.debug(f"saving map: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(map, cf )
            
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()
    
    #s = next(iter(map[0]))
    #logging.debug(f"map e.g. :\n{map} ")
    logging.info(f"saving map: to {cachefile}")
    SPECMAPS = map
    return map

def get_specmaps(config):
    sm = build_specmaps(config, usecache=True)
    return sm



def parse_uniprot_dat(config):
        """
        Parses uniprot/sprot DAT file, returns as list of dicts, with sub-dicts...
        [ {'proteinid': '4CLL9_ARATH', 
           'protein': '4CLL9', 
           'species': 'ARATH', 
           'proteinacc': 'Q84P23', 
           'taxonid': '3702', 
           'goterms': {'GO:0005777': ['cc', 'IDA'], 
                       'GO:0005524': ['mf', 'IEA'], 
                       'GO:0004321': ['mf', 'IDA'], 
                       'GO:0016874': ['mf', 'IEA'], 
                       'GO:0009695': ['bp', 'IDA'], 
                       'GO:0031408': ['bp', 'IEA']}
           },
           .
           .
           .
    
        """
        logging.debug(f"using config with sections: {config.sections()}...")
        try:
            filepath = os.path.expanduser(config.get('uniprot','datfile'))
            logging.debug(f" attempting to open '{filepath}'")
            filehandle = open(filepath, 'r')
        except FileNotFoundError:
            logging.error(f"No such file {filepath}")                
        
        allentries = []
        current = None
        sumreport = 1
        suminterval = 10000
        repthresh = sumreport * suminterval
        try:
            while True:
                line = filehandle.readline()
                if line == '':
                    break

                if line.startswith("ID "):
                    # ID   001R_FRG3G              Reviewed;         256 AA.
                    #      <prot_name>_<prot_spec>
                    proteinid = line[5:16].strip()
                    current = defaultdict(dict)
                    current['proteinid'] = proteinid
                    (protein, species) = proteinid.split('_')
                    current['protein'] = protein
                    current['species'] = species
                    #logging.debug("Handling ID. New entry.")                
                
                elif line.startswith("AC "):
                    # AC   Q6GZX4;
                    # AC   Q91896; O57469;
                    #logging.debug("Handling AC.")
                    accession = line[5:11].strip()
                    current['proteinacc'] = accession

                elif line.startswith("OX   "):
                    #OX   NCBI_TaxID=654924;
                    #logging.debug("Handling OX.")
                    taxonid = ""
                    val = line[5:]
                    fields = val.split('=')
                    if fields[0] == 'NCBI_TaxID':
                        taxonid = fields[1].strip().replace(';','')
                    current['taxonid'] = taxonid
                    
                elif line.startswith("DR   GO;"):
                    # DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
                    # P biological process, C cellular component, F molecular function.  
                    #logging.debug("Handling DR.")
                    fields = line.split(';')
                    goterm = fields[1].strip()
                    goinfo = fields[2]
                    aspcode = goinfo.split(':')[0].strip()
                    goaspect = UPASPECTMAP[aspcode]
                    goevsrc = fields[3]
                    (goevidence, evsrc) = goevsrc.split(':') 
                    goevidence = goevidence.strip()
                    current['goterms'][goterm] = goevidence

                elif line.startswith("SQ   SEQUENCE"):
                    #logging.debug("Handling SQ:  XXX")
                    # line = filehandle.readline()
                    current['seqlength'] = int(line.split()[2])
                    current['sequence'] = ""
                    seqlen = current['seqlength']
                    aaread = 0
                    while aaread < seqlen:
                        line = filehandle.readline()
                        lineseq = line.strip().replace(" ","")
                        current['sequence'] = "%s%s" % (current['sequence'], lineseq)
                        aaread += len(lineseq) 

                elif line.startswith("GN   "):
                    # Examples:
                    #  GN   ABL1 {ECO:0000303|PubMed:21546455},
                    #  GN   Name=BRCA1; Synonyms=RNF53;
                    #  GN   ORFNames=T13E15.24/T13E15.23, T14P1.25/T14P1.24;
                    #logging.debug("Handling GN.")
                    val = line[5:]
                    
            
                elif line.startswith("//"):
                    #logging.debug("End of entry.")                  
                    allentries.append(current)
                    #logging.debug(f"All entries list now {len(allentries)} items... ")
                    if len(allentries) >= repthresh:
                        logging.info(f"Processed {len(allentries)} entries... ")
                        sumreport +=1
                        repthresh = sumreport * suminterval
                    current = None
                
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        if filehandle is not None:
            filehandle.close()      
        
        logging.info(f"Parsed file with {len(allentries)} entries" )
        logging.debug(f"First entry:  {allentries[0]}")
        return allentries


def parse_obo(config):
    """
    creates dict of dicts. key is goterm, contents is dict of 
       
       goterm  ""
       goname  ""
       goasp   ""
       godef   ""
       goisa   [] of goterms
       gohasa  [] of goterms
    
    
    """
    obofile = os.path.expanduser(config.get('ontology','obofile'))
    filehandle = open(obofile)
    godict = {}
    altids = {}
    current = None
    logging.info(f"Parsing file {obofile}")
    try:
        for line in filehandle:
            if line.startswith("[Typedef]"):
                godict[current['goterm']]= current
                break
            elif line.startswith("[Term]"):     
                if current is not None:
                    godict[current['goterm']]= current
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
                
                current['goasp'] = GOASPECTMAP[asp]

            # must create a separate mapping from alt_ids that maps to
            # primary, so it can be added to gomatrix properly
            elif line.startswith("alt_id: "):
                #current['alt_id'] = line[8:18].strip()
                #current['goasp'] = GOASPECTMAP[asp]            
                altid = line[8:18].strip()
                altids[altid] = current['goterm'] 
            
            #elif line.startswith("def: "):
            #    current['godef'] = line[5:].strip()

            #elif line.startswith("synonym: "):
            #    current.synonym.append(line[9:].strip())

            elif line.startswith("is_a: "):
                current['is_a'].append(line[6:16].strip())
            
            elif line.startswith("relationship"):
                if "part_of" in line:
                    current['part_of'].append(line[22:32])
                                 
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    logging.info(f"Parsed file with {len(godict)} terms and {len(altids)} alt terms.")    
    return (godict, altids)


def parse_speclist(config, filepath):
    '''
    Parses uniprot speclist.txt    https://www.uniprot.org/docs/speclist.txt
    to local .CSV
    
    taxonid   species   lineanname       commonname
    72259      ABANI    Abaeis nicippe   Sleepy orange butterfly
                                         
    OXYMO E  475340: N=Oxytenis modestia
                     C=Costa Rica leaf moth
                     S=Dead-leaf moth

    returns 
    
    '''
    logging.debug("Opening species map file %s" % filepath)

    
    try:
        fh = open(filepath, 'r')
    except FileNotFoundError:
        logging.error("No such file %s" % filename)                
   
    species = None
    kingdom = None
    taxonid = None
    lineanname = None
    commonname = None
    
    columnnames = ['species','kingdom','taxonid','lineanname','commonname']
    datalist = []
    # list of tuples
          
    try:
        for line in fh:
            #logging.debug("handling line %s" % line)
            if 'N=' in line and not line.startswith('Code')  :
                #logging.debug("handling N= line. taxonid is %s" % taxonid)
                if species is not None:
                    tup = (species, kingdom, taxonid, lineanname, commonname)
                    #logging.debug("Adding tuple: %s" % str(tup))
                    datalist.append( tup )
                    # reset all varaiables
                    species = kingdom = taxonid = lineanname = commonname = None
                species = line[:5].strip()
                kingdom = line[6]
                taxonid = line[7:15].strip()
                lineanname = line[19:].strip()
                #logging.debug("handling N= line. taxonid is %s" % taxonid)         
            elif 'C=' in line :
                commonname = line[19:].strip()
            elif 'S=' in line :
                 pass
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    finally:
        fh.close()    
    logging.debug("Parsed file with %d terms" % len(datalist) )
    return datalist


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cshl-work/etc/fastcafa.conf"))
    return cp


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


    parser.add_argument('-a', '--aspect',
                        action="store", 
                        dest='goaspect', 
                        default=None,
                        help='Generate/score only for specific GO aspect. [bp|cc|mf] ')

    parser.add_argument('-C', '--usecache',
                        action='store_true', 
                        dest='nocache',
                        default=False, 
                        help='Use cached information.' )
    
    
    subparsers = parser.add_subparsers( dest='subcommand',
                                        help='sub-command help.')

################################ run phmmer ########################################
   
    parser_phmmer = subparsers.add_parser('phmmer',
                                          help='run phmmer and output prediction')
    
    parser_phmmer.add_argument('-n','--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence files')

    parser_phmmer.add_argument('-o','--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')

################################ run prior ########################################

    parser_prior = subparsers.add_parser('prior',
                                          help='calculate prior and output prediction')
    
    parser_prior.add_argument('-n','--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence files')

    parser_prior.add_argument('-o','--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')

    parser_prior.add_argument('-s','--species', 
                               metavar='species', 
                               type=str,
                               default=None, 
                               help='species-specific prior, otherwise global')

######################### build and cache ontology ####################################

    parser_buildontology = subparsers.add_parser('build_ontology',
                                          help='build and cache GO ontology')    

    parser_builduniprot = subparsers.add_parser('build_uniprot',
                                          help='build and cache uniprot info')

######################### build and cache prior info ####################################
    
    parser_buildprior = subparsers.add_parser('build_prior',
                                          help='build and cache uniprot prior')
    
    parser_buildprior.add_argument('-s', '--species',
                        action="store", 
                        dest='species', 
                        default=None,
                        help='Limit to species where relevant.')

######################### build and cache species maps ####################################

    parser_buildspecies = subparsers.add_parser('build_species',
                                          help='build and cache NCBI species maps')             

######################### evaluate prediction for known ####################################    
        
    parser_evaluate = subparsers.add_parser('evaluate',
                                          help='evaluate prediction against known. output stats.')
    
    parser_evaluate.add_argument('-p', '--predictcsv', 
                               metavar='predictcsv', 
                               type=str, 
                               help='a .csv prediction file')        

    parser_evaluate.add_argument('-a', '--aspect', 
                               metavar='goaspect', 
                               type=str,
                               default=None, 
                               help='GO aspect to limit evaluation to.') 
          
    parser_evaluate.add_argument('-o', '--outcsv', 
                               metavar='outcsv', 
                               type=str, 
                               help='a .csv output file with stats')    


    parser_builduniprot_test = subparsers.add_parser('build_uniprot_test',
                                          help='build and cache uniprot test source info')

    parser_testset = subparsers.add_parser('testset',
                                          help='generate a set of input .TFA files from species')

    parser_testset.add_argument('-n','--numseq', 
                               metavar='numseq',
                               dest='numseq', 
                               default=10,
                               type=int, 
                               help='Number of sequences to select.')   
    
    parser_testset.add_argument('-s', '--species', 
                               metavar='species',                                
                               dest='species', 
                               required=True,
                               type=str, 
                               help='a .fasta sequence file with exp. annotated proteins')   
    
    parser_testset.add_argument('-o', '--out',  
                               metavar='outfile',
                               dest='outfile',
                               required=True, 
                               type=str, 
                               help='a DF .csv prediction file')
    

    args= parser.parse_args()
    
    # default to INFO
    logging.getLogger().setLevel(logging.INFO)
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
    
    logging.debug(f"args: {args}")
    
    if args.subcommand == 'phmmer':
        do_phmmer(cp, args.infile, args.outfile, usecache=True )
    
    if args.subcommand == 'expression':
        do_expression(cp, args.infile, args.outfile, usecache=True )
    
    if args.subcommand == 'prior':
        do_prior(cp, args.infile, args.outfile, usecache=True)
       
    if args.subcommand == 'build_ontology':
        build_ontology(cp, usecache=False)

    if args.subcommand == 'build_uniprot':
        build_uniprot(cp, usecache=False)

    if args.subcommand == 'build_species':
        build_specmaps(cp,  usecache=False)

    if args.subcommand == 'build_prior':
        build_prior(cp, args.species )

    if args.subcommand == 'build_uniprot_test':
        build_uniprot_test(cp, usecache=False )

    if args.subcommand == 'testset':
        do_testset(cp, args.numseq, args.species, args.outfile )

    if args.subcommand == 'evaluate':
        run_evaluate_pr(cp, args.predictcsv, args.outcsv, args.goaspect)
    