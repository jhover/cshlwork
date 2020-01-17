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
# pest      Probability estimate for prediction.                       0.68  [.01-1.0]  
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
from collections import defaultdict
from configparser import ConfigParser
import logging
import pickle
import pprint as pp
import subprocess
import tempfile
import traceback

import pandas as pd
import numpy as np
np.set_printoptions(threshold=400)
from scipy import sparse

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
SPECMAPS = None

# filled in by build_ontology()
GOTERMIDX = None
ALTIDDICT = None
GOMATRIX = None
GOTERMLIST = None

# filled in by calc_prior()
GOPRIOR = None


def dorun(config, filename, runname, usecache, species):


    logging.debug(f"filename={filename}, runname={runname},usecache={usecache},species={species}")
    logging.info("starting...")
    
    logging.info("building ontology...")
    gomatrix = build_ontology(config, usecache)
    logging.info(f"got ontology matrix:\n{matrix_info(gomatrix)} ")
  
    #logging.info("getting uniprot/sprot info..")
    #ubt = get_uniprot_byterm(config, usecache)
    #logging.info(f"got ubt list:\n{ubt[0:20]}")

    #logging.info("getting uniprot/sprot info df..")
    #ubtdf = get_uniprot_byterm_df(config, usecache)
    #logging.info(f"got ubtdf:\n{ubtdf}")

    #logging.info("calculating prior..")
    #priormatrix = calc_prior(config, usecache, species=None)
    #logging.info(f"priormatrix = {matrix_info(priormatrix)}")
    #priordf = get_prior_df(config, usecache)
    #logging.info(f"priordf:\n{priordf}")


    logging.info("running phmmer")
    pdf = get_phmmer_df(config, filename)
    logging.info(f"got phmmer df:\n{pdf}")

    logging.info("making phmmer prediction...")
    out = calc_phmmer_prediction(config, pdf, usecache)
    logging.debug(f"prediction=\n{out}")

    #logging.info("generating test targetset...")
    #testfile = generate_targetset(config, species)
    #logging.info(f"done generating targetset in file {testfile}")
      
       
    #logging.info("evaluate test prediction...")
       
       
    #logging.info("creating species maps...")
    #map = build_specmaps(config, usecache)
    #logging.info("got species map.")    
    

    
    
    #priormatrix = calc_prior(config, usecache, species)    
    ###priormatrix = priormatrix.transpose()
    #logging.debug(f"priormatrix shape: {priormatrix.shape}")
    #gotermlist=list(GOTERMIDX)
    #logging.debug(f"terms in goterm list: {len(gotermlist)}")
    #df = pd.DataFrame(priormatrix, index=gotermlist, columns=['pest'])
    #df = df.sort_values(by='pest', ascending=False)
    
    logging.info("done.")


def run_phmmer(config, filename):
    (outfile, exclude_list) = execute_phmmer(config, filename)
    phdict = parse_phmmer(config, outfile, exclude_list )
    return phdict

def get_phmmer_dict(config, filepath):
    logging.info("running phummer")
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
    
    
    """
    exclude_list = None 
    exclude_list = get_tfa_geneids(filename)
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
    return (outfile, exclude_list)

def get_tfa_geneids(filename):
    """
    Get the geneids of the TargetFiles in order to exclude these results from, e.g. 
    phmmer runs. This wouldn't necessarily cause a problem with the real files, because they 
    will be un-annotated. But for validation testing, the targets *will* be annotated (and
    need to be excluded). 
        
    """
    gids = []
    
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
            gids.append(geneid)
    logging.debug(f"got {len(gids)} geneids to exclude from results.")
    return gids
            

def parse_phmmer(config, filename, excludelist):
    """
    Read phmmer tbl out. Return dict  
    
    excludelist = [ '<pids of inbound target geneid.>'] e.g. 2A5D_ARATH
     
    remove_self_hits = True|False
        shouldn't be necessary in standard case. 
    
    topx_threshold=<int>
        only include top X hits 
    
    """
    logging.info(f"Reading {filename}")
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
           
    topx = config.getint('phmmer','topx_threshold')
    if topx is not None:
        df = df.groupby('cid').head(topx).reset_index(drop=True) 
    
    
    dict = df.to_dict(orient='index')
    
    idxtodel = []
    for idx in dict.keys():
         (db,pacc,pid) = dict[idx]['target'].split('|')
         if pid in excludelist:
             logging.debug(f"Found exluded pid {pid} in exclude list. Removing...")
             idxtodel.append(idx)
         else:
             (protein, species) = pid.split('_')
             dict[idx]['pacc'] = pacc
             dict[idx]['pid'] = pid
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


def calc_phmmer_prediction(config, dataframe, usecache):
    """
    Takes phmmer df PDF: 
                cid           eval  score  bias    pacc          pid
1     T100900000001  4.100000e-155  518.4   7.7  P35213    1433B_RAT
2     T100900000001  5.400000e-155  518.0   7.2  A4K2U9  1433B_PONAB
    
    Gets uniprot_byterm_df UBTDF:
         pacc    speecies      goterm   goev
0        Q6GZX4   FRG3G     GO:0046782   IEA
1        Q6GZX3   FRG3G     GO:0033644   IEA
    
    Gets gomatrix:
type: <class 'numpy.ndarray'> shape: (47417, 47417) dtype: bool     

    Gets gotermidx map:
      47417    

    Algorithm:
    
    CIDGV = govector(0)
    for each cid in PDF:
        GV  = govector(0)
        for each row of PDF (where cid == target):
            for each (goterm) row of UBTDF (where pacc == pacc):
                  UGV = gomatrix[gotermidx[goterm]]
                  GV = GV + UGV
            #  if we want to weight by score, rather than topX...
            # GV * PDF.score (for that pacc)
            CIDGV = CIDGV + GV
    [ for a target we now have a goterm vector with aggregate scores ]
    CIDGV to DF:   cid  goterm  score
    sort by cid, score   
    [ normalize score to .01 - 1.0 ]
    
    create dataframe with 
      cid         goterm       pest
   
    """
    logging.debug("getting uniprot_byterm_df..")
    ubtdf = get_uniprot_byterm_df(config, usecache)
    logging.debug("getting gomatrix...")
    gomatrix = get_ontology_matrix(config, usecache)
    logging.debug("Local goterm index for lookup...")
    gotermidx = GOTERMIDX
    example = gomatrix[gotermidx['GO:0005737']]
    logging.debug(f"example: term GO:0005737: {example}")
    logging.debug(f"example: sum={example.sum()} min={example.sum()} max={example.max()}")
    example = example
    logging.debug(f" ubtdf:\n{ubtdf}\ngomatrix: {matrix_info(gomatrix)}\ngotermidx: {len(gotermidx)} items.")
    
    max_goterms = config.getint('global','max_goterms')
    
    #phmmer dataframe inbound. 
    pdf = dataframe
    cidlist = list(pdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")
    
    gtarray = np.array(list(GOTERMIDX))
    
    topdf = pd.DataFrame(columns=['cid','goterm','pest'])
    
    for cid in cidlist:
        gv = np.zeros(len(gotermidx))
        cdf = pdf[pdf.cid == cid]
        #logging.debug(f"one cid df:\n{cdf}")
        for (i, ser) in cdf.iterrows():
            #logging.debug(f"pacc is {ser.pacc}")
            udf = ubtdf[ubtdf.pacc == ser.pacc]
            for (j, upser) in udf.iterrows():
                #logging.debug(f"goterm is {upser.goterm} gv={gv} addgv={gomatrix[gotermidx[upser.goterm]].astype(np.int8)}")
                gv = gv + gomatrix[gotermidx[upser.goterm]].astype(np.int64)
            #logging.debug(f"sum is {gv.sum()} ")
        # we now have a govector with all goterms indicated by all targets.
        # each entry is sum of number of times that goterm was indicated (as annoated or
        # parent of an annotation term).
        # gv = array([123,   0,   3,   7, 345, 0])
        # handle this cid
        
        # gv has been coerced from array to matrix, back to array??
        #gv = 
        logging.debug(f"gv: {matrix_debug(gv)}")
        gvnz = gv > 0
        logging.debug(f"gvnz: {matrix_debug(gvnz)}")
        logging.debug(f"gtarray:length: {len(gtarray)} type:{type(gtarray)}\n{gtarray}")
        gotermar = gtarray[gvnz]
        #logging.debug(f"goterms with val > 0: {gotermar}")

        govalar = gv[gvnz]
        #logging.debug(f"values for element with val> 0: {govalar}")
        cidar = np.full( len(govalar), fill_value=cid)      
        df = pd.DataFrame({ 'cid': cidar, 'goterm': gotermar, 'pest' : govalar })
        df = df.nlargest(max_goterms, 'pest')
        #df.sort_values(by='pest', ascending=False)
        logging.debug(f"made dataframe for cid {cid}:\n{df}")
        topdf = topdf.append(df, ignore_index=True)
        
    logging.debug(f"made dataframe for all:\n{topdf}")

    return topdf



def generate_targetset(config, name,  species):
    """
    
    Creates random target set from annotated swissprot proteins for given species. 
    
    Creates FASTA files exactly like CAFA TargetFiles
    <name>.<taxonid>.tfa
    
        
    :param    type    name:   desc
    :return    testfile  Path to test file generated in fasta format.
    :rtype
    :raise:   Xerror if blah. 
    
    """


def get_ontology_matrix(config, usecache):
    gmtx = build_ontology(config, usecache)
    mi = matrix_debug(gmtx)
    logging.debug(f"got ontology matrix: {mi} ")
    return gmtx


def build_ontology(config, usecache=False):
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
    global GOMATRIX
    global GOTERMIDX
    global ALTIDDICT
    global GOTERMLIST
    
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('ontology','cachedir'))
    cachefile = f"{cachedir}/ontology.npy"
    termidxfile = f"{cachedir}/gotermindex.pickle"
    altiddictfile = f"{cachedir}/altiddict.pickle"
    include_partof = config.getboolean('ontology','include_partof')
    
    gomatrix = None
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing matrix...")
        gomatrix = np.load(cachefile)
        logging.debug(f"loaded matrix: {matrix_info(gomatrix)} from {cachefile}")
        
        f = open(termidxfile, 'rb')
        termidx = pickle.load(f)
        f.close()
        GOTERMIDX = termidx
        GOTERMLIST = list(GOTERMIDX)
        
        f = open(altiddictfile, 'rb')
        altiddict = pickle.load(f)
        f.close()
        ALTIDDICT = altiddict
        
        logging.debug(f"goterm index, e.g. : \n{list(termidx)[0]} :  {termidx[list(termidx)[0]]} ")
    
    else:
        (godict, altids) = parse_obo(config)
        logging.debug(f"storing {len(altids)} altids to global.")
        ALTIDDICT = altids
        f = open(altiddictfile, 'wb')
        pickle.dump(altids, f)
        f.close()        
        # get keys from dict
        gotermlist = list(godict)
        logging.debug(f"parsed obo with {len(gotermlist)} entries. ")
        logging.debug(f"example entry:\n{gotermlist[0]}")
        logging.debug("sorting goterms")
        gotermlist.sort()
        logging.debug(f"sorted: e.g. {gotermlist[0:5]} ")
        logging.debug("creating goterm index dict.")
        termidx = {}
        i = 0
        for gt in gotermlist:
            termidx[gt] = i
            i = i + 1
        GOTERMIDX = termidx
        GOTERMLIS = list(GOTERMIDX)
        
        f = open(termidxfile, 'wb')   
        pickle.dump(termidx, f )
        f.close()

        logging.debug(f"creating zero matrix of dimension {len(gotermlist)}")
        shape = (len(gotermlist), len(gotermlist))
        gomatrix = np.zeros( shape, dtype=bool )
        logging.debug(f"filling in parent matrix for all goterms...")
        for gt in godict.keys():
            for parent in godict[gt]['is_a']:
                    gomatrix[termidx[parent]][termidx[gt]] = True
        if include_partof:
            logging.debug("Including part_of relationships as is_a")
            for gt in godict.keys():
                for parent in godict[gt]['part_of']:
                        gomatrix[termidx[parent]][termidx[gt]] = True    
             
        logging.debug("Calculating sparsity...")
        logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")
        logging.debug("converting to sparse matrix.")
        gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
        logging.debug("converging matrix...")
        gomatrix = converge_sparse(gomatrix)
        logging.info(f"got converged matrix:\n{matrix_info(gomatrix)} ")
        logging.debug(f"converged matrix sum={gomatrix.sum()}")
        logging.debug("Calculating sparsity...")
        #sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
        #logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")        
        gomatrix = gomatrix.todense()
        # gomatrix = gomatrix.asA()
        gomatrix = np.asarray(gomatrix)
        logging.debug(f"saving matrix: type: {matrix_info(gomatrix)} to {cachefile}")
        np.save(cachefile, gomatrix)

    GOMATRIX = gomatrix
    return gomatrix


def converge_sparse(matrix):
    logging.debug(f"starting matrix: \n{matrix_info(matrix)}")
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



def get_prior_df(config, usecache=False, species = None):
    if GOTERMIDX is None:  
        build_ontology(config, usecache)
    
    priormatrix = calc_prior(config, usecache, species)    
    priormatrix = priormatrix.transpose()
    logging.debug(f"priormatrix shape: {priormatrix.shape}")
    gotermlist=list(GOTERMIDX)
    logging.debug(f"terms in goterm list: {len(gotermlist)}")
    df = pd.DataFrame(priormatrix, index=gotermlist, columns=['pest'])
    df = df.sort_values(by='pest', ascending=False)
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
    global GOMATRIX
    global GOTERMIDX
    global GOPRIOR
 

    logging.debug(f"building sprot. species={species}")
    sprot = get_uniprot_byterm(config, usecache)
    logging.debug(f"sprot, e.g.:\n{pp.pformat(sprot[0:5])} ... ")
   
    freqarray = None
    
    logging.debug(f"usecache={usecache}")
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
        df = pd.DataFrame(sprot, columns=['proteinacc','species', 'goterm', 'goevidence'])
        if species is not None: 
            logging.debug(f"species {species} specified. converted to df with {df.shape[0]} rows: {df}")
            df = df[df.species == species ]
            logging.debug(f"removed other species. {df.shape[0]} rows left.")
        
        sprot = df.to_numpy().tolist()
        logging.debug(f"got uniprot by term, e.g.:\n{pp.pformat(sprot[0:5])} ")
        
        gomatrix = GOMATRIX
        gomatrix = gomatrix.astype(np.int)
        logging.debug(f"gomatrix e.g. {gomatrix[0:5]}")
        gtidx = GOTERMIDX   
        #logging.debug(f"gtidx: {gtidx}")
        shape = (len(gomatrix))
        sumarray = np.zeros(shape, dtype=np.int)
        logging.debug(f"made array {sumarray}")
        i = 0
        e = 0 
        elist = []
        logging.debug("summing over all terms and parents...")
        for r in sprot:
            gt = r[2]
            try:
                sumarray = sumarray + gomatrix[gtidx[gt]]
                i += 1
            except KeyError:
                #logging.debug(f"got missing key: {gt} look up"  )
                primeid = ALTIDDICT[gt]
                #logging.debug(f"primeid is: {primeid}, using...")
                sumarray = sumarray + gomatrix[gtidx[ primeid ]]
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
    GOPRIOR = freqarray
    return freqarray


def get_uniprot_df(config, usecache=True):
    lod = build_uniprot(config, usecache)
    df = pd.DataFrame(lod )
    logging.debug(f"Built dataframe:\n{df}")    
    return df    


def get_uniprot_byterm_df(config, usecache=True):
    lol = get_uniprot_byterm(config, usecache)
    df = pd.DataFrame(lol,columns=['pacc','species','goterm','goev'])
    logging.debug(f"Built dataframe:\n{df}")    
    return df 


def get_uniprot_byterm(config, usecache):
    """
    [ {'proteinid': '3AHD_EGGLE',
      'protein': '3AHD',
      'species': 'EGGLE',
      'proteinacc': 'C8WMP0',
      'taxonid': '479437',
      'goterms': {'GO:0047043': ['mf', 'IEA'],
       'GO:0006694': ['bp', 'IEA']}}),
      ...
    ]
    
    Generate non-redundant uniprot with a row for every goterm:
       pacc   species   goterm       goev
    0  '3AHDP'      'RUMGV'   'GO:0016491'  'IEA'
    1  '3AHDP'      'RUMGV'   'GO:0006694'  'TAS'
    
    
    """    
    lod = build_uniprot(config, usecache)
    ubt = []
    for p in lod:
        for gt in p['goterms'].keys():
            item = [ p['proteinacc'],
                     p['species'],
                     gt, 
                     p['goterms'][gt][1] 
                    ]
            ubt.append(item)
    logging.debug(f"created uniprot_byterm with {len(ubt)} entries.")          
    return ubt          


def build_uniprot(config, usecache):
    """
    
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
            
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()        
    
    #df = pd.DataFrame(listofdicts)
    #logging.debug(f"Built dataframe {df}")
    return listofdicts


def build_specmaps(config, usecache=False):
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
                    current['goterms'][goterm] = [ goaspect, goevidence]

                elif line.startswith("SQ   SEQUENCE"):
                    #logging.debug("Handling SQ:  XXX")
                    # line = filehandle.readline()
                    pass

                elif line.startswith("GN   "):
                    # Examples:
                    #  GN   ABL1 {ECO:0000303|PubMed:21546455},
                    #  GN   Name=BRCA1; Synonyms=RNF53;
                    #  GN   ORFNames=T13E15.24/T13E15.23, T14P1.25/T14P1.24;
                    #   
                    
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
                species = line[:5]
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
    #logging.debug("")
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
    
    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='a .fasta sequence files')
    


    parser.add_argument('-n', '--name',
                        action="store", 
                        dest='runname', 
                        default='default',
                        help='Run-specific identifier to use in file output.')

    parser.add_argument('-s', '--species',
                        action="store", 
                        dest='species', 
                        default=None,
                        help='Limit to species where relevant.')

    parser.add_argument('-a', '--aspect',
                        action="store", 
                        dest='goaspect', 
                        default=None,
                        help='Generate/score only for specific GO aspect. [bp|cc|mf] ')

    parser.add_argument('-C', '--usecache',
                        action='store_true', 
                        dest='usecache',
                        default=False, 
                        help='Use any cached information to speed processing.')
                    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
           
    dorun(cp, args.infile, args.runname, args.usecache, args.species)
        