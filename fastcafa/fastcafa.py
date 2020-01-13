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
from collections import defaultdict
from configparser import ConfigParser
import logging
import pickle
import subprocess
import tempfile
import traceback


import pandas as pd
import numpy as np
np.set_printoptions(threshold=1000)
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


def dorun(config, filename, runname, usecache):
    logging.debug(f"filename={filename}, runname={runname},usecache={usecache}")
    logging.info("starting...")
    #logging.info("running phummer")
    #outfile = runphmmer(config, filename)
    #phdict = parsephmmer(config, outfile )
    #logging.debug(phdict)
    logging.info("building ontology...")
    gomatrix = build_ontology(config, usecache)
    logging.info(f"got ontology matrix:\n{matrix_info(gomatrix)} ")

    
    logging.info("creating species maps...")
    map = build_specmaps(config, usecache)
    logging.info("got species map.")    
    
    logging.info("calculating prior..")
    priormatrix = calc_prior(config, species=None)

    
    logging.info("done.")

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


def build_ontology(config, usecache=False):
    """
    obofile=~/data/go/go.obo
    
    { 'GO:2001315': 
         {'is_a': ['GO:0009226', 'GO:0046349', 'GO:2001313'], 
         'part_of': [], 
         'goterm': 'GO:2001315', 
         'goname': 'UDP-4-deoxy-4-formamido-beta-L-arabinopyranose biosynthetic process', 
         'goasp': 'bp', 
       ...  
    }
    cachedir = ~/play/cafa4   
     
    """
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('ontology','cachedir'))
    cachefile = f"{cachedir}/ontology.npy"
    gomatrix = None
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing matrix...")
        gomatrix = np.load(cachefile)
    else:
        godict = parse_obo(config)
        gotermlist = list(godict)
        logging.debug(f"parsed obo with {len(gotermlist)} entries. ")
        logging.debug("sorting goterms")
        gotermlist.sort()
        logging.debug("creating goterm index dict.")
        termidx = {}
        i = 0
        for gt in gotermlist:
            termidx[gt] = i
            i = i + 1   
        logging.debug(f"creating zero matrix of dimension {len(gotermlist)}")
        shape = (len(gotermlist), len(gotermlist))
        #gomatrix = np.zeros( shape, dtype=np.int8 )
        gomatrix = np.zeros( shape, dtype=bool )
        logging.debug(f"filling in parent matrix for all goterms...")
        for gt in godict.keys():
            for parent in godict[gt]['is_a']:
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
        logging.debug(f"saving matrix: {matrix_info(gomatrix)} to {cachefile}")
        np.save(cachefile, gomatrix)
    return gomatrix


def converge_sparse(matrix):
    logging.debug(f"starting matrix: \n{matrix_info(matrix)}")
    sh = matrix.shape
    logging.debug("Creating ones matrix as mask...")
    ones = np.ones(sh, dtype=matrix.dtype)
    #logging.debug("Converting mask to lil_matrix")
    #ones = sparse.lil_matrix(ones)
    logging.debug("Summing inbound matrix...")
    oldval = matrix.sum()
    newval = 0
    logging.debug("Beginning convergence loop.")
    while oldval != newval:
        logging.debug(f"Inbound matrix:\n{matrix_info(matrix)}")
        logging.debug(f"oldval={oldval} newval={newval}")
        oldval = newval
        logging.debug("If not lil_matrix, convert...")
        if not isinstance(matrix,  sparse.lil.lil_matrix): 
            logging.debug(f"{type(matrix)} is not scipy.sparse.lil.lil_matrix, converting... ")
            matrix = sparse.lil_matrix(matrix, dtype=np.bool)
        else:
            logging.debug("matrix already lil_matrix...")
        logging.debug("Multiplying...")
        mult = matrix @ matrix
        logging.debug("Adding back original...")
        matrix = mult + matrix
        #logging.debug("Taking minimum with mask...")
        #matrix = matrix.minimum(ones)
        logging.debug("Getting new sum...")
        newval = matrix.sum()
        logging.debug(f"matrix:\n{matrix_info(matrix)}")
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
    return f"type: {type(matrix)} shape: {matrix.shape} dtype: {matrix.dtype} sum: {matrix.sum()}\n{matrix}"
    
    
def calc_prior(config,species=None):
    """
    take ontology matrix. goterms x goterms
    build protein x goterm matrix. 
    for each
    
    if species is specified, calculate prior for species only, otherwise global
    
    """
    logging.debug("building sprot...")
    sprot = build_uniprot_df(config)
    logging.debug(f"{sprot}")
    if species is None:
        pass
    else:
        pass
        
    

def build_uniprot_df(config):
    """
    
    """
    listofdicts = parse_uniprot_dat(config)
    df = pd.DataFrame(listofdicts)
    logging.debug(f"Built dataframe {df}")
    return df


def build_specmaps(config, usecache=False):
    """
    builds three maps in form of list of dicts:
    
       taxid   -> speccode, linnean
       speccode  -> taxonid, linnean 
       linnean   -> taxonid, speccode
    
    caches as pickle object. 
    
    """    
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
        
    logging.debug(f"map: {map}")
    logging.info(f"saving map: to {cachefile}")
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
        try:
            filepath = os.path.expanduser(config.get('uniprot','datfile'))
            filehandle = open(filepath, 'r')
        except FileNotFoundError:
            logging.error("No such file %s" % filepath)                
        
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
                
                current['goasp'] = GOASPECTMAP[asp]
            
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
    
    logging.info(f"Parsed file with {len(goidx)} terms")    
    return goidx


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
           
    dorun(cp, args.infile, args.runname, args.usecache)
    
    