#!/usr/bin/env python
# 
# Handles interaction with uniprot data.
#
#  uniprot vertebrates. 6.7 millino entries. 2 minutes from cache. vs. 14minutes.  
#   
import argparse
import logging
import os
import pickle
import sys
import traceback

from collections import defaultdict
from configparser import ConfigParser

import pandas as pd


ASPECTMAP = { 'C': 'cc',
              'F': 'mf',
              'P': 'bp'
            }

def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cshlwork/etc/uniprot.conf"))
    return cp

def parse_uniprot_dat(config, filepath):
        """
        Parses uniprot/sprot DAT file, returns dictionary of dicts 
        using primary and secondary accession codes as keys.  
                
        { 'Q84P23' :
            {    'proteinid': '4CLL9_ARATH', 
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
        allentries = None
        paccidx = {}
       
        cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
        filename = os.path.basename(os.path.expanduser(filepath))
        (filebase, e) = os.path.splitext(filename)
        cachefile =f"{cachedir}/{filebase}.allbypacc.pickle"
        
        if os.path.exists(cachefile):
            logging.info("Cache hit. Using existing data...")
            f = open(cachefile, 'rb')
            paccidx = pickle.load(f)
            logging.info("Loaded from cache...")
            f.close()      
        else:
            logging.debug(f"opening datfile={filepath}")
            try:
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
            
                    
                    elif line.startswith("AC "):
                        # AC   Q6GZX4;
                        # AC   A0A023GPJ0; 
                        # AC   Q91896; O57469;
                        rest = line[5:]
                        acclist = rest.split(';')
                        current['proteinacc'] = acclist[0]
                        for c in acclist:
                            if len(c) > 2:
                                paccidx[c] = current
                                                
                    elif line.startswith("OX   "):
                        #OX   NCBI_TaxID=654924;
                        taxonid = ""
                        val = line[5:]
                        fields = val.split('=')
                        if fields[0] == 'NCBI_TaxID':
                            taxonid = fields[1].strip().replace(';','')
                        current['taxonid'] = taxonid
                        
                    elif line.startswith("DR   GO;"):
                        # DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
                        # P biological process, C cellular component, F molecular function.  
                        fields = line.split(';')
                        goterm = fields[1].strip()
                        goinfo = fields[2]
                        aspcode = goinfo.split(':')[0].strip()
                        goaspect = ASPECTMAP[aspcode]
                        goevsrc = fields[3]
                        (goevidence, evsrc) = goevsrc.split(':') 
                        goevidence = goevidence.strip()
                        current['goterms'][goterm] = [ goaspect,  goevidence ]
    
                    elif line.startswith("SQ   SEQUENCE"):
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
                        val = line[5:]
                        if val.startswith("Name="):
                            fields = val.split()   # by whitespace
                            (n, gname) = fields[0].split("=")
                            gname = gname.upper()
                            current['gene'] = gname.replace(';','') 
                
                    elif line.startswith("//"):          
                        allentries.append(current)
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
            logging.info(f"Created prim/sec index. {len(paccidx)} entries ")
            logging.debug(f"Some entries:  {allentries[10:15]}")
            logging.info(f"Caching paccindex to {cachefile}...")
            f = open(cachefile, 'wb')   
            pickle.dump(paccidx, f )
            f.close()
            logging.debug("Done.")
        return paccidx


def index_by_acc(dictlist):
    accidx = {}
    n = 0
    for p in dictlist:
        acc = p['proteinacc']
        accidx[acc] = p
        n += 1
    logging.debug(f"indexed by acc {n} entries...")
    return accidx



def write_tfa_fromlist(pacclist, paccidx, outfile):
    '''
    writes only select protein sequences to tfa outfile. 
    fixes proteinacc inside in case the list entry is an alternate code. 
    
    '''
    newdlist = []
    total = 0
    found = 0
    missing = 0
    missinglist =[]
    for pacc in pacclist:
        total += 1
        try:
            pdict = paccidx[pacc]
            pdict['proteinacc'] = pacc
            newdlist.append(pdict)
            found += 1
        except KeyError:
            missing += 1
            missinglist.append(pacc)
    logging.info(f"total={total} found={found} missing={missing}")            
    logging.info(f"missing list: {missinglist}")
    logging.debug(f"Made shorter dictlist length={len(newdlist)} writing to file={outfile}.")
    write_tfa_file(newdlist, outfile)
      

def write_tfa_file(dictlist, outfile):
    '''
    defaultdict(<class 'dict'>, 
        {'proteinid': '1433Z_XENTR', 
         'protein': '1433Z', 
         'species': 'XENTR', 
         'proteinacc': 'Q6P4Z5', 
         'gene': 'YWHAZ', 
         'taxonid': '8364', 
         'goterms': {'GO:0005737': ['cc', 'IEA'], 
                     'GO:0070372': ['bp', 'ISS']}, 
         'seqlength': 245, 
         'sequence': 'MDKNELVQKAKL...  ...EGGEN'})
    
    '''
    s=""
    snum = 1
    header=""
    x = 60
    for p in dictlist:
        header = f">{p['proteinacc']}\t{p['protein']}\t{p['species']}\t{p['gene']}"
        header = header.replace('{}','')  # remove missing values. 
        sequence =  p['sequence']
        s += f"{header}\n"
        chunklist = [ sequence[y-x:y] for y in range(x, len(sequence)+x, x) ] 
        for c in chunklist:
            s += f"{c}\n"
        snum += 1
    logging.debug(s)   
    
    try:
        f = open(outfile, 'w')
        f.write(s)
        logging.debug(f"Wrote TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()




        

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
                        nargs='?',
                        default=os.path.expanduser("~/data/uniprot/uniprot_all_vertebrates.dat"), 
                        help='Uniprot .dat file') 

    parser.add_argument('-t','--tfafile', 
                        metavar='tfafile', 
                        type=str,
                        required=False,
                        default=None, 
                        help='Fasta .tfa output file') 
    
    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    
    config = get_default_config()
    logging.info(f"Requesting info from {args.infile} ...")
    uplist = parse_uniprot_dat(config, args.infile)
    print(f"outlist is length = {len(uplist)}")
    if args.tfafile is not None:
        write_tfa_file(uplist, os.path.expanduser(args.tfafile))    
    
    logging.info("Generating acc index..")
    idx = index_by_acc(uplist)
    logging.info("done.")


