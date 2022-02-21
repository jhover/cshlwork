#!/usr/bin/env python
# 
# Handles interaction with uniprot data.
#
#  uniprot vertebrates. 6.7 million entries. 2 minutes from cache. vs. 14minutes.  
#   
#
#   Uniprot entries  have 1 ID, but sometimes multiple AC s. 
#   Parsing should index by both. 
#
#
#  Test CLI:
#    read/parse/cache uniprot.dat
#
#    .dat as arg -> use, otherwise .dat in config
#
#    list of proteins + tfa:   write protein sequences to tfa file. 
#    protein(s) as arg     :   write protein sequences to stdout 
#    tfa file only         :   write all protein sequences in .dat to tfa
#    nothing               :   write all to stdout

import argparse
import logging
import os
import pickle
import sys
import traceback

from collections import defaultdict
from configparser import ConfigParser

import pandas as pd

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from protlib.utils import *

ASPECTMAP = { 'C': 'cc',
              'F': 'mf',
              'P': 'bp'
            }

def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cshlwork/etc/uniprot.conf"))
    return cp

def parse_uniprot_dat(config, datfile=None):
        """
        Parses uniprot/sprot DAT file, returns a list of dicts. 
        
        [   raw-entries ,    # all the entries as raw dicts. 
            pididx,          # dictionary indexed by pid -> entry
            accidx          # dictionary indexed by accession -> entry
        ]
        
        # raw entries            
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
        }

        """
        uplist = ()      #   tuple: entrylist, pididx, accidx  
        entrylist = []   #     list of dicts
        pididx = {}        
        accidx = {}
       
        cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
        filepath = os.path.expanduser(config.get('uniprot','datfile'))
        if datfile is not None:
            filepath = datfile
        filebase = os.path.basename(filepath)
        (filebase, e) = os.path.splitext(filebase)
        cachefile =f"{cachedir}/{filebase}.allbypacc.pickle"
        
        if os.path.exists(cachefile):
            logging.info("Cache hit. Using existing data...")
            f = open(cachefile, 'rb')
            uniprot_data = pickle.load(f)
            logging.info("Loaded from cache...")
            f.close()      
        
        else:
            logging.debug(f"opening datfile={filepath}")
            try:
                logging.debug(f" attempting to open '{filepath}'")
                filehandle = open(filepath, 'r')
            except FileNotFoundError:
                logging.error(f"No such file {filepath}")
                raise                
            
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
                        val = line[5:]
                        fields = val.split() 
                        proteinid = fields[0].strip().replace(';','')
                        current = defaultdict(dict)
                        current['proteinid'] = proteinid
                        (protein, species) = proteinid.split('_')
                        current['protein'] = protein
                        current['species'] = species
                        pididx[proteinid] = current
            
                    elif line.startswith("AC "):
                        # AC   Q6GZX4;
                        # AC   A0A023GPJ0; 
                        # AC   Q91896; O57469;
                        rest = line[5:]
                        acclist = rest.split(';')
                        current['proteinacc'] = acclist[0].strip()
                        for c in acclist:
                            if len(c) > 2:
                                c = c.strip()
                                accidx[c] = current
                                                
                    elif line.startswith("OX   "):
                        #OX   NCBI_TaxID=654924;
                        taxonid = ""
                        val = line[5:]
                        fields = val.split('=')
                        if fields[0] == 'NCBI_TaxID':
                            tfields = fields[1].split()
                            taxonid = tfields[0].strip().replace(';','')
                            #taxonid = fields[1].strip().replace(';','')
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
                        # standard example
                        #  GN   Name=GRF1;
                        
                        #  w/ extra info, but all on one line per protein
                        #  GN   Name=BRCA1; Synonyms=RNF53;
                        #  GN   Name=GRF10; OrderedLocusNames=At1g22300; ORFNames=T16E15.8;
                        #  GN   Name=GRF2; Synonyms=GF14; OrderedLocusNames=At1g78300; ORFNames=F3F9.16;

                        #  non Name= info only
                        #  GN   ORFNames=OsI_006296;
                                               
                                                                      
                        # multiple lines in one protein (with Name=)
                        # GN   Name=GRF12 {ECO:0000303|PubMed:11553742};
                        # GN   OrderedLocusNames=At1g26480 {ECO:0000312|Araport:AT1G26480};
                        # GN   ORFNames=T1K7.15 {ECO:0000312|EMBL:AAF98570.1};
                        
                        # multi-line, no key(s)
                        # GN   Name=matK {ECO:0000256|HAMAP-Rule:MF_01390,
                        # GN   ECO:0000313|EMBL:ACK76147.1};
                        
                        # multiple lines in one protein (no Name=)
                        #  GN   OrderedLocusNames=Os02g0224200, LOC_Os02g13110;
                        #  GN   ORFNames=OsJ_005772, P0470A03.14;
                        val = line[5:]
                        if val.startswith("Name="):
                            try:
                                fields = val.split()   # by whitespace
                                (n, gname) = fields[0].split("=")
                                gname = gname.upper()
                                gname = gname.replace(';','')
                                gname = gname.replace(':','')
                                current['gene'] = gname
                            except ValueError as ve:
                                logging.error(f"val= {val}")
                                traceback.print_exc(file=sys.stdout)
                                current['gene'] = '' 
                                 
                        #elif val.startswith("")
                        else:
                            current['gene'] = '' 
                
                    elif line.startswith("//"):          
                        entrylist.append(current)
                                             
                        if len(entrylist) >= repthresh:
                            logging.info(f"Processed {len(entrylist)} entries... ")
                            sumreport +=1
                            repthresh = sumreport * suminterval
                        current = None
                    
            except Exception as e:
                traceback.print_exc(file=sys.stdout)                
            
            if filehandle is not None:
                filehandle.close()      

                        
            uniprot_data = (entrylist, pididx, accidx)
            
            logging.info(f"Parsed file with {len(entrylist)} entries" )
            logging.info(f"Created prim/sec index. {len(accidx)} entries ")
            logging.debug(f"Some entries:  {entrylist[10:15]}")
            logging.info(f"Caching uniprot list to {cachefile}...")
            f = open(cachefile, 'wb')   
            pickle.dump(uniprot_data, f )
            f.close()
            logging.debug("Done.")
        return uniprot_data


def uniprot_to_df(cp):
    """
    
    """
    COLUMNS = ['pacc', 'proteinid', 'protein', 'species', 'proteinacc', 'gene', 'taxonid']
    XCOLS = ['proteinid', 'protein', 'species', 'proteinacc', 'gene', 'taxonid']
    pidx = parse_uniprot_dat(cp) 
    (entries, pididx, accidx) = parse_uniprot_dat(cp)
    logging.debug('Done parsing uniprot .dat file. Building LOL.')
    lol = []
    for k in pididx.keys():
        e = pidix[k]
        flist = [k]
        for col in XCOLS:
            v = e[col]
            if len(v) == 0:
                v = ''
            flist.append(v)
        lol.append(flist)
    logging.debug(f'made lol with {len(lol)} entries. Making DF...')         
    df = pd.DataFrame(lol, columns=COLUMNS)
    logging.debug(f'completed DF.')
    return df


def write_df_tsv(dataframe, filepath):
    """
    
    
    """
    dataframe.to_csv(filepath, sep='\t', index=False)
    


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
      

def get_fasta(pdict, keylist=None):
    """
    generate fasta entry string from uniprot dict indexed by pid or other. 
    
    """
    s=""
    snum = 1
    header=""
    x = 60
    if keylist is None:
        keylist = pdict.keys()
        
    for k in keylist:
        logging.debug(f'retrieving {k}')
        try:
            p = pdict[k]
            logging.info(f'p = {p}')
            header = f">{p['proteinacc']}\t{p['protein']}\t{p['species']}\t{p['gene']}"
            header = header.replace('{}','')  # remove missing values. 
            sequence =  p['sequence']
            s += f"{header}\n"
            chunklist = [ sequence[y-x:y] for y in range(x, len(sequence)+x, x) ] 
            for c in chunklist:
                s += f"{c}\n"
                snum += 1
            logging.debug(f'\n{s}')
        except KeyError:
            logging.error(f'no entry found for {k}')
            traceback.print_exc(file=sys.stdout) 
       
    return s

def write_tfa(s, outfile=None):
       
    try:
        if outfile is not None:
            logging.debug(f'opening {outfile}...')
            f = open(outfile, 'w')
        else:
            logging.debug(f'outfile is None. Writing to stdout')
            f = sys.stdout
        f.write(s)
        logging.debug(f"Wrote TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()    


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
    
    parser.add_argument('-u', '--uniprot', 
                        metavar='infile', 
                        type=str,
                        nargs='?',
                        default=os.path.expanduser("~/data/uniprot/uniprot_all_vertebrates.dat"), 
                        help='A Uniprot .dat file') 

    parser.add_argument('-t','--tfafile', 
                        metavar='tfafile', 
                        type=str,
                        required=False,
                        default=None, 
                        help='Fasta .tfa output file') 

    parser.add_argument('-p','--pfile', 
                        metavar='pfile', 
                        type=str,
                        required=False,
                        default=None, 
                        help='Protein list file. ')

    parser.add_argument('-i','--identifier', 
                        metavar='identifier', 
                        type=str,
                        required=False,
                        default='accession', 
                        help='Identifier type accession|protein  (B4FCB1 | B4FCB1_MAIZE ')    
    
    
    parser.add_argument('protlist' ,
                        metavar='protlist', 
                        type=str,
                        nargs='*',
                        default=None, 
                        help='UPID [UPID UPID ...] ')

    
    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    
    config = get_default_config()
    logging.info(f"Parsing uniprot .dat file={args.uniprot} ...")
    (entries, pididx, accidx) = parse_uniprot_dat(config, args.uniprot)
    logging.debug(f"uniprot length = {len(entries)}")
    
    tfastr = ""
    
    if args.identifier == 'accession':
        upidx = accidx
    elif args.identifier == 'protein':
        upidx = pididx
    else:
        logging.error(f'invalid identifier type: {args.identifier}')
    
    
    if len(args.protlist) > 0:
        logging.info(f'protlist={args.protlist}')
        tfastr = get_fasta(upidx, args.protlist)
        logging.debug(f'writing to {args.tfafile}')
        write_tfa(tfastr, args.tfafile)
        
    elif args.pfile is not None:
        logging.info(f'fpile={args.pfile}')
        idlist = read_identifier_file(args.pfile)
        logging.debug(f'prot list={idlist}')
        tfastr = get_fasta(upidx, idlist)
        logging.debug(f'writing to {args.tfafile}')
        write_tfa(tfastr, args.tfafile)                
    else:
        logging.debug("no proteins our outfile specified. write all to stdout. ")
        tfastr = get_fasta(upidx)
        write_tfa(tfastr, args.tfafile)                        

    logging.info("done.")


