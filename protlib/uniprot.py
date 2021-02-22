#!/usr/bin/env python
# 
# Handles interaction with uniprot data.
#
import argparse
import logging
import os
import sys
import traceback

from collections import defaultdict

ASPECTMAP = { 'C': 'cc',
              'F': 'mf',
              'P': 'bp'
            }



#def parse_uniprot_dat(config, version='2019'):
def parse_uniprot_dat(filepath):
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
                    #logging.debug("Handling ID. New entry.")                
                
                elif line.startswith("AC "):
                    # AC   Q6GZX4;
                    # AC   A0A023GPJ0; 
                    # AC   Q91896; O57469;
                    #logging.debug("Handling AC.")
                    rest = line[5:]
                    accession = rest.split(';')[0]
                    #accession = line[5:11].strip()
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
                    goaspect = ASPECTMAP[aspcode]
                    goevsrc = fields[3]
                    (goevidence, evsrc) = goevsrc.split(':') 
                    goevidence = goevidence.strip()
                    current['goterms'][goterm] = [ goaspect,  goevidence ]

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
                    # logging.debug(f"Handling GN. {line}")
                    val = line[5:]
                    if val.startswith("Name="):
                        fields = val.split()   # by whitespace
                        (n, gname) = fields[0].split("=")
                        gname = gname.upper()
                        #logging.debug(f"Gene name is {gname} ")
                        current['gene'] = gname.replace(';','') 
            
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
        logging.debug(f"Some entries:  {allentries[10:15]}")
        return allentries

def write_tfa_file(dlist, outfile):
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
    
    #for index, row in df.iterrows():
    for p in dlist:
        # to be compatible (db, pacc, pid) = fields[0].split('|')
        header = f">up|{p['proteinacc']}|{p['protein']}_{p['species']} {p['gene']} "
        sequence =  p['sequence']
        s += f"{header}\n"
        chunklist = [ sequence[y-x:y] for y in range(x, len(sequence)+x, x) ] 
        #logging.debug(f"chunklist {chunklist}")
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
    
    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str,
                        nargs='?',
                        default=os.path.expanduser("~/data/uniprot/uniprot_all_vertebrates.dat"), 
                        help='Uniprot .dat file') 
    
    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    
    logging.info(f"Parsing {args.infile} ...")
    uplist = parse_uniprot_dat(args.infile)
    print(f"outlist is length = {len(uplist)}")    
