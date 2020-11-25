#!/usr/bin/env python
#
#
#   netid2fasta.py  <species>_[meta|prio]AggNet.hdf5  <species>_info.csv <uniprot-all.fasta> 
#    
#   "EntrezID","EnsemblID","GeneSymbol","LocusTag","Synonyms","dbXrefs",
#        "Chromosome","Type","UniProtID","UniProtID2","OrthoID","NetworkIDs"
##
#
#  >sp|Q8CJG1|AGO1_MOUSE  Protein argonaute-1 OS=Mus musculus OX=10090 GN=Ago1 PE=1 SV=2
#
#

import argparse
import logging

import pandas as pd
import h5py

def parse_tfa_file(infile):
    """
    Reads .tfa file, determines species, target ids, geneids. 
    
    returns dataframe:
    pacc pid sequence
    
    """
    listoflists = []
    try:
        f = open(infile, 'r')
        currentheader = None
        currentseq = None
                
        for line in f:
            #  tr|B5X0T5|B5X0T5_SALSA Glutamate dehydrogenase OS=Salmo salar OX=8030 GN=DHE3 PE=2 SV=1
            if line.startswith(">"):
                # new entry, deal with old:
                if currentseq is not None:
                    logging.debug(f"header: {currentheader}")
                    logging.debug(f"sequence: {currentseq}")
                    currentheader.append(currentseq)
                    listoflists.append(currentheader)
                    currentseq = None
                # start new entry
                currentheader = line[1:].split()
                currentheader = currentheader[0].split('|')[1:]
            else:
                line = line.strip()
                if currentseq is not None:
                    currentseq = f"{currentseq}{line}"
                else:
                    currentseq = line
    except FileNotFoundError:
        logging.error(f"file not readable {infile} ")
        
    logging.info(f"got {len(listoflists)} entries...")
    logging.debug(f"listoflists: {listoflists}") 
    df = pd.DataFrame(listoflists, columns=['pacc','pid','sequence'])
    #df = pd.DataFrame(priormatrix, index=ontobj.gotermlist, columns=['score'])
    #df.set_index('pacc', inplace=True)
    map = pd.Series(df.sequence.values, index=df.pacc).to_dict()
    
    #logging.debug(f"dimension:  {df.shape}")
    return map


def parse_info_file(filename):
    '''
    return dictionary of netids to UniprotIDs for protein-coding gene/proteins in network. 
    '''
    
    df = pd.read_csv(filename)
    logging.debug(f"initial shape={df.shape}")
    df = df[df.Type == "protein-coding"]
    logging.debug(f"protein-coding shape={df.shape}")
    try:
        logging.debug(f"has column UniProtID.y")
        df = df[df["UniProtID.y"].notna()]
        logging.info(f"length after NA filtering {len(df)}")
        map = pd.Series(df['UniProtID.y'].values, index=df.NetworkIDs).to_dict()
    except:
        logging.debug(f"has column UniProtID")
        df = df[df["UniProtID"].notna()]
        logging.info(f"length after NA filtering {len(df)}")
        map = pd.Series(df['UniProtID'].values, index=df.NetworkIDs).to_dict()

    return map


def parse_expression_hd5(filename):
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
    return df 

def make_uidlist(netidlist, uidmap):
    outlist = []
    nummissing = 0
    numfound = 0
    for id in netidlist:
        try:
            out = uidmap[id]
            tup = (id, out)
            outlist.append(tup)
            numfound += 1
        except KeyError:
            nummissing += 1
            logging.debug(f"missing value for networkid {id}")
    outlist.sort()
    logging.info(f"Generated UniprotID list. {numfound} found. {nummissing} missing")
    return outlist

def write_tfa_file(outlist, seqmap, outfile):
    '''
    outlist is tuple: ( netid, uid )
    '''
    
    s=""
    x = 60
    nummissing = 0
    numfound = 0
    
    for (netid, uid ) in outlist:    
        try:
            sequence = seqmap[uid]
            numfound += 1
            s += f">{netid} {uid}\n"
            chunklist = [ sequence[y-x:y] for y in range(x, len(sequence)+x, x) ] 
            #logging.debug(f"chunklist {chunklist}")
            for c in chunklist:
                s += f"{c}\n"
        except KeyError:
            logging.warning(f"no sequence found for uid {uid} !")
            nummissing += 1
    
    if nummissing > 10:
        logging.warning(f"{nummissing} sequences missing from DB!")
    logging.info(f"found {numfound} sequences. {nummissing} missing")
    try:
        f = open(outfile, 'w')
        f.write(s)
        logging.info(f"Wrote TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()


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

    parser.add_argument('netfile', 
                        metavar='netfile', 
                        type=str, 
                        help='<species>_[prio|meta]AggNet.hdf5 file')
    
    parser.add_argument('infofile', 
                        metavar='infofile', 
                        type=str, 
                        help='<species>_info.csv file')

    parser.add_argument('fastafile', 
                        metavar='fastafile', 
                        type=str, 
                        help='uniprot fasta file')

    parser.add_argument('outfile', 
                        metavar='outfile', 
                        type=str, 
                        help='output fasta file')

   
    args= parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    logging.debug(f"netfile={args.netfile} infofile={args.infofile}  fastafile={args.fastafile} ")
    
    edf = parse_expression_hd5(args.netfile)
    logging.debug(f"{edf.head()}")

    mapdict = parse_info_file(args.infofile)    
    logging.debug(f"got mapdict={mapdict}")

    outlist = make_uidlist(list(edf.columns), mapdict)
    logging.debug(f"got outlist = {outlist}")
    
    smap = parse_tfa_file(args.fastafile)
    logging.debug(f"got smap={smap}")
    
    write_tfa_file(outlist, smap, args.outfile)
        