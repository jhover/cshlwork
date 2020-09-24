#!/usr/bin/env python
import os
import sys
import logging
import traceback

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

from fastcafa.fastcafa import *

# Test files
#dupepairs = os.path.expanduser("~/play/hamsini/mouse_dup_pairs_uniprot_10.txt")
#dupefasta = os.path.expanduser("~/play/hamsini/dupe_pairs_10.tfa")
#targetfasta = os.path.expanduser("~/play/hamsini/dupe_targets_10.tfa")
#phmmerdf = os.path.expanduser("~/play/hamsini/dupe_targets_10_phmmer.csv")


# Real files
dupepairs = os.path.expanduser("~/play/hamsini/mouse_dup_pairs_uniprot.txt")
dupefasta = os.path.expanduser("~/play/hamsini/dupe_pairs.tfa")
targetfasta = os.path.expanduser("~/play/hamsini/dupe_targets.tfa")
phmmerdf = os.path.expanduser("~/play/hamsini/dupe_targets_phmmer.csv")

uniprot_fasta='/Users/jhover/data/uniprot/uniprot_mouse_all.fasta'

def indexbypacc(lod):
    logging.debug("indexing uniprot list of dicts...")
    upbypacc = {}
    for p in lod:
        pacc = p['proteinacc']
        seq = p['sequence']
        upbypacc[pacc] = seq
    return upbypacc

def parse_dupepairs():
    f = open(dupepairs)
    lines = f.readlines()
    dupelist = []
    for line in lines:
        (p1, p2) = line.split(',')
        p1 = p1.strip()
        p2 = p2.strip()
        dupelist.append( (p1, p2) )
    logging.debug(f"dupelist[1] = {dupelist[1]}")
    return dupelist

def write_sequences(dupelist, upbypacc):
    outfile = dupefasta
    snum = 0
    tnum = 0
    mnum = 0
    x = 60
    s = ""
    t = ""
    #for p in spl:
    for tup in dupelist:
        (p1, p2) = tup
        try:
            seq = upbypacc[p1]
            #r = naspec[naspec.pacc == p].reset_index().iloc[0]
            s += f">{p1} {p2}\n"
            chunklist = [ seq[y-x:y] for y in range(x, len(seq)+x, x) ] 
            for c in chunklist:
                s += f"{c}\n"
            snum += 1

            seq2 = upbypacc[p2]
            #r = naspec[naspec.pacc == p].reset_index().iloc[0]
            t += f">{p2}\n"
            chunklist = [ seq2[y-x:y] for y in range(x, len(seq2)+x, x) ] 
            for c in chunklist:
                t += f"{c}\n"
            tnum += 1
        
        except KeyError:
            mnum += 1
            logging.warning(f"Key missing: {p1}")
        
        
    logging.debug(f"handled {snum} dupe sequences. {mnum} missing. ")
    try:
        f = open(outfile, 'w')
        f.write(s)
        logging.debug(f"Wrote TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()    

    outfile = targetfasta
    try:
        f = open(outfile, 'w')
        f.write(t)
        logging.debug(f"Wrote TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()           

def parse_uniprot_fasta(infile):
    '''
    parses fasta
    returns
    
    '''
    
    try:
        f = open(infile)
        lines = f.readlines()
        lod = []
        
        seq = None
        current = None   
        
        for line in lines:
            if line.startswith("#"):
                pass
            
            elif line.startswith(">"):
                if current is not None:
                    current['sequence'] = seq
                    lod.append(current)
                current = {}                   
                seq = ""
                fields = line[1:].split()
                (db, pacc, pid) = fields[0].split('|')
                desc = " ".join(fields[1:])
                current['proteinacc'] = pacc
            else:
                seq += line.strip()
    except IOError:
        logging.error(f"could not read file {infile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()        
    
    logging.debug(f"got lod len: {len(lod)}, e.g. {lod[0]}")
    return lod


if __name__=='__main__':

    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)
    #logging.getLogger().setLevel(logging.DEBUG)
    logging.debug("dupe_similarity...")
    
    config = get_default_config()
    up = parse_uniprot_fasta(uniprot_fasta)   

    logging.debug(f"up len: {len(up)}")
    upbypacc = indexbypacc(up)
    logging.debug(f"upbypacc len: {len(upbypacc)}")

    pairlist = parse_dupepairs()
    write_sequences( pairlist, upbypacc )
    
    outfile, exclude_list, cidgidmap = execute_phmmer(config, dupefasta, version='current')
    logging.info(f"wrote phmmer output to {outfile}")

    df = get_phmmer_df(config, dupefasta)
    logging.debug(f"df: {df}")

    df.to_csv(phmmerdf)
    logging.debug(f"Wrote phmmer DF to {phmmerdf}")
        




