#!/usr/bin/env python
#
# run phmmer against comma separated list of Uniprot IDs. 
# produce csv of pairwise match alignment. 
# 
# 
#
#
#
import os
import sys
import logging
import traceback

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)
from fastcafa.fastcafa import *

#dupepairs = os.path.expanduser("~/project/hamsini2/mouse_sgd_genes.txt")
#dupefasta = os.path.expanduser("~/project/hamsini2/mouse_sgd_genes.tfa")
#targetfasta = os.path.expanduser("~/project/hamsini2/dupe_targets.tfa")
#phmmerdf = os.path.expanduser("~/project/hamsini2/dupe_targets_phmmer_20210207.csv")
#uniprot_fasta=os.path.expanduser('~/data/uniprot/uniprot_mouse_all.fasta')
#uniprot_altcodes = os.path.expanduser("~/project/hamsini2/uniprot_all_rodents_altcodes.txt")



def indexbypacc(lod):
    logging.debug(f"indexing uniprot list of dicts len: {len(lod)}")
    upbypacc = {}
    for p in lod:
        pacc = p['proteinacc']
        #if pacc == "A0A0J9YTW6":
        #    logging.debug("Indexing later missing pacc! A0A0R4J0X7")
        seq = p['sequence']
        upbypacc[pacc] = p
    logging.debug(f"produced indexed dict len: {len(upbypacc)}")   
    return upbypacc


def parse_dupepairs(filename):
    f = open(filename)
    lines = f.readlines()
    dupelist = []
    lnum = 0
    knum = 0
    for line in lines:
        (p1, p2) = line.split(',')
        p1 = p1.strip()
        p2 = p2.strip()
        if p2 != "NA":
            dupelist.append( (p1, p2) )
        else:
            knum += 1
            #logging.debug("skipping NA target. ")
            
        lnum += 1
    logging.debug(f" processed {lnum} lines. skipped {knum} NAs. produced {len(dupelist)} items in dupelist[1] = {dupelist[1]}")
    #logging.debug(f"dupelist: {dupelist}")
    return dupelist

def write_sequences(dupelist, upbypacc, pairtfa, targettfa):
    
    #tp = upbypacc['A0A0J9YTW6']
    #logging.debug(f"A0A0J9YTW6 still in upbypacc in write_sequences. ")

    outfile = pairtfa
    qnum = 0
    tnum = 0
    qmnum = 0
    tmnum = 0
    x = 60
    s = ""
    t = ""
    #for p in spl:
    for tup in dupelist:
        (p1, p2) = tup
        try:
            seq = upbypacc[p1]['sequence']
            #r = naspec[naspec.pacc == p].reset_index().iloc[0]
            s += f">{p1} {p2}\n"
            chunklist = [ seq[y-x:y] for y in range(x, len(seq)+x, x) ] 
            for c in chunklist:
                s += f"{c}\n"
            qnum += 1

        except KeyError:
            tmnum += 1
            logging.warning(f"Query key missing: {p1}")
        
    logging.debug(f"handled {qnum} dupe queries. {qmnum} missing. ")
    
    try:
        f = open(outfile, 'w')
        f.write(s)
        logging.debug(f"Wrote query TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()    

    for tup in dupelist:
        (p1, p2) = tup
        p1 =p1.strip()
        p2 = p2.strip()
        #if p2 == "A0A0J9YTW6":
        #    logging.debug(f"found p2 in duplest: A0A0J9YTW6 in p2. checking. ")
        #    pseq = upbypacc[p2]
        #    logging.debug("found p OK!! ")
        
        try:
            seq2 = upbypacc[p2]['sequence']
            #r = naspec[naspec.pacc == p].reset_index().iloc[0]
            t += f">{p2}\n"
            chunklist = [ seq2[y-x:y] for y in range(x, len(seq2)+x, x) ] 
            for c in chunklist:
                t += f"{c}\n"
            tnum += 1
        
        except KeyError:
            tmnum += 1
            logging.warning(f"Target key missing: {p2}")
        
    logging.debug(f"handled {tnum} dupe queries. {tmnum} missing. ")    
    
    outfile = targettfa
    try:
        f = open(outfile, 'w')
        f.write(t)
        logging.debug(f"Wrote target TFA sequences to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()           

def parse_uniprot_fasta(infile):
    '''
    parses fasta
    returns list of dicts:
    
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
                if pacc == "A0A0J9YTW6":
                    logging.debug("Found later missing pacc! A0A0J9YTW6")
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

def add_altcodes(upbypacc, infile):
    '''
      upbypacc   { <pacc> : { 'proteinacc' : <pacc>,
                              'sequence' : <seq> }
                              ,
                              ,
                              ,
                              }
        altcodes:
        
        cat <uniprot>.dat | grep "^AC" > <altcodes>.txt
        
        AC   Q9CQV8; O70455; Q3TY33; Q3UAN6;
        AC   P35213;
        AC   P62259; P29360; P42655; Q63631;

    '''
    logging.debug(f"len upbypacc before: {len(upbypacc)}")
    nadded = 0
    nmissing = 0
    try:
        f = open(infile)
        lines = f.readlines()    
        for line in lines:
            # remove leading AC
            fields = line.split()[1:]
            #logging.debug(f"fields: {fields}")
            if len(fields) > 1:
                #logging.debug("more than one field.")
                ecode = fields[0].replace(';','')
                try:
                    entry = upbypacc[ecode]
                    for alt in fields[1:]:
                        alt = alt.replace(';','')
                        upbypacc[alt] = entry
                        #logging.debug(f"added alt {alt} for entry code {ecode}")
                        nadded += 1
                except KeyError:
                    #logging.warn(f"entry {ecode} not found in upbypacc.")
                    nmissing += 1

    except IOError:
        logging.error(f"could not read file {infile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()    

    logging.debug(f"len ubypacc after: {len(upbypacc)} {nadded} alts added. {nmissing} missing.")


def parse_filebase(filepath):
    '''
    gives back filepath minus the last dot extension, or the
    same filepath if there is not extension. 
    '''
    return os.path.splitext(filepath)[0]


def run_phmmer(pairlist, uniprot_fasta, uniprot_altcodes, pairtfa, targettfa):
    config = get_default_config()
    up = parse_uniprot_fasta(uniprot_fasta)   
    logging.debug(f"up len: {len(up)}")
    upbypacc = indexbypacc(up)   
    add_altcodes(upbypacc, uniprot_altcodes)
    
    logging.debug(f"upbypacc len: {len(upbypacc)}")
    
    write_sequences( pairlist, upbypacc, pairtfa, targettfa  )   
    outfile, exclude_list, cidgidmap = execute_phmmer(config, pairtfa, version='current')
    logging.info(f"wrote phmmer output to {outfile}")
    
    df = get_phmmer_df(config, pairtfa)
    logging.debug(f"df: {df}")
    return df


def get_match(query, target, df):
    logging.debug(f"query={query} target={target}")
    qdf = df[df.cid == query]
    row = qdf[qdf.pacc == target]
    if len(row) > 1 :
        logging.warning(f'multiple matches for query={query} target={target} ')
        return None
    elif len(row) == 1:
        r = row.iloc[0]
        eval = r['eval']
        pscore =r['pscore']
        bias = r['bias']
        return (eval, pscore, bias)
    else:
        logging.warning(f'no matches for query={query} target={target} ')    
        return None

def make_evaltable(pdf, pairlist, evalfile   ):
    #config = get_default_config()
    #pdf  = pd.read_csv(phmmerdf, index_col=0)
    pdf.drop_duplicates(inplace=True,ignore_index=True)   
    #dupelist = parse_dupepairs()

    lod = []
        
    for tup in pairlist:
        (p1, p2) = tup
        logging.debug(f"looking for {p1} -> {p2}")
        rv = get_match(p1, p2, pdf)
        if rv is not None:
            (eval, pscore, bias ) = rv
            lod.append( { 'query' : p1,
                          'target' : p2,
                          'eval' : eval,
                          'pscore' : pscore,
                          'bias' : bias,  
                            }
            )
    logging.debug(f"dupelist length: {len(pairlist)}")
    logging.debug(f"matchlist length: {len(lod)}")
    edf = pd.DataFrame(lod)
    edf.to_csv(evalfile)
    logging.debug(f"wrote match df to {evalfile}")    







if __name__=='__main__':

    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARNING)

    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('pairfile', 
                        metavar='pairfile', 
                        type=str, 
                        help='')

    parser.add_argument('uniprottfa', 
                        metavar='uniprottfa', 
                        default=os.path.expanduser('~/data/uniprot/uniprot_mouse_all.fasta'),
                        type=str, 
                        help='')

    parser.add_argument('uniprotalt', 
                        metavar='uniprotalt', 
                        default=os.path.expanduser("~/project/hamsini2/uniprot_all_rodents_altcodes.txt"),
                        type=str, 
                        help='')

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    fbase = parse_filebase(args.pairfile)
    pairtfa = f"{fbase}.tfa"
    targettfa = f"{fbase}_targets.tfa"
    phdf = f"{fbase}_phdf.csv"
    evalfile = f"{fbase}_scores.csv"
        
    logging.debug("fbase={fbase} pairtf={pairtfa} targettffa={targettfa} phdf={phdf}")
    logging.debug("uniprottfa={args.uniprottfa} uniprotalt={args.uniprotalt")

    pairlist = parse_dupepairs(args.pairfile)
    df = run_phmmer(pairlist, args.uniprottfa, args.uniprotalt, pairtfa, targettfa)
    df.to_csv(phdf)
    logging.debug(f"Wrote phmmer DF to {phdf}")
    
    make_evaltable(df, pairlist, evalfile )
    
    
    
    
    