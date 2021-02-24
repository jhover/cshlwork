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

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)
from protlib.uniprot import *
from protlib.phmmer import *


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
    logging.debug(f" processed {lnum} lines. skipped {knum} NAs. produced {len(dupelist)} items in dupelist[0] = {dupelist[0]}")
    #logging.debug(f"dupelist: {dupelist}")
    return dupelist




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

def make_evaltable(pdf, pairlist, evalfile ):
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


def split_pairlist(pairlist):
    qlist = []
    tlist = []
    for (q, t) in pairlist:
        qlist.append(q)
        tlist.append(t)
    return (qlist, tlist)


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
    querytfa = f"{fbase}.tfa"
    targettfa = f"{fbase}_targets.tfa"
    phdf = f"{fbase}_phdf.csv"
    evalfile = f"{fbase}_scores.csv"
        
    logging.debug(f"fbase={fbase} querytfa={querytfa} targettffa={targettfa} phdf={phdf}")
    logging.debug(f"uniprottfa={args.uniprottfa} uniprotalt={args.uniprotalt}")

    pairlist = parse_dupepairs(args.pairfile)
    (querylist, targetlist) = split_pairlist(pairlist)
    logging.debug(f"qlist[:2] = {querylist[:2]} tlist[:2] = {targetlist[:2]} ")
    

#    df = run_phmmer(pairlist, args.uniprottfa, args.uniprotalt, querytfa, targettfa)
#    df.to_csv(phdf)
#    logging.debug(f"Wrote phmmer DF to {phdf}")
#    make_evaltable(df, pairlist, evalfile )
    
    
    
    
    