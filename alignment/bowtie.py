#!/usr/bin/env python
# 
import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser
from collections import defaultdict

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import numpy as np
import pandas as pd

from cshlwork.utils import run_command_shell, NonZeroReturnException, setup_logging

#
# Possibly set up to do pipline internally?
# See:  https://www.youtube.com/watch?v=_cTbADrGLCQ  
#
#    ConfigParser entries:
#
#    [bowtie]
#    threads = 10
#
#    [bowtie2]
#    threads = 10
#

BOWTIE_1_COLS=['name_read', 'strand','name_align','offset','seq','quals','ceil','mm_desc']

#
#  first row is mandatory fields. 
#   second row are optional expected fields. 
BOWTIE_2_COLS=['name_read', 'flagsum', 'name_align','offset','qual', 'cigar', 'mate', 'mate_offset', 'fraglen', 'seq', 'quals', 
                  'score', 'next', 'n_amb', 'n_mismatch', 'n_gap', 'n_gapext', 'distance', 'md', 'yt' ]
BOWTIE_2_INT_COLS= ['flagsum', 'offset', 'qual', 'mate_offset', 'fraglen', 'score', 'n_amb', 'n_mismatch', 'n_gap', 'n_gapext', 'distance', ]

BOWTIE_OPT_COLS= ['score', 'next', 'n_amb', 'n_mismatch', 'n_gap', 'n_gapext', 'distance', 'md', 'yt' ]
OPT_MAP = { 'score'     : 'AS',
            'next'      : 'XS',
            'n_amb'     : 'XN',
            'n_mismatch': 'XM',
            'n_gap'     : 'XO',
            'n_gapext'  : 'XG',
            'distance'  : 'NM', 
            'md'        : 'MD',
            'yt'        : 'YT'    
            }



def fix_columns_int(df, columns):
    '''
    forces column in dataframe to be an integer. NaNs become '0'
    Only floating points can be NaN. No good solution for integers...
    '''   
    for col in columns:
        try:
            logging.debug(f'trying to fix col {col}')
            fixed = np.array(df[col], np.int16)
            logging.debug(f'fixed=\n{fixed}')
            df[col] = fixed
                
        except ValueError:
            logging.debug(f'invalid literal in {col}')
    return df



def run_bowtie(config, qfile, rfile, outfile, tool='bowtie', force=False):
    '''
    bowtie-build -q BC1.seq.fasta indexes/BC1.bt 
    bowtie -v 3 -p 10 -f --best -a indexes/BC1.bt BC1_seq.fasta BC1.bt.algn
   
    
   
    '''
    logging.debug(f'running bowtie {qfile} -> {rfile} => {outfile}')
     
    qfilepath = os.path.abspath(qfile)    
    qdirname = os.path.dirname(qfilepath)
    qfilename = os.path.basename(qfilepath)
    (qbase, qext) = os.path.splitext(qfilename)   
    logging.debug(f'handling query {qfilepath}')

    rfilepath = os.path.abspath(rfile)    
    rdirname = os.path.dirname(rfilepath)
    rfilename = os.path.basename(rfilepath)
    (rbase, rext) = os.path.splitext(rfilename)   
    logging.debug(f'handling reference {rfilepath}')    
    
    if tool == 'bowtie':
        nidx = ''
    elif tool == 'bowtie2':
        nidx = '2'
    
    idxdir = os.path.abspath(f'{rdirname}/bt{nidx}idx')
    os.makedirs(idxdir, exist_ok = True )
    idxpfx = f'{idxdir}/{rbase}'

    # don't run second if first fails. 
    build_ok = False

    threads = config.get(tool, 'threads')
   
    if tool == 'bowtie':
        idxoutput = f'{idxpfx}.1.ebwt'
        output_exists = os.path.exists(idxoutput)
        logging.info(f'output {idxoutput} exists={output_exists}')
        if not output_exists or force:
            cmd = ['bowtie-build',
                   #'-q',
                   '--threads', threads , 
                   rfile,
                   idxpfx, 
                   ]
            logging.debug(f'running {cmd}')
            try:
                run_command_shell(cmd)
                build_ok = True
                
            except NonZeroReturnException as nzre:
                logging.error(f'problem with infile {qfile}')
                logging.error(traceback.format_exc(None))
                raise     
        
            logging.debug(f'bowtie-build done.')
        else:
            build_ok = True

        max_mismatch = config.get(tool, 'max_mismatch')
        cmd = ['bowtie',
               '-v', max_mismatch,
               '--threads', threads ,  # # threads
               '-f',      # -f query input files are (multi-)FASTA .fa/.mfa
               '--best',
               '-a',      # -a/--all report all alignments; very slow, MAPQ not meaningful
               idxpfx,
               qfile,
               outfile
               ]
        
        if build_ok: 
            logging.debug(f'running {cmd}')

            try:
                run_command_shell(cmd) 
            except NonZeroReturnException as nzre:
                logging.error(f'problem with infile {qfile}')
                logging.error(traceback.format_exc(None))
                raise
        else:
            logging.warn(f'Unable to run, bowtie{nidx}-build not OK. ')     
        logging.debug(f'bowtie{nidx} done.')
    
   
    elif tool == 'bowtie2':
        idxoutput = f'{idxpfx}.1.bt2'
        output_exists = os.path.exists(idxoutput)
        logging.info(f'output {idxoutput} exists={output_exists}')
        
        if not output_exists or force:
            cmd = ['bowtie2-build',
               '--threads', threads , 
               rfile,
               idxpfx, 
               ]
    
            logging.debug(f'running {cmd}')
            try:
                run_command_shell(cmd)
                build_ok = True
                
            except NonZeroReturnException as nzre:
                logging.error(f'problem with infile {qfile}')
                logging.error(traceback.format_exc(None))
                raise     
        else:
            build_ok = True

        if build_ok:
            cmd = ['bowtie2',
                   #'-N', '3',
                   '--threads', threads ,    # # threads
                   '-f',       # -f query input files are (multi-)FASTA .fa/.mfa
                   #'--best',
                   '--all',   # -a/--all report all alignments; very slow, MAPQ not meaningful
                   '-x', idxpfx,
                   '-S', outfile, 
                   qfile,
                   ]
            logging.debug(f'running {cmd}')
            try:
                run_command_shell(cmd)
            
            except NonZeroReturnException as nzre:
                logging.error(f'problem with infile {infile}')
                logging.error(traceback.format_exc(None))
                raise                
            logging.debug(f'bowtie{nidx} done.')
        else:
            logging.warn(f'unable to run, bowtie{nidx}-build not OK. ')
        
    return outfile



def make_bowtie_df(infile, max_mismatch=3):
    with open(infile) as f:
        line=f.readline()
    if line.startswith('@HD'):
        logging.info('Detected bowtie2 input.')
        df = make_bowtie2_df(infile)
        logging.debug(f'df before max_mismatch =< {max_mismatch}: {df}')
        df = df[df['n_mismatch'] <= max_mismatch]
        # alignments to *other* sequences only

        logging.debug(f'df after max_mismatch < {max_mismatch} =\n{df}')
    else:
        logging.debug('Detected bowtie1 input.')
        df = make_bowtie1_df(infile)
        df['n_mismatch'] = max_mismatch
    return df


def make_bowtie1_df( infile):
    '''
    parse standard bowtie output and create pandas dataframe with standard columns
    
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)    
    bdf = pd.read_csv(filepath, sep='\t',header=None, names=BOWTIE_1_COLS)
    return bdf


def make_bowtie2_df(infile):
    '''
    bt1: name_read  strand  name_align  offset  seq  quals   ceil    mm_desc
    
    bt2: name_read         name_align  offset  seq  quals   
            algn_score  next_score  n_amb  n_mm  n_gaps n_gapext  distance      
    
    parse bowtie2 output format.
    produce dataframe.     
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)
    outfile = os.path.join(dirname, f'{base}.btdf')
    logging.debug(f'handling base={base}  ->  {outfile}')    

    try:
        logging.debug(f" attempting to open '{infile}'")
        filehandle = open(infile, 'r')
    except FileNotFoundError:
        logging.error(f"No such file {infile}")
        raise                
    
    current = 0
    sumreport = 1
    suminterval = 1000000
    repthresh = sumreport * suminterval
    
    # list of lists to hold data
    lol = []
    
    #
    def def_value(): 
        return "None"
    
    # BOWTIE_2_COLS=['name_read', 'flagsum', 'name_align','offset', 'qual', 'cigar', 'mate', 'mate_offset', 'fraglen', 'seq', 'quals', 
    # 'score', 'next_score', 'n_amb', 'n_mismatch', 'n_gaps', 'n_gapext', 'distance','md','yt' ]
    try:
        while True:
            line = filehandle.readline()
            if line == '':
                break
            if line.startswith("@HD"):
                pass
            elif line.startswith("@SQ"):
                pass                   
            elif line.startswith("@PG"):
                pass
            else:
                allfields = [x.strip() for x in line.split('\t')]
                mands = allfields[0:11]
                flist = mands
                optfields = allfields[11:]
                #logging.debug(f'allfields length={len(allfields)} mands len={len(mands)} opts len={len(optfields)}')
                #logging.debug(f'mands=\n{mands}')
                #logging.debug(f'opts=\n{optfields}')
                optdict = defaultdict(def_value)
                for of in optfields:
                    ofields = of.split(':')
                    key = ofields[0]
                    val = ofields[2]
                    optdict[key] = val 
                    #logging.debug(f'handling optkey {key} = {val}')
                for colname in BOWTIE_OPT_COLS:
                    mapid = OPT_MAP[colname]
                    val = optdict[mapid]
                    #logging.debug(f'for colname={colname} map={mapid} val={val}')
                    flist.append(val)
                logging.debug(f'appending flist={flist}')
                
                lol.append(flist)
                current += 1                        
                if current >= repthresh:
                    logging.info(f"Processed {current} entries... ")
                    sumreport +=1
                    repthresh = sumreport * suminterval
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    if filehandle is not None:
        filehandle.close()          

    logging.debug(f'making dataframe from LOL len={len(lol)}')
    df = pd.DataFrame(data=lol, columns=BOWTIE_2_COLS )
    df = df[df['name_align'] != '*']
    df.reset_index(inplace=True, drop=True)
    logging.debug(f'made initial dataframe={df}')
    df = fix_columns_int(df, BOWTIE_2_INT_COLS)
    return df



def make_adjacency_df(bowtiedf):
    '''
    consume custom bowtie dataframe (from all x all alignment) and
    create adjacency matrix dataframe
    '''
    labels = np.unique(btdf[['name_read','name_align']])
    sdf = btdf.filter( ['name_read','name_align'], axis=1 )
    sdf['val'] = 1
    mdf = sdf.pivot(index = 'name_read', 
                    columns='name_align', 
                    values='val').reindex(columns=labels, index=labels, fill_value=0)
    mdf.fillna(0, inplace=True)
    return mdf


    

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
    
    parser.add_argument('-c','--config', 
                        metavar='config',
                        required=False,
                        default=os.path.expanduser('~/git/cshlwork/etc/bowtie.conf'),
                        type=str, 
                        help='config file.')    

    parser.add_argument('-a','--aligner', 
                    metavar='aligner',
                    required=False,
                    default='bowtie2', 
                    type=str, 
                    help='aligner tool  [bowtie | bowtie2]')

    parser.add_argument('-m','--max_mismatch', 
                        metavar='max_mismatch',
                        required=False,
                        default=3,
                        type=int, 
                        help='Max mismatch for aligner read collapse.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Output bowtie DF file')  

    parser.add_argument('query',
                        type=str,
                        help='fasta query file')

    parser.add_argument('reference',
                        type=str,
                        help='fasta reference file')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'query={args.query}  ref={args.reference}')  

    #  config, qfile, rfile, outfile, tool='bowtie'
    filepath = os.path.abspath(args.query)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    btfile = os.path.join(dirname, f'{base}.{args.aligner}') 
       
    btoutfile = run_bowtie(cp, args.query, args.reference, btfile, tool=args.aligner)
    logging.info(f'produced {btoutfile}')
    btdf = make_bowtie_df(btoutfile, max_mismatch= args.max_mismatch)
    logging.info(f'produced DF=\n{btdf}')
    outfile = args.outfile
    if outfile is None:
        outfile= os.path.join(dirname, f'{base}.btdf.tsv')
    btdf.to_csv(outfile, sep='\t') 
    
    

    
