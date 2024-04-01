#!/usr/bin/env python
#

import argparse
import logging
import os
import pprint
import re
import sys
import traceback

from configparser import ConfigParser

import pandas as pd
import scipy

from scipy.io import loadmat


gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)


def list_vars(infile, config=None):
    annots = loadmat(infile)
    alist = list(annots.keys())
    varnames = []
    for a in alist:
        if not a.startswith('__'):
            varnames.append(a)
    return varnames

def handle_chars(data):
    # pull sequences for index. 
    charlist = []
    for row in data:
        s = ''.join([chr(item) for item in row])
        charlist.append(s)
    return charlist


def get_vars(infile, config=None):
    '''

    ''' 
    annots = loadmat(infile)
    alist = list(annots.keys())
    varnames = []
    for a in alist:
        if not a.startswith('__'):
            varnames.append(a)
    
    annot_dict = {}
    
    for v in varnames:
        data = annots[v]
        logging.debug(f'var={v} shape={data.shape}')
        
        if len(data.shape) < 3:
            # assume 2-dimensional shape
            lol = [ ]
            yrows, xcols = data.shape
            for y in range(0, yrows):
                rlist = []
                for x in range(0, xcols):
                    rlist.append( data[y][x][0] )
                lol.append(rlist)
            df = pd.DataFrame(lol)
                    
            annot_dict[v] = df 
        else:
            logging.warning(f"don't know how to dump 3+ dimensional data.")
    return annot_dict
        
        
def write_var(df, varname, outdir):
    outfile = f'{outdir}/{v}.tsv'
    df.to_csv(outfile, sep='\t')
    logging.debug(f'wrote DF to {outfile}')




def snippets():
    
    outfile = f'{outdir}/refbarcodes.txt'
    writelist(outfile, refseqlist)
    rbcdf = pd.DataFrame(refseqlist, columns=['sequence'])
    outfile = f'{outdir}/refbarcodes.fasta'
    write_fasta_from_df(rbcdf, './refbarcodes.fasta')
    
    barcodematrix = annots['barcodematrix']
    logging.debug(f'bcmatrix shape={barcodematrix.shape}')    
    bdf = pd.DataFrame(barcodematrix, index=refseqlist)
    logging.debug(f'bcmatrix columns={bdf.columns} converting to ')
    bdf.columns = list( [str( x + 1)  for x in list(bdf.columns)] )
   

def format_config(cp):
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    s = pprint.pformat(cdict, indent=4)
    return s


def check_outdir(infile, outdir=None):
    if outdir is not None:
        outdir = os.path.abspath(outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        filepath = os.path.abspath(infile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    return outdir
   
   
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
                        default=os.path.expanduser('~/git/cshlwork/etc/matlab.conf'),
                        type=str, 
                        help='config file.')    

    subparsers = parser.add_subparsers( dest='subcommand',
                                        help='sub-command help.')

################################ subparsers  ########################################

    parser_list = subparsers.add_parser('list',
                                          help='list contents of .mat')
    
    parser_list.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single MATLAB .mat file.')

    parser_dump = subparsers.add_parser('dump',
                                          help='dump contents of .mat')
    
    parser_dump.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single MATLAB .mat file.')

    parser.add_argument('-V','--varlist', 
                        metavar='varlist',
                        required=False,
                        default=None,
                        type=str, 
                        help='comma separated list of variable names to dump. ')

    parser_dump.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = format_config(cp)    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infile={args.infile}')
       
    if args.subcommand == 'list':
        varnames = list_vars(args.infile, config=cp)
        for v in varnames:
            print(v)

    if args.subcommand == 'dump':
        outdir = check_outdir(args.infile, args.outdir)
        annot_dict = get_vars(args.infile, config=cp)
        logging.debug(f'got annot_dict len={len(annot_dict)} keys={list(annot_dict.keys())}')
        
        if args.varlist is None:
            varlist = list(annot_dict.keys())
        else:
            # accept comma or quoted whitespace-separated list elements on command line
            varlist = [ x.strip() for x in re.split(',|\s',args.varlist ) ]
            varlist = [ x for x in varlist if len(x) > 0 ]
        
        for v in varlist:
            logging.debug(f'handling variable {v}')
            df = annot_dict[v]
            write_var(df, v, outdir )
            
            
            
      
       
    # create and handle 'real' 'spikein' and 'normalized' barcode matrices...
    #dump_matlab(cp, args.infile, sampdf, outdir, expid=args.expid, label=args.label )
    