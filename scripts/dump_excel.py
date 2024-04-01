#!/usr/bin/env python
#

import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd
import scipy

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

FILTERMAP = {
         "&" : "and",
         ":" : ".",
         ")" : "",
         "(" : "",
         ";" : ".",
         "," : "",
         "`" : "",
         '"' : "",
         "'" : "",
         ' ' : '_',
         '_' : '_',
         '[' : "",
         ']' : "",
         #'-' : ".",
         '..': '.',
         '+': '_',
         }

def normalize_name(aname):
    '''
    Fix name to make nice columns
    '''

    for k in FILTERMAP.keys():
        aname = aname.replace( k, FILTERMAP[k])
    aname = aname.lower()
    return aname


def dump_xlsx(infile, outdir, merge=False, normalize=True):
    '''
    Create TSVs for all sheets in infile Excel file in the outdir. 
    
    Merge assumes all sheets have same columns. Adds column with sheet source. 
    
    ''' 
    infile = os.path.abspath(infile)
    dirname = os.path.dirname(infile)
    filename = os.path.basename(infile)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'dirname={dirname} filename={filename} base={base} ext={ext}')
    
    if not os.path.exists(outdir):
        outdir = os.path.abspath(outdir)
        os.makedirs(outdir, exist_ok=True)
        logging.debug(f'made outdir={outdir}')

    logging.debug(f'dumping {infile} to {outdir}')
    
    alldfs = []
    if infile.endswith('.xlsx'):
        edf = pd.read_excel(infile, sheet_name=None)
        slist = list(edf.keys())
        if len(slist) < 2:
            # one-sheet file, name sheet after filename. 
            sdf = edf[slist[0]]
            if normalize:
                sdf.columns = [ normalize_name(x) for x in sdf.columns ]            
            outfile = f'{outdir}/{base}.tsv'
            sdf.to_csv(outfile, sep='\t')
            logging.info(f'wrote single-sheet file to {outdir}/{base}.tsv')
        # Note merge is ignored for single-sheet files. 
                
        else:
            # 2 or more sheets in file.
            for sheet in slist:
                sdf = edf[sheet]
                if normalize:
                    sdf.columns = [ normalize_name(x) for x in sdf.columns ]
                outfile = f'{outdir}/{sheet}.tsv'
                sdf.to_csv(outfile, sep='\t')
                logging.info(f'wrote {outdir}/{sheet}.tsv')            
                if merge:
                    sdf['sheet'] = normalize_name(sheet) 
                    alldfs.append(sdf)

            if merge:
                mergedf = pd.concat(alldfs, axis=0)
                mergedf.reset_index(inplace=True, drop=True)
                outfile = f'{outdir}/{base}.merged.tsv'
                mergedf.to_csv(outfile, sep='\t')
                logging.info(f'wrote {outfile}') 
    else:
        logging.warn(f'{infile} does not end xlsx.')





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

    parser.add_argument('-m', '--merge',
                        default=False , 
                        action="store_true", 
                        dest='merge', 
                        help='concatenate all sheets, assuming columns are the same.')

    parser.add_argument('-n', '--normalize',
                        default=True , 
                        action="store_true", 
                        dest='normalize', 
                        help='normalize column names (lowercase, remove spaces, special chars.')


    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=True,
                    default=None, 
                    type=str, 
                    help='outdir. built from .xlsx name if not given.')     

    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single .xlsx file.')
       

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   


    logging.debug(f'infile={args.infile}')
       
    outdir = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infile
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
        
    # create and handle 'real' 'spikein' and 'normalized' barcode matrices...
    dump_xlsx(args.infile, outdir=args.outdir, merge=args.merge)
    