#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import logging
import time

from configparser import ConfigParser

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import run_command_shell, NonZeroReturnException

#  Example gcloud storage URL
#      gs://fc-secure-f2e71e6a-3081-4816-98dc-a21c7e498971/c8d604ec-8e89-4ee6-8bd5-0e546481fcaf/Optimus/a47f9e7e-4aa3-44b7-afca-68404cf714b4/call-MergeStarOutputs/pBICCNsMMrACAACAiF019d210630A4_sparse_counts_row_index.npy
# 


def download_gcloud_url(gsurl, outdir=None, force = False):
    '''
    
    '''
    f_exists = False    
    if outdir is None:
        outdir = os.path.abspath('./')
    outfile = gsurl.split('/')[-1:][0]
    logging.debug(f'downloading {gsurl} -> {outdir}')        
    outpath = f'{outdir}/{outfile}'
    if os.path.exists():
        logging.warning(f'file {outpath} exists. ')
        f_exists = True
    else:
        f_exists = False
        logging.debug(f'file {outpath} does not exist. ')
    
    cmd = [ 'gsutil','cp',
            gsurl,
            outdir   
           ]
    if not f_exists or force:
        logging.debug(f'running gsutil command...')
        time.sleep(5)
        try:
            run_command_shell(cmd)
        except NonZeroReturnException as nzre:
            logging.error(f'problem with gsurl {gsurl}')
            logging.error(traceback.format_exc(None))
            raise    

def handle_sample_set(setid, urllist, outdir):
    '''
     downloads all files in urllist in subdirectory {outdir}/{setid}/
    '''
    outdir = os.path.expanduser(outdir)
    outdir = os.path.abspath(outdir)
    outpath = f'{outdir}/{setid}'
    logging.debug(f'handling {len(urllist)} files -> {outpath}')
    os.makedirs(outpath, exist_ok = True)
    for url in urllist:
        download_gcloud_url(url, outpath)


def parse_sample_text(infile):
    lines = []
    with open(infile) as fh:
        for line in fh.readlines():
            fields = line.split()
            print(fields)
            lines.append(fields)    
    logging.debug(f'got {len(lines)} samples in file {infile}')
    return lines
    



# to do: include drivers for 10x and ss alignments
if __name__ == "__main__":


    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
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

    parser.add_argument('-i', '--infiles',
                        metavar='infiles',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='infile with one line per sample <setid> <gsurl> [<gsurl> ... ]')

    parser.add_argument('-o', '--outdir',
                        metavar='outdir',
                        type=str,
                        default=os.path.abspath("./"),
                        help='Outdir. ')

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    logging.debug(f"args: {args}")
    
    if args.infiles is not None:
        for infile in args.infiles:
            lol = parse_sample_text(infile)
            for flist in lol:
                handle_sample_set(flist[0], flist[1:], args.outdir)
            
        logging.info('Done')
    

    
