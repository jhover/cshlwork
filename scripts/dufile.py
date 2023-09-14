#!/usr/bin/env python
#
#
# Provides total storage usage of all input files
#
#
import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import traceback



def sum_storage(infiles, human_readable=False):
    retval = None
    total_bytes = 0
    for infile in infiles:
        basename = os.path.basename(infile)
        logging.debug(f"checking '{infile}' ...")
        size_bytes = os.path.getsize(infile)
        total_bytes = total_bytes + size_bytes
    
    if human_readable:
        #units = ['','K','M','G','T']
        units = ['T','G','M','K','']
        
        size_b = size_bytes
        size_k = size_bytes/1024
        size_m = size_bytes/(1024*1024)
        size_g = size_bytes/(1024*1024*1024)
        size_t = size_bytes/(1024*1024*1024*1024)

        #size_list = [size_b,size_k,size_m,size_g,size_t ]
        size_list = [size_t,size_g,size_m,size_k,size_b ]
        logging.debug(f"{size_b} {size_k}K {size_m}M {size_g}G {size_t}T")
        for i, size in enumerate(size_list):
            #print(f'{i} -> {size}')
            if size > 1:
                retval = f'{size:.0f}{units[i]}'
                break
    else:
        retval = total_bytes

    return retval
 

if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    parser = argparse.ArgumentParser(add_help=False)  
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('-h', '--human', 
                        action="store_true", 
                        dest='human', 
                        help='human readable')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='out file. stdout if not given.')  
    
    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='TSV [ TSV TSV ...] ')

    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)       
    
    print(sum_storage(args.infiles, human_readable = args.human))
        