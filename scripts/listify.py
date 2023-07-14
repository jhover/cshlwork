#!/usr/bin/env python
#
# Turn input into a JSON/Python list. 
# Input is stdin or file
# Output has no spaces, for use in bash
#
#
import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)


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
        
    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='*',
                        default=None, 
                        help='infile [ infile2, infile3] ')

    args= parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)           
    
    #print(args.infiles)
    
    outlist = []
    if len(args.infiles) < 1:
        # no infiles, use stdin
        for line in sys.stdin:
            fields = [x.strip() for x in line.split()]
            for f in fields:
                outlist.append(f)
            
    else:
        for fn in args.infiles:
            with open(fn) as fh:
                for line in fh.readlines():
                    fields = [x.strip() for x in line.split()]
                    for f in fields:
                        outlist.append(f)                    
    #print(outlist)
    outstr = str(outlist)
    outstr = outstr.replace(" ","")
    print(outstr) 
    
