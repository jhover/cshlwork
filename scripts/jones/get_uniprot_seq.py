#!/usr/bin/env python
#
#   URL   https://www.uniprot.org/uniprot/?query=organism:39947&format=fasta&compress=no
#
#

import argparse
import logging
import requests


def read_idfile(infile):
    idlist = []
    f = open(infile, 'r')
    for line in f:
        id=line.strip()
        if len(id) > 1:
            idlist.append(id)
    f.close()
    logging.debug(f"read file of {len(idlist)} items")
    return idlist

def down_seq(orgid, format='fasta' ):
    logging.debug(f"getting sequences for id {orgid} ... ")
    url = f"https://www.uniprot.org/uniprot/?query=organism:{orgid}&format={format}&compress=no"
    r = requests.get(url, allow_redirects=True)
    #print(r.content)
    logging.debug(f"writing to file {orgid}_sequences.{format} ")
    f = open(f'{orgid}_sequences.{format}','wb')
    f.write(r.content)
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

    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='input file')
   
    args= parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    idlist = read_idfile(args.infile)
    for id in idlist:
        down_seq(id)
    
