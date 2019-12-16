#!/usr/bin/env python
#
# Prints one or more specified GO terms for OBO file to stdout. 
# usage gocat.py  <GOTERM>
#

import argparse
import os


GOFILE=os.path.expanduser('~/data/go/go.obo') 


def main(gotermlist):
    filehandle = open(GOFILE, 'r')
    lines = filehandle.readlines()
    print("read in %d lines" % len(lines))
    print("got arg %s" % gotermlist)
    for gt in gotermlist:
        found = False
        for line in lines:
            if line.startswith("id: %s" % gt):
                found = True
                print('[Term]')
                #print(line.strip())
            if line.startswith("[Term]"):
                found = False
            if found and not line.startswith("[Term]"):
                print(line.strip()) 
        


if __name__ == '__main__':       
    parser = argparse.ArgumentParser()

    parser.add_argument('goterms', 
                        metavar='goterms', 
                        type=str, 
                        help='one or more space-separated goterms GO:XXXXXXXX' ,
                        nargs='*'
                   )
    args= parser.parse_args()

main(args.goterms)