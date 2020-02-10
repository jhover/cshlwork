#!/usr/bin/env python
#
#  Takes filenames on input and renames them, replacing species codes with NCBI taxids.
#
#
import os
import sys

gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

from fastcafa.fastcafa import *

replacelist = [ 'ARATH' ,'BACCR','BOVIN','CAEEL','DANRE','DICDI','DROME','ECOLI',
               'HUMAN','MAIZE','MOUSE','MYCGE','PIG','PSEAI','RAT','SALTY','SCHPO','YEAST' ]

config = get_default_config()
smo = get_specmap_object(config)
replacedict = {}
for scode in replacelist:
    tc = smo.get_taxonid(scode)
    print(tc)
    replacedict[scode] = tc 

#print(smo)

flist = sys.argv[1:]
for filepath in flist:
    #print(filepath)
    dirpath = os.path.dirname(filepath)
    fixname = os.path.basename(filepath)
    for s in replacelist:
        fixname = fixname.replace(s, replacedict[s], 1)
        #print(f"fixed filename: {fixname}")
    #print(f"oldname is {filepath}")
    newpath = f"{dirpath}/{fixname}"
    #print(f"newpath is {newpath}")
    print(f"{filepath} -> {newpath}")
    os.rename(filepath, newpath)

    