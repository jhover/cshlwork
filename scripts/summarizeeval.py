#!/usr/bin/env python
#
#  Takes eval file and outputs various useful info...
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

import pandas as pd
import numpy as np

from fastcafa.fastcafa import *

if not len(sys.argv) > 1:
    print("need input file arg")
    sys.exit()
    
infile = sys.argv[1]
# print(f"{infile}")
try:
    df = pd.read_csv(infile, index_col=0, comment="#")
    #print(f"{df}")
except FileNotFoundError:
    print(f"no such file {infile}")
    sys.exit()
    
totaltrue = len(df[df.correct == True])
totalnum = df.shape[0] 
accuracy = totaltrue / totalnum
numcids = len(df.cid.unique() )
meanf1max = df.groupby('cid')['f1max'].max().mean()  


print(f"num cafaids: {numcids}")
print(f"meanf1max: {meanf1max}")
print(f"accuracy: {accuracy}")
