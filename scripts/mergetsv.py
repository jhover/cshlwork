#!/usr/bin/env python
#
#  Takes N .tsv file and outputs merged tsv. 
#
import os
import sys
#gitpath=os.path.expanduser("~/git/cshlwork")
#sys.path.append(gitpath)

import pandas as pd
import numpy as np
import traceback

if not len(sys.argv) > 1:
    print("need input file arg(s)")
    sys.exit()

topdf = None  
  
infiles = sys.argv[1:]
for infile in infiles:
    basename = os.path.basename(infile)
    # print(f"{infile}")
    try:
        df = pd.read_csv(infile, index_col=None, sep='\t', comment="#")
        if topdf is None:
            topdf = pd.DataFrame(columns=list(df.columns))   
        topdf = topdf.append(df, ignore_index=True, sort=True)
    
    except:
        print('something went wrong')

topdf.drop_duplicates(inplace=True)
topdf = topdf.reset_index(drop=True)
topdf.to_csv(sys.stdout, sep='\t', index=False)


