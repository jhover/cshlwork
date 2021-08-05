#!/usr/bin/env python
#
#  Merge all GTEx tissue-specific eqtl datasets into one dataframe, with 
# additional column for tissue.  
#
#  Typical ~30 minutes on 12G data. 

import logging
import os
import sys
import traceback

import pandas as pd
import numpy as np

if not len(sys.argv) > 1:
    print("need input file arg(s)")
    sys.exit()


topdf = None  

dflist = []


outfile = 'all.signif_variant_gene_pairs.csv'

dirname = "/"  

infiles = sys.argv[1:]
for infile in infiles:
    basename = os.path.basename(infile)
    dirname = os.path.dirname(infile)
    print(f'filename = {dirname}/{basename}')
    tissue = basename.split('.')[0]
    colname = tissue.lower()
    print(tissue)
    try:
        df = pd.read_table(infile, header=0)
        df['tissue'] = colname
        logging.debug(df)
        if topdf is None:
            topdf = df   
        else:
            topdf = topdf.append(df, ignore_index=True, sort=True)
    
    # Gracefully exit on Ctrl-C
    except(KeyboardInterrupt):
        print("\nCaught Ctrl-C. Bye.")
        sys.exit(0)
    
    except Exception as e:
            traceback.print_exc(file=sys.stdout)

print('splitting variant id...')
new = topdf['variant_id'].str.split('_',2, expand=True)
print('creating new columns...')
topdf['chr'] = new[0] 
topdf['var_pos'] = new[1]
print('splitting gene_id...')
new = topdf['gene_id'].str.split('.',1,expand=True)
topdf['gene'] = new[0]
topdf['gene_ver'] = new[1]
print(f'writing final to file {outfile}')
            
topdf.to_csv(outfile)
    
    
