#!/usr/bin/env python

import logging
import pandas as pd
import os 
import sys

FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
logging.basicConfig(format=FORMAT)
logging.getLogger().setLevel(logging.DEBUG)

# gene_upmap
d1 = [ { 'db'  : 'tr', 
         'acc' : 'A0A1S3RID5',
         'uid' : 'A0A1S3RID5_SALSA',
         'gn'  : 'LOC106602976', 
        },
        { 'db'  : 'tr', 
         'acc' : 'B5XFF4',
         'uid' : 'B5XFF4_SALSA',
         'gn'  : 'WRB', 
        },
        { 'db'  : 'tr', 
         'acc' : 'B9EPD1',
         'uid' : 'B9EPD1_SALSA',
         'gn'  : 'MLP3C', 
        }       
      ]

# seq_upmap

d2 = [ { 'seqid' : 'G803000050759',
         'uid'   : 'A0A1S3RID5_SALSA',
       },
       { 'seqid' : 'G803000050760',
         'uid'   : 'B5XFF4_SALSA',
       },
       { 'seqid' : 'G803000050761',
         'uid'   : 'B9EPD1_SALSA',
       },
    ]

# predout
d3 = [ { 'seqid'    : 'G803000050759',
         'goterm'   : 'GO:0005667',
         'prob'     : 0.11
       },
       { 'seqid' : 'G803000050759',
         'goterm'   : 'GO:0043966',
         'prob'     : 0.12
       },
       { 'seqid' : 'G803000050761',
         'goterm'   : 'GO:0045893',
         'prob'     : 0.13
       },
    ]

geneupmapdf = pd.DataFrame(d1)
logging.debug(f"\n{geneupmapdf}")
sequpmapdf = pd.DataFrame(d2)
logging.debug(f"\n{sequpmapdf}")
predoutdf = pd.DataFrame(d3)
logging.debug(f"\n{predoutdf}")

#cdf = pd.merge(cdf1, cdf2 , how='outer', on=['cid','cgid','goterm'] )
fixedpredout = pd.merge(predoutdf, sequpmapdf,  how='left', on=['seqid'])
fixedpredout.drop(['seqid'], inplace=True, axis=1)
logging.debug(f"fixed pred out is \n{fixedpredout}")





