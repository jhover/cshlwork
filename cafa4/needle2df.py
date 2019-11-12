#!/usr/bin/env python
# generate dataset/network graph from needle sequence similarity output
# 
# Reads in multiple pairwise files, creates a single pandas dataframe, stores to disk
#
#      GENE1     GENE2    identity  similarity  score
# 0   <NAME1> x <NAME2>    .016      .020       9.0
# 1   <NAME1> x <NAME2>    .321      .392       32.0
# 2   <NAME1> x <NAME2>    .121      .185       29.2
#
# Length: 665
# Identity:      25/665 ( 3.8%)
# Similarity:    42/665 ( 6.3%)
# Gaps:         572/665 (86.0%)
# Score: 24.0
#
#
#  To load saved data:
#   import pandas as pd
#   import pyarrow as pa
#   import pyarrow.parquet as pq
#
#   t = pq.read_table('<filepath>')
#   df = t.to_pandas()
#
#


import argparse
import logging
import os
import sys
import traceback

import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


class NeedleParse(object):
    
    def __init__(self, filelist, outfile):
        self.log = logging.getLogger()
        self.filelist = filelist
        self.outfile = os.path.expanduser(outfile)
        self.outfile = os.path.abspath(self.outfile)
        self.data=[]


    def handlefiles(self):     
        for filename in self.filelist:
            filename = os.path.abspath(filename)
            # 003R_FRG3Gx002L_FRG3G.needle
            base = os.path.basename(filename)
            desc = os.path.splitext(base)[0]                        
            try:
                (alabel, blabel) = desc.split('x')
                self.log.debug("parsed out a and b: %s %s" % (alabel,blabel))
                self.log.debug("opening file %s" % filename)
                filehandle = open(filename, 'r')
                dtup = self.parsefile(filehandle, alabel, blabel)
                self.data.append(dtup)
                filehandle.close()
            except ValueError:
                self.log.error("Problem parsing base filename: %s " % base )
            except FileNotFoundError:
                self.log.error("No such file %s" % filename)   
    
    
    def parsefile(self, filehandle, alabel, blabel):
        current = None
        length = None
        ident = None 
        simil = None 
        score = None 
        
        try:
            for line in filehandle:
                if line.startswith("# Length: "):
                    length = int(line.split()[2])
                    self.log.debug("length is %d" % length)
                
                if line.startswith("# Identity: "):    
                    (n,d) = line.split()[2].split('/')
                    ident = int(n) / int(d)
                    self.log.debug("identity is %f" % ident)
                
                if line.startswith("# Similarity: "):
                    (n,d) = line.split()[2].split('/')
                    simil = int(n) / int(d)
                    self.log.debug("similarity is %f" % simil)
                #if line.startswith("# Gaps: "):
                
                if line.startswith("# Score: "):
                    score = float(line.split()[2])
                    self.log.debug("score is %f" % score)   
            ntup = (alabel, blabel, length, ident, simil, score)
            self.log.debug("ntup=%s" % str(ntup))
            return ntup
        except Exception as e:
            traceback.print_exc(file=sys.stdout)        

    def getdf(self):
        '''
        Return pandas dataframe of all data. 
        '''
        df = pd.DataFrame(self.data, columns=['genea' ,'geneb','length','identity','similarity','score'])
        return df

    def save(self):
        '''
        Save as pandas dataframe in msgpack format to outfile. 
        
        Load it back in via:
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
table = pq.read_table('run1.par')
df = table.to_pandas()
        
        '''
        self.log.debug("Saving file...")
        #df.to_msgpack(args.outfile)
        df = self.getdf()
        table = pa.Table.from_pandas(df)
        pq.write_table(table, self.outfile)
        

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s')
    
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
                        dest='infiles', 
                        type=str,
                        required=False, 
                        nargs='+',
                        help='a list of .fasta sequence files')
       
    parser.add_argument('-w', '--workdir', 
                        action="store", 
                        dest='workdir', 
                        default='~/work/cafa4-play/seqout',
                        help='run-specific workdir [~/work/cafa4-play/seqout]')

    parser.add_argument('-o', '--outfile', 
                        action="store", 
                        dest='outfile', 
                        help='outfile')
    
    parser.add_argument('-L', '--filelist', 
                        action="store", 
                        dest='filelist', 
                        help='file containing listing of files to process (to avoid directory/shell limits)')
    
                   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if args.filelist is not None:
        logging.debug("Found filelist opt. ")
        f = open(args.filelist, 'r')
        args.infiles = [x.strip() for x in f]
        f.close()

    n2d = NeedleParse(args.infiles, args.outfile)
    n2d.handlefiles()
    df = n2d.getdf()
    print(df)
    n2d.save()
    
    
    
