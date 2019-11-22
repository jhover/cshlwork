import logging
import os
import logging
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess


def read_phmmer_table(filename):
    df = pd.read_table(filename, 
                     names=['target','t-acc','query','q-acc',
                            'e-value', 'score', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                            'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                     skip_blank_lines=True,
                     comment='#',
                     index_col=False,
                     skiprows=3,
                     engine='python', 
                     sep='\s+')
    df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'], axis=1)
    df['db'] = df.apply(lambda row: row.target.split('|')[0], axis=1)
    df['tacc'] = df.apply(lambda row: row.target.split('|')[1], axis=1)
    df['prot_spec'] = df.apply(lambda row: row.target.split('|')[2], axis=1)
    df['protein'] =   df.apply(lambda row: row.prot_spec.split('_')[0], axis=1)
    df['species'] =   df.apply(lambda row: row.prot_spec.split('_')[1], axis=1)
    return df



def run_phmmer_files(filelist):
#
#  time phmmer --tblout 7955.phmmer.2.txt 
#              --cpu 16 
#              --noali 
#              ~/data/cafa4/TargetFiles/sp_species.7955.tfa 
#              ~/data/uniprot/uniprot_sprot.fasta 
#              > 7955.phmmer.console.out 2>&1
    
    dbase = "/data/hover/data/uniprot/uniprot_sprot.fasta"
    for file in filelist:
        cmd = _make_phmmer_cmdline(file, dbase)
        cp = subprocess.run(cmd, 
                            shell=True, 
                            universal_newlines=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        logging.debug("Ran cmd='%s' returncode=%s " % (cmd, cp.returncode))

        
def _make_phmmer_cmdline(filename, database):
    outpath = os.path.dirname(filename)
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outfile = "%s.phmmertbl.txt" % filebase
    #self.log.debug("outfile=%s" % outfile)
    cmdlist = ['time', 'phmmer']
    cmdlist.append( '--tblout  %s ' % outfile )
    cmdlist.append('--noali' )
    cmdlist.append(' %s ' % filename )
    cmdlist.append(' %s ' % database )
    cmd = ' '.join(cmdlist).strip()
    #self.log.debug("command is '%s'" % cmd)
    return (cmd, outfile)
    
    
    
    