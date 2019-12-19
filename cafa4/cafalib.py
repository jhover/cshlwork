#!/usr/bin/env python
#
#     Top level  Runplugins      Infoplugins
#     CAFA4Run
#                PhmmerPlugin                         
#                OrthologPlugin                       looks up ortholog info
#                                QuickGOPlugin        gets GO info
#                                UniprotGOPlugin      gets GO info
#                Expression
#
# 
#      ************************CANONICAL DATAFRAME COLUMNS*********************************
#   
# COLUMN        DESCRIPTION               MAPPINGS                             EXAMPLES
# cafaid        cafa4 target identifier   N/A                                  T2870000001
# cafaprot      cafa4 target protein id                                        1433B
# cafaspec      cafa4 target protien species                                   MOUSE
# proteinid     UniProtKB:entry/accession  quickgo:gene product_id             P63103
# protein       all caps name                                                  1433B
# gene          Free-text gene name.                                           Lrrk2  Ywahb
# geneid        Gene name+species.                                             LRRK2_MOUSE     
# taxonid      NCBI taxon id                                                  9606                 
# species       all caps code                                                  MOUSE   PONAB
# goterm        annotated term                                                 GO:0005634
# goaspect      biological process|molecular function|cellular component       cp    bp
# goevidence    evidence codes for GO annotation.                              IEA 
# evalue        BLAST/HMMER/PHMMER expect statistic                                               1.000000e-126
# bias          Adjustement to score for char prevalence                       3.5
# score         BLAST/HMMER/PHMMER bit-score                                   400.3
# db            
# probest       Probability estimate for prediction.                           0.68  [.000001-1.0]  
#
#

__author__ = "John Hover"
__copyright__ = "2019 John Hover"
__credits__ = []
__license__ = "Apache 2.0"
__version__ = "0.99"
__maintainer__ = "John Hover"
__email__ = "hover@cshl.edu"
__status__ = "Testing"

import argparse
from configparser import ConfigParser
import logging
import os
import sys
import tempfile

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

gitpath=os.path.expanduser("~/git/cshl-work/cafa4")
sys.path.append(gitpath)
from ontology import GeneOntologyGOInfoPlugin
import quickgo
import uniprot

class CAFA4Run(object):
    
    def __init__(self, config, targetfile, name=None):
        '''
        Embodies all the processing for a single run against all targets.
        Overall input is a set of Target sequence files. 
        Overall output is a properly-formatted CAFA4 prediction file.   
        
        '''
        self.config = config
        self.targetfile = targetfile
        if name is None:
            self.name = config.get('global','name')
        else:
            self.name = name
            self.config.set('global','name',self.name)
        self.log = logging.getLogger('CAFA4Run')
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.author = config.get('global','author')
        self.modelnumber = config.get('global','modelnumber')
        self.keywords = config.get('global','keywords')
        self.profile = config.get('global','profile')
        self.pipeline = [ x.strip() for x in config.get( self.profile,'pipeline').split(',')]
        self.outbase = config.get( self.profile,'outbase')
    
    def __repr__(self):
        s = "CAFA4Run:"
        for atr in ['name', 'outdir','targetfile','pipeline']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s


    def cafafile(self, dataframe):
        '''
        Produce properly-formated CAFA submission file:
        
        E.g. filename:    gillislab_1_10090_go.txt
                          gillislab_1_78_go.txt   
        
        HOVER          GILLIS_LAB
        MODEL          1
        KEYWORDS       orthologs, phmmer
        ACCURACY       1  PR=0.86; RC=0.30
        T100900000001  GO:0042203    .85
        T100900000002  GO:0003998    .03
        .
        .
        .
        
        '''
        cafafile = "%s/gillislab_1_%s_go.txt" % (self.outdir,  "287")
        self.log.debug("Opening cafafile=%s" % cafafile)
        f = open( cafafile, 'w' )
        
        s = ""
        s += "%s\t%s\n" % (self.author, 'GILLIS_LAB')
        s += "MODEL\t%s\n" % 1
        s += "KEYWORDS\tortholog, gene expression\n"
        s += "ACCURACY\t1\tPR=%f;\tRC=%f\n" % (1.00, 1.00) 
        self.log.debug("dataframe columns=%s" % dataframe.columns )
        for row in dataframe.iterrows():
            target = row[1]['cafaid']
            goterm = row[1]['goterm']
            probest = float(row[1]['probest'])  # 1.00 for call/no-call ->  precision/recall *point*. 
            s += "%s\t%s\t%f\n" % (target, goterm, probest)
        s+="END\n"
        
        f.write(s)
        f.close()
        return s
    
    def execute(self):
        self.log.info("Begin run...")
        
        phm = PhmmerPlugin(self.config, self.targetfile)
        self.log.debug(phm)
        df = phm.execute()
        
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/%s-%s-phmmer.csv" % (self.outdir, self.name,  self.outbase))
        
        ortho = OrthologPlugin(self.config)
        self.log.debug(ortho)
        df = ortho.execute(df)
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/%s-%s-ortho.csv" % (self.outdir, self.name, self.outbase))
        
        
        go = GOPlugin(self.config)        
        df = go.execute(df)
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/%s-%s-go.csv" % (self.outdir, self.name, self.outbase))        
        
        # not needed if we use quickgo. contains both evidence codes and go_aspect
        #go = GOPlugin(self.config)
        #df = go.execute(df)
        #self.log.info("\n%s" % str(df))
        
        cfstr = self.cafafile(df)
        lines = cfstr.split('\n')
        summary = '\n'.join(lines[:10] + ['.', '.','.']  + lines[-10:]) 
        self.log.debug("\n%s" % summary)
        self.log.info("Ending run...")


class CAFAPlugin(object):
    '''
    A plugin that normally takes a Pandas dataframe as input, executes, and produces another DF.
    Returns a DF from disk cached if regen=False
    '''
    
    REPR_ATTRS=['outdir']
    
    def __init__(self, config):
        self.config = config
        self.kname = self.__class__.__name__
        self.log = logging.getLogger(self.kname)
        self.outdir = os.path.expanduser( config.get('global','outdir') )

    def __repr__(self):
        s = "%s:" % self.kname
        for atr in self.__class__.REPR_ATTRS:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s 

    def execute(self, dataframe):
        self.log.warn("You have called a base Plugin. No change.")
        return dataframe
 
 

class PhmmerPlugin(CAFAPlugin):
    '''
    Pipeline object. Takes list of Fasta files to run on, returns pandas dataframe of 
    similar sequences with score. 
    Input:   List of FASTA files 
    Output:  Pandas DataFrame
    '''
    
    REPR_ATTRS=['targetfile','database', 'score_threshold','cpus']

    def __init__(self, config, targetfile):
        super(PhmmerPlugin, self).__init__(config)
        self.targetfile = targetfile
        self.configname = self.__class__.__name__.lower()
        self.database = os.path.expanduser( config.get(self.configname, 'database') )
        self.score_threshold = config.get(self.configname,'score_threshold')
        self.cpus = config.get(self.configname,'cpus')
        

    def execute(self):
        self._get_targetinfo()
        outfile = self.run_phmmer_file(self.targetfile , self.database)
        self.log.debug("phmmer outfile=%s" % outfile)
        df = self.read_phmmer_table(outfile)
        #self.log.debug("dflist is %s" % outdfs)
        #df = pd.concat(outdfs, ignore_index=True)
        self.log.debug(str(df))
        return df
            
    def _get_targetinfo(self):
        '''
        Reads self.targetlist files, pulls out 
        target entry info 
        Returns dict of list.
           <targetid> : [<cafaprotein>, <cafaspecies> ]
        
        '''
        self.targetinfo={}
        filehandle = None
        try:
            filehandle = open(self.targetfile, 'r')
            for line in filehandle:
                if line.startswith('>'):
                    (targetid , prot_spec) = line.split()
                    targetid = targetid[1:]
                    (cafaprot, cafaspec) = prot_spec.split('_')
                    self.targetinfo[targetid] = [cafaprot, cafaspec]
                else:
                    pass
               
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        finally:
            if filehandle is not None:
                filehandle.close()
        
        self.log.debug("Parsed file(s) with %d targets" % len(self.targetinfo) )
        return self.targetinfo
    
    
    def run_phmmer_files(self, targetlist, database="/data/hover/data/uniprot/uniprot_sprot.fasta"):
    #
    #  time phmmer --tblout 7955.phmmer.2.txt --cpu 16  --noali 
    #              ~/data/cafa4/TargetFiles/sp_species.7955.tfa 
    #              ~/data/uniprot/uniprot_sprot.fasta 
    #              > 7955.phmmer.console.out 2>&1
        
        dbase = database
        outlist = []
        for file in targetlist:
            if not os.path.exists(file):
                raise FileNotFoundError()
            (cmd, outfile) = self._make_phmmer_cmdline(file)
            self.log.debug("Running cmd='%s' outfile=%s " % (cmd, outfile))
            cp = subprocess.run(cmd, 
                                shell=True, 
                                universal_newlines=True, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            
            outlist.append(outfile)
            self.log.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd,outfile, cp.returncode))
        return outlist

    def run_phmmer_file(self, targetfile, database="/data/hover/data/uniprot/uniprot_sprot.fasta"):
      
        dbase = database
        if not os.path.exists(targetfile):
            raise FileNotFoundError()
        (cmd, outfile) = self._make_phmmer_cmdline(targetfile)
        self.log.debug("Running cmd='%s' outfile=%s " % (cmd, outfile))
        cp = subprocess.run(cmd, 
                            shell=True, 
                            universal_newlines=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        self.log.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd,outfile, cp.returncode))
        return outfile

            
    def _make_phmmer_cmdline(self, filename):
        outpath = os.path.dirname(filename)
        filebase = os.path.splitext(os.path.basename(filename))[0]
        outfile = "%s/%s.phmmer.tbl.txt" % (self.outdir, filebase)
        #self.log.debug("outfile=%s" % outfile)
        cmdlist = ['time', 'phmmer']
        cmdlist.append( '--tblout  %s ' % outfile )
        cmdlist.append('--noali' )
        cmdlist.append('--cpu %s ' % self.cpus)
        if self.score_threshold is not None:
            cmdlist.append("-T %s " % self.score_threshold)
        cmdlist.append(' %s ' % filename )
        cmdlist.append(' %s ' % self.database )
        cmd = ' '.join(cmdlist).strip()
        #self.log.debug("command is '%s'" % cmd)
        return (cmd, outfile)

    def read_phmmer_table(self, filename):
        self.log.debug("Reading %s" % filename)
        df = pd.read_table(filename, 
                         names=['target','t-acc','cafaid','q-acc',
                                'evalue', 'score', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                                'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                         skip_blank_lines=True,
                         comment='#',
                         index_col=False,
                         skiprows=3,
                         engine='python', 
                         sep='\s+')
        self.log.debug(str(df))
        self.log.debug("Dropping unneeded columns..")
        df = df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 
                 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 
                 'description'] , axis=1)
        self.log.debug("Parsing compound fields to define new columns...")
        self.log.debug("Splitting first field for db")
        df['db'] = df.apply(lambda row: row.target.split('|')[0], axis=1)
        self.log.debug("Splitting first field for target accession")
        df['proteinid'] = df.apply(lambda row: row.target.split('|')[1], axis=1)
        self.log.debug("Splitting first field for prot_species")
        df['prot_spec'] = df.apply(lambda row: row.target.split('|')[2], axis=1)
        self.log.debug("Splitting protein_species to set protein")
        df['protein'] =   df.apply(lambda row: row.prot_spec.split('_')[0], axis=1)
        self.log.debug("Splitting protein_species to set species")
        df['species'] =   df.apply(lambda row: row.prot_spec.split('_')[1], axis=1)
        self.log.debug("Dropping split columns...")
        df.drop(columns=['target','prot_spec'], axis=1, inplace=True)
        self.log.debug("Adding target metadata from input files")
        df['cafaprot']= df.apply(lambda row:  self.targetinfo[row.cafaid][0], axis=1)      
        df['cafaspec']= df.apply(lambda row:  self.targetinfo[row.cafaid][1], axis=1)          
        
        return df
    

class OrthologPlugin(CAFAPlugin):
    '''
    Pipeline object. Takes Pandas dataframe
    looks up orthologs and GO codes.  
    Input:  Pandas DataFrame
    Output: Pandas DataFrame
    '''
    REPR_ATTRS=['outdir','backend']

    def __init__(self, config):
        '''
        
        '''
        super(OrthologPlugin, self).__init__(config)
        self.configname = self.__class__.__name__.lower()
        self.backend = config.get(self.configname ,'backend').strip()
        self.exc_ec_list = [i.strip() for i in config.get(self.configname, 'excluded_evidence_codes').split(',')] 
        self.uniprot = uniprot.UniProtGOPlugin(self.config)
        self.quickgo = quickgo.QuickGOPlugin(self.config)

        
    def execute(self, dataframe):
        '''
        for each row of dataframe, look up ortholog in uniprot and for each GO code
        add a new row with gene, goterm, gocategory
        
        iterate input df fully, putting new info in new df. 
        merge old + new df, return resulting dataframe
        
        '''
        self.log.info("Looking up each ortholog")

        newdf = None
        
        if self.backend == 'uniprot': 
            self.log.debug("Calling uniprot back end.")
            newdf = self.uniprot.get_df(dataframe)
            
        if self.backend == 'quickgo':
            self.log.debug("Querying QuickGo for %d unique entries" % len(entries))
            udf = self.quickgo.get_df(entries)
            self.log.debug(qdf)
            qdf.to_csv("%s/quickgo.csv" % self.outdir)
        
        self.log.debug("\n%s" % str(newdf))
        return newdf
            
    
class GOPlugin(CAFAPlugin):
    '''
    Pipeline object. 
    Takes Pandas dataframe 
        1) For each proteinid, looks up GO aspect and evidenceid.
        2) Walks up GO tree, fanning out each line with different GOterm, *increasing* likelihood/prob estimate.  

    Input:  Pandas DataFrame
    Output: Pandas DataFrame
        Added rows
        Added columns: goaspect  goevidence cafapest

    '''
    REPR_ATTRS=['outdir']

    def __init__(self, config):
        '''
        
        '''
        super(GOPlugin, self).__init__(config)
        self.configname = self.__class__.__name__.lower()
        self.go = GeneOntologyGOInfoPlugin(self.config)


    def execute(self, dataframe):
        gdf = self.go.get_df()
        self.log.debug("\n%a" % str(gdf))
        
        dataframe['probest'] = 1.00
        igdf = gdf.set_index('goterm')
        gdict = igdf.to_dict('index')
        # now indexed by goterm   gdict[goterm] -> {'name': 'osteoblast differentiation', 'namespace': 'bp'}
        # rd['GO:0001649']['namespace']  
                
        dataframe['goaspect'] = dataframe.apply(
            lambda row: gdict[row.goterm]['namespace'], 
            axis=1)
        self.log.debug(str(dataframe))
        return dataframe
        

#def get_plugin(klassname):
#    return getattr(sys.modules[__name__], klassname)
   
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
                        help='a .fasta sequence files')
    
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default='~/etc/cafa4.conf',
                        help='Config file path [~/etc/cafa4.conf]')

    parser.add_argument('-n', '--name', 
                        action="store", 
                        dest='runname', 
                        default='default',
                        help='Run-specific identifier to use in file output.')    
                    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
           
    c4run = CAFA4Run(cp, args.infile, args.runname)
    logging.debug(c4run)
    c4run.execute()
