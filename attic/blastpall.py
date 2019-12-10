#!/usr/bin/env python
# 
#
#  swissprot blastp ~3 seconds per sequence. 
#
#  blastp -db "swissprot" 
#         -query <infile> 
#         -evalue 1.00e-58 -outfmt "6 qseqid sseqid pident length evalue bitscore"
#

import argparse
import logging
import os
import subprocess
import sys
import threading
import time
import traceback


class CommandRunner(threading.Thread):
    
    def __init__(self, overwrite=True, *args, **kwrds):
        super(CommandRunner, self).__init__(*args, **kwrds)
        self.log = logging.getLogger()
        # self.commands is a list of tuples  (commandstr, outfilepath)
        self.commands = []
        self.overwrite = overwrite
        
          
    def run(self):
        for (cmd, of) in self.commands: 
            self.log.info("Running cmd='%s' outfile=%s " % (cmd, of))
            if os.path.exists(of) and not self.overwrite:
                self.log.debug("Outfile %s exists. Skipping..." % of)
            else:
                cp = subprocess.run(cmd, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.log.debug("Ran cmd='%s' returncode=%s " % (cmd, cp.returncode))
            
    def __repr__(self):
        s = ""
        s += "[%s]CommandRunner with %d commands. " % (self.name, len(self.commands))
        return s

class BlastpRun(object):
   
    def __init__(self, filelist, workdir, overwrite=False, nthreads=1 ):
        self.log = logging.getLogger()
        self.filelist = filelist
        self.threadlist = []
        self.overwrite = overwrite
        self.nthreads = int(nthreads)
        self.program = program
        self.workdir = os.path.abspath(os.path.expanduser(workdir))
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
            self.log.info("Created workdir %s" % self.workdir)
        self.log.debug("Created BlastpRun workdir=%s" % self.workdir)




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

    parser.add_argument('-t', '--threads', 
                        action="store", 
                        dest='nthreads', 
                        default=1,
                        help='number of threads to use.')

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
    
    logging.info("Got arguments...")      
    run = BlastpRun(args.infiles, args.workdir, args.overwrite, args.nthreads)
    run.dosearch()