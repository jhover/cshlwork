#!/usr/bin/env python
#
# Process all .fasta files in a directory, producing
# pairwise similarity statistics to an output file
#   
# Needleman-Wunsch global sequence alignment
# Comparing homology/proteins with similar function  in similar-sized 
# sequenes
#
#  time needle -gapopen 10.0 
#              -gapextend 0.5   
#              -asequence 001R_FRG3G.fasta 
#              -bsequence 002L_FRG3G.fasta 
#              -outfile 001R_FRG3Gx002L_FRG3G.water

#
#  Finding conserved domains or motifs. 
#
#  time water -gapopen 10.0 
#             -gapextend 0.5   
#             -asequence 001R_FRG3G.fasta 
#             -bsequence 002L_FRG3G.fasta 
#             -outfile 001R_FRG3Gx002L_FRG3G.water

#
#  Esprit http://www.ijbcb.org/ESPRITPIPE/php/download.php
#  Won tests of PSA tools. ICBR Interdisciplinary Centor for Biotechnology Research
#    University of Florida. 
#  ??
# https://www.majordifferences.com/2016/05/difference-between-global-and-local.html
#
#
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
            self.log.debug("Running cmd='%s' outfile=%s " % (cmd, of))
            if os.path.exists(of) and not self.overwrite:
                self.log.debug("Outfile %s exists. Skipping..." % of)
            else:
                cp = subprocess.run(cmd, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.log.debug("Ran cmd='%s' returncode=%s " % (cmd, cp.returncode))

            
    def __repr__(self):
        s = ""
        s += "[%s]CommandRunner with %d commands. " % (self.name, len(self.commands))
        return s


class PairwiseRun(object):
    
    def __init__(self, filelist, workdir, overwrite=False ):
        self.log = logging.getLogger()
        self.filelist = filelist
        self.threadlist = []
        self.overwrite = overwrite
        self.workdir = os.path.abspath(os.path.expanduser(workdir))
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
            self.log.info("Created workdir %s" % self.workdir)
        self.log.debug("Created PairwiseRun workdir=%s" % self.workdir)


    def makewatercommand(self, f1, f2):
        #self.log.debug("water: comparing file %s to file %s" % ( f1, f2))
        f1base = os.path.splitext(os.path.basename(f1))[0]
        f2base = os.path.splitext(os.path.basename(f2))[0]
        outfile = "%s/%sx%s.water" % (self.workdir, f1base, f2base)
        #self.log.debug("outfile=%s" % outfile)
        cmdlist = ['time', 'water']
        cmdlist.append( '-gapopen 10.0' )
        cmdlist.append('-gapextend 0.5' )
        cmdlist.append('-asequence %s' % f1 )
        cmdlist.append('-bsequence %s' % f2 )
        cmdlist.append('-outfile %s' % outfile )
        #self.log.debug("cmdlist=%s" % cmdlist)
        cmd = ' '.join(cmdlist).strip()
        #self.log.debug("command is '%s'" % cmd)
        return (cmd, outfile)
        #self.log.info("Running %s against %s" % (f1, f2) )
        #cp = subprocess.run(cmd, check=True, shell=True)
        #self.log.debug("Completed generating %s" % outfile)
        
    
    def makeneedlecommand(self, f1, f2):  
        #self.log.debug("water: comparing file %s to file %s" % ( f1, f2))
        f1 = os.path.abspath(f1)
        f2 = os.path.abspath(f2)
        
        f1base = os.path.splitext(os.path.basename(f1))[0]
        f2base = os.path.splitext(os.path.basename(f2))[0]
        outfile = "%s/%sx%s.needle" % (self.workdir, f1base, f2base)
        #self.log.debug("outfile=%s" % outfile)
        cmdlist = ['time', 'needle']
        cmdlist.append( '-gapopen 10.0' )
        cmdlist.append('-gapextend 0.5' )
        cmdlist.append('-asequence %s' % f1 )
        cmdlist.append('-bsequence %s' % f2 )
        cmdlist.append('-outfile %s' % outfile )
        #self.log.debug("cmdlist=%s" % cmdlist)
        cmd = ' '.join(cmdlist).strip()
        #self.log.debug("command is '%s'" % cmd)
        return (cmd, outfile)

    def makecommands(self):
        #     
        # Take list of files, run water pairwise: 
        #
        commandlist = []
        for i in range(0,len(self.filelist)):
            f1 = os.path.relpath(os.path.expanduser(self.filelist[i]))
            for j in range(i + 1,len(self.filelist)):
                os.path.relpath(os.path.expanduser(self.filelist[j]))
                f2 = self.filelist[j]
                #self.log.debug("comparing file %s to file %s" % ( f1, f2))
                c = self.makeneedlecommand(f1, f2)
                commandlist.append(c)
        self.log.debug("commandlist of %d commands made" % len(commandlist))
        
        numthreads = 4
        for i in range(0,numthreads):
            t = CommandRunner(name=str(i), overwrite=self.overwrite)
            self.threadlist.append(t)
        self.log.debug("Made %d Runners to run %d commands" % (len(self.threadlist), len(commandlist)))
        
        for i in range(0, len(commandlist)):           
            touse = i % numthreads
            c = commandlist[i]
            usethread = self.threadlist[touse] 
            usethread.commands.append(c)
        
        s = ""
        for i in range(0, numthreads):
            #s+= "thread [%d]: %d commands "% (i, len(threadlist[i].commands))
            s+= str(self.threadlist[i])
        self.log.debug("%s" % s)

    def runcommands(self):
        for t in self.threadlist:
            t.start()
        
        for t in self.threadlist:
            t.join()


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

    parser.add_argument('infiles', 
                        metavar='infiles', 
                        type=str, 
                        nargs='+',
                        help='a list of .fasta sequence files')
    
    parser.add_argument('-o', '--overwrite', 
                        action="store_true", 
                        dest='overwrite',
                        default=False, 
                        help='redo commands that have already created output')
    
    parser.add_argument('-w', '--workdir', 
                        action="store", 
                        dest='workdir', 
                        default='~/work/cafa4-play/seqout',
                        help='run-specific workdir [~/work/cafa4-play/seqout]')
                   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    
    run = PairwiseRun(args.infiles, args.workdir, args.overwrite)
    run.makecommands()
    run.runcommands()
    
    