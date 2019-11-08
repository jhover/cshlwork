#!/usr/bin/env python
#
# Process all .fasta files in a directory, producing
# pairwise similarity statistics to an output file
#
#
#  time water -gapopen 10.0 
#             -gapextend 0.5   
#             -asequence 001R_FRG3G.fasta 
#             -bsequence 002L_FRG3G.fasta 
#             -outfile 001R_FRG3Gx002L_FRG3G.water

#  time needle -gapopen 10.0 
#              -gapextend 0.5   
#              -asequence 001R_FRG3G.fasta 
#              -bsequence 002L_FRG3G.fasta 
#              -outfile 001R_FRG3Gx002L_FRG3G.water

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
import traceback



class CommandRunner(threading.Thread):
    
    def __init__(self, commandlist=[]):
        self.commandlist = commandlist
          
    def run(self):
        pass
            




class PairwiseRun(object):
    
    def __init__(self, filelist, workdir ):
        self.log = logging.getLogger()
        self.filelist = filelist
        self.workdir = os.path.abspath(os.path.expanduser(workdir))
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
            self.log.info("Created workdir %s" % self.workdir)
        self.log.debug("Created PairwiseRun workdir=%s" % self.workdir)


    def makewatercommand(self, f1, f2):
        self.log.debug("water: comparing file %s to file %s" % ( f1, f2))
        f1base = os.path.splitext(os.path.basename(f1))[0]
        f2base = os.path.splitext(os.path.basename(f2))[0]
        outfile = "%s/%sx%s.water" % (self.workdir, f1base, f2base)
        self.log.debug("outfile=%s" % outfile)
        cmdlist = ['time water']
        cmdlist.append( '-gapopen 10.0' )
        cmdlist.append('-gapextend 0.5' )
        cmdlist.append('-asequence %s' % f1 )
        cmdlist.append('-bsequence %s' % f2 )
        cmdlist.append('-outfile %s' % outfile )
        self.log.debug("cmdlist=%s" % cmdlist)
        cmd = ' '.join(cmdlist).strip()
        self.log.debug("command is '%s'" % cmd)
        return cmd
        #self.log.info("Running %s against %s" % (f1, f2) )
        #cp = subprocess.run(cmd, check=True, shell=True)
        #self.log.debug("Completed generating %s" % outfile)
        
    
    def makeneedlecommand(self, f1, f2):
        self.log.debug("needle: comparing file %s to file %s" % ( f1, f2))    


    def handlefiles(self):
        #     
        # Take list of files, run water pairwise: 
        #
        commandlist = []
        for i in range(0,len(filelist)):
            f1 = os.path.relpath(os.path.expanduser(filelist[i]))
            for j in range(i + 1,len(filelist)):
                os.path.relpath(os.path.expanduser(filelist[j]))
                f2 = filelist[j]
                self.log.debug("comparing file %s to file %s" % ( f1, f2))
                c = self.makewatercommand(f1, f2)
                commandlist.append(c)
        self.log.debug("list of %d commands made" % len(commandlist))
        
        numthreads = 7
        threadlist = []
        for i in range(0,numthreads):
            t = CommandRunner()
            threadlist.append(t)
        self.log.debug("Made %d threads to run %d commands" % (len(threadlist), len(commandlist)))
        
        for i in range(0, len(commandlist)):
            touse = i % numthreads 
            threadlist[touse].commandlist.append(commandlist[i])
        
        s = ""
        for i in range(0,numthreads):
            s+= "thread [%d]: %d commands "% (i, len(threadlist[i].commandlist))
        self.log.debug("\n%s" % s)

        
        


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
    
    filelist = args.infiles 
    run = PairwiseRun(filelist, args.workdir)
    run.handlefiles()