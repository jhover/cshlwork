#!/usr/bin/env python
#
# splits multi-fasta file into single, with 

import argparse
import logging
import os
import sys
import traceback

class NumSeqReachedException(Exception):
    pass


class Sequence(object):

    def __init__(self, accession, geneid, sequence=""):
        self.accession = accession
        self.geneid = geneid
        self.sequence = sequence

    def addsequence(self, sequence):
        self.sequence += sequence
        
    def asfasta(self):
        s = ">%s|%s\n%s" % ( self.accession, self.geneid, self.sequence)
        return s


class ProcessingRun(object):
    
    def __init__(self, filelist, workdir, numseq):
        self.log = logging.getLogger()
        self.filelist = filelist
        self.workdir = os.path.expanduser(workdir)
        self.workdir = os.path.relpath(self.workdir)
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
            self.log.info("Created workdir %s" % self.workdir)
        self.numseq = int(numseq)
        self.numoutput = 0

    def outputlast(self, seq):
        if seq is not None:
            logging.debug(seq.asfasta())
            try:
                filename = '%s/%s.fasta' % (self.workdir, seq.geneid)
                fh = open(filename, 'w')
                fh.write(seq.asfasta())
                self.numoutput += 1
            except IOError as ioe:
                print("error opening file %s" % filename)
            
            finally:
                fh.close()
    
    
    def parsefile(self, filehandle):
        current = None
        try:
            for line in filehandle:
                if line.startswith(">"):
                    if current is not None:
                        current.addsequence("\n")
                        self.outputlast(current)
                        if self.numoutput >= self.numseq:
                            raise NumSeqReachedException
                    
                    fields = line.split("|")
                    accession = fields[1]
                    fields2 = fields[2].split()
                    geneid = fields2[0]
                    logging.debug('header: accessno=%s id=%s' % (accession, geneid))
                    current = Sequence(accession, geneid)
                elif line.strip().startswith("#"):
                    pass
                else:
                    s = line.strip()
                    current.addsequence(s)        
        except NumSeqReachedException:
            raise
            
        
        except Exception as e:
            traceback.print_exc(file=sys.stdout)        

    
    def handlefiles(self):     
        for filename in self.filelist:
            filename = os.path.relpath(filename)
                        
            try:
                self.log.debug("opening file %s" % filename)
                filehandle = open(sys.argv[1], 'r')
                self.parsefile(filehandle)
                filehandle.close()
            
            except FileNotFoundError:
                self.log.error("No such file %s" % filename)                
            except NumSeqReachedException:
                self.log.info("Desired number of sequences reached %d" % self.numseq)
    

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
                        default='~/work/cafa4-play/seqwork',
                        help='run-specific workdir [~/work/cafa4-play/seqwork]')

    parser.add_argument('-n', '--numseq', 
                        action="store", 
                        dest='numseq', 
                        default=100,
                        help='number of sequences to process. ')    
    
                  
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    filelist = args.infiles 
    run = ProcessingRun(filelist, args.workdir, args.numseq)
    run.handlefiles()


    