#!/usr/bin/env python
#   Usage: runendcheck.py <basefile> <inputdir> <outputdir> 
#$ -N werner3
#$ -wd $HOME/project/$JOB_NAME
#$ -pe threads 2
#$ -l m_mem_free=2G
#$ -l gpu=0 
#$ -o  $HOME/project/$JOB_NAME/logs/$JOB_NAME.o$JOB_ID.$TASK_ID
#$ -e  $HOME/project/$JOB_NAME/logs/$JOB_NAME.e$JOB_ID.$TASK_ID

import argparse
import logging
import os
import subprocess
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from utils import elzar

FILEBLOCK = 10

def run_samview(infile):
	'''
	Runs sam view. Returns single return line. 
	'''
	command = f"samtools view -c -f 1 {infile}"
	commandlist = command.split()
	logging.debug(f"command is {command}")
	try:
		p = subprocess.run(commandlist, check=True, stdout=subprocess.PIPE, universal_newlines=True)
		output = p.stdout
		lines = output.split('\n')
	except subprocess.CalledProcessError:
		logging.warning(f"Problem with {infile}")
	return output.strip()


def parse_basefile(basefile):
	f = open(basefile)
	flist = f.readlines()
	f.close() 
	return flist
 

def run_all(basefile, indir, outdir):
	taskid = elzar.get_taskid()
	logging.debug(f"taskid is {taskid}")

	lowidx = (taskid - 1) * FILEBLOCK 
	highidx = taskid * FILEBLOCK
	logging.debug(f"lowidx={lowidx} highidx={highidx} ")
	blist = parse_basefile(basefile)
	logging.debug(f"basefile raw len={len(blist)}")	
	blist = blist[lowidx:highidx]
	logging.debug(f"basefile processed len={len(blist)}")
	donelist = []
	for bf in blist:
		bf = bf.strip()
		infile=f"{indir}/{bf}"
		try:
			out = run_samview(infile)
			logging.debug(f"got {out} for {bf}. next...")
			donelist.append( (bf, out))
		except:
			logging.warning(f"got exception from infile {infile}")
	logging.debug("done. returning donelist.")
	return donelist

def write_output( outfile, tuplelist ):
	logging.debug(f"writing out to {outfile}")
	s = ""
	for (a, b) in tuplelist:
		s += f"{a}  {b}\n"
	f = open(outfile, 'w')
	f.write(s)
	logging.debug("done writing.")
	f.close()
		

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
	
	parser.add_argument('basefile', 
                        metavar='basefile', 
                        type=str, 
                        help='list of filenames')
	
	parser.add_argument('indir',
					metavar='indir', 
                    type=str, 
                    help='input dir')
	
	parser.add_argument('outdir',
					metavar='outdir', 
                    type=str, 
                    help='output dir')
	
	args= parser.parse_args()
	
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	if args.verbose:
		logging.getLogger().setLevel(logging.INFO)
		
	outlist = run_all(args.basefile, args.indir, args.outdir)
	taskid = elzar.get_taskid()
	outfile = f"{args.outdir}/{taskid}_ends.out"
	write_output(outfile, outlist)
	#print(outlist)    
  




