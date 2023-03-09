#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import logging

from configparser import ConfigParser

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *
from sra.utils import *

if __name__ == "__main__":
    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug',
                        action="store_true",
                        dest='debug',
                        help='debug logging')

    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        dest='verbose',
                        help='verbose logging')

    parser.add_argument('-c', '--config',
                            action="store",
                            dest='conffile',
                            default=None,
                            help='Config file path [~/git/cshlwork/etc/sra.conf]')

    parser.add_argument('-f', '--fasterq',
                        metavar='fastqrun',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download args with fasterq-dump. e.g. SRR14584407')

    parser.add_argument('-p', '--projects',
                        metavar='proj_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download project files, e.g. SRP151064')

    parser.add_argument('-e', '--experiments',
                        metavar='exp_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download experiment files, e.g. SRX2798918')


    parser.add_argument('-r', '--runs',
                        metavar='run_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download run SRA file, e.g. SRR14584407')

    parser.add_argument('-s', '--samples',
                        metavar='samp_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download project files, e.g. SRP151064')

    parser.add_argument('-o', '--outdir',
                        metavar='outdir',
                        type=str,
                        default=None,
                        help='Outdir [ <cwd> ] ')

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.conffile is not None:
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile)) 
    else:
        cp = get_default_config()
        
    cs = get_configstr(cp)
    logging.debug(f"got config: {cs}")

    logging.debug(f"args: {args}")
    
    if args.runs is not None:
        for run_id in args.runs:
            logging.debug(f'Downloading SRA for run {run_id}')
            download_run_sra(run_id, args.outdir)
        logging.info('Done')


    if args.samples is not None:
        for samp_id in args.samples:
            logging.debug(f'Querying SRA for sample {samp_id}')
            #download_run_sra(run_id, args.outdir)
        logging.info('Done')

    if args.projects is not None:
        for proj_id in args.projects:
            logging.debug(f'Querying SRA for project {proj_id}')
            df = query_project_metadata(proj_id)
            print(df)
            #download_run_sra(run_id, args.outdir)
        logging.info('Done')        











###############################
#   OLD
###############################

if __name__ == 'NONONO':    

    parser.add_argument('-q', '--query',
                        metavar='project_id',
                        type=str,
                        nargs='+',
                        default=None,
                        help='Perform standard query on supplied proj_id')


    parser.add_argument('-p', '--download',
                        metavar='proj_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download (Runs) within project args with wget. e.g. SRR14584407')

    parser.add_argument('-m', '--metadata',
                        metavar='project_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download metadata for args. ')


    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.conffile is not None:
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile)) 
    else:
        cp = get_default_config()
        
    cs = get_configstr(cp)
    logging.debug(f"got config: {cs}")

    logging.debug(f"args: {args}")

    if args.setup:
        s = setup(cp)
        # s.execute()

    if args.query is not None:
        q = Query(cp)
        for pid in args.query:
            q.execute(pid)

    if args.download is not None:
        # start a queue
        dq = Queue() 
        # loop through each project id 
        for proj_id in args.download:
            # get the runs associated with that project
            # run_ids = get_runs_for_project(cp, proj_id)
            
                # download the SRA binary file for the run
            fq = Download(cp)
            fq.execute(proj_id)
            # dq.put(fq.execute(proj_id))

        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(cp.get('download', 'max_downloads'))

        # limit number of jobs 
        for n in range(md): 
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')


    if args.fasterq is not None:
        dq = Queue()
        for proj_id in args.fasterq:
            srr_ids = get_runs_for_project(cp, proj_id)
            for srr in srr_ids:
                fq = FasterqDump(cp, srr)
                dq.put(fq)
            logging.debug(f'created queue of {dq.qsize()} items')
            md = int(cp.get('analyze', 'max_jobs'))

        for n in range(md):
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')


    elif args.metadata is not None:
        for srp in args.metadata:
            df = query_project_metadata(srp)
            exps = list(df['Experiment'].unique())
            logging.debug(f"Got list of {len(exps)} experiments")
            for e in exps:
                print(e)


    
