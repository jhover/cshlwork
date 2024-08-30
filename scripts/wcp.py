#!/usr/bin/env python
#
#   wrap wget and provide cleaner  rwget  <src> <dst> 
#   to make it more like scp, cp, or mv
#
import argparse
import os
import sys
import time
import logging
import subprocess
from urllib.parse import urlparse

import datetime as dt

from configparser import ConfigParser
from _testcapi import datetime_check_delta

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

class NonZeroReturnException(Exception):
    """
    Thrown when a command has non-zero return code. 
    """


def run_command(cmd):
    """
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    @return 
    
    """
    cmdstr = " ".join(cmd)
    logging.debug(f"command: {cmdstr} running...")
    start = dt.datetime.now()
    cp = subprocess.run( cmd, 
                    text=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT
                    )
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.warn(f"got stderr: {cp.stderr}")
    if cp.stdout is not None:
        logging.debug(f"got stdout: {cp.stdout}")   
    if str(cp.returncode) == '0':
        #logging.debug(f'successfully ran {cmdstr}')
        logging.debug(f'got rc={cp.returncode} command= {cmdstr}')
    else:
        logging.warn(f'got rc={cp.returncode} command= {cmdstr}')
       
        #raise NonZeroReturnException(f'For cmd {cmdstr}')
    return cp

def run_command_shell(cmd):
    """
    maybe subprocess.run(" ".join(cmd), shell=True)
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    logging.debug(f"running command: {cmdstr} ")
    start = dt.datetime.now()
    cp = subprocess.run(" ".join(cmd), 
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)

    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.warn(f"got stderr: {cp.stderr}")
        pass
    if cp.stdout is not None:
        #logging.debug(f"got stdout: {cp.stdout}")
        pass
    if str(cp.returncode) == '0':
        #logging.debug(f'successfully ran {cmdstr}')
        logging.debug(f'got rc={cp.returncode} command= {cmdstr}')
    else:
        logging.warn(f'got rc={cp.returncode} command= {cmdstr}')
        raise NonZeroReturnException(f'For cmd {cmdstr}')
    return cp


  
def handle_urlroot(urlroot, dest):
    dest=os.path.abspath(dest)
    start = dt.datetime.now()
    logging.debug(f'downloading {urlroot} to {dest}')
    o = urlparse(urlroot, scheme='', allow_fragments=True)
    scheme = o.scheme
    hostname = o.hostname
    port = o.port
    urlpath = o.path
    logging.debug(f'scheme={scheme} host={hostname} path={urlpath}')
    
    sleeptime = 15
    keep_going = True    
    progress_score = 0
    total_bytes = 0
    runcount = 0
    #min_run = 3  # run at least 3 times to be sure.
    
    pfields = urlpath.split('/')
    pfields = [x for x in pfields if x != '']
    n_cut = len(pfields) 
    topdir = pfields[-1]
    logging.debug(f'topdir={topdir} n_cut={n_cut} sleeptime={sleeptime}')
    
    outroot = f'{dest}/{topdir}'
    outroot = os.path.abspath(outroot)
    
    if not os.path.exists(outroot):
        os.makedirs(outroot, exist_ok=True)
        logging.debug(f'confirmed outroot={outroot}')
    else:
        logging.info(f'outroot exists. Check output...')
        progress_score = score_dest(outroot)
            
    while keep_going:
        runcount += 1
        logging.info(f'[run {runcount}] running wget to {dest}')
        
        copy_cmd = ['wget',
               '--recursive',
               '--continue',
               f'--cut-dirs={n_cut}',  
               '--no-parent',
               '-l', '10',
               '-w', '0.5',   
               '--no-host-directories',
               '-R', 'index.html*',
               f'--directory-prefix={outroot}', 
               urlroot ]

        try:
            logging.info(f'Running wget. run number: {runcount} ')
            logging.debug(f"command: {' '.join(copy_cmd)}")
            proc = subprocess.Popen(copy_cmd,    
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    close_fds=True)                                    
        except NonZeroReturnException as nzre:
            logging.error(f'problem with {infile}')
            logging.error(traceback.format_exc(None))        

        process_running = True
        while process_running:
            logging.debug(f'sleeping {sleeptime} seconds... ')
            time.sleep(sleeptime)                  
            rc = proc.poll()
            if rc is not None:
                logging.debug(f'process has exited rc={rc} ')
                logging.info(f'wget has finished. Checking output progress...')
                old_progress_score = progress_score
                (progress_score, size_kbytes) = score_dest(outroot)
                if old_progress_score < progress_score:
                    logging.info(f'New output found! old:{old_progress_score} < new:{progress_score} Run again...')
                    logging.debug(f'{size_kbytes} downloaded so far.')                
                else:
                    logging.info(f'No new output found. old:{old_progress_score} == new:{progress_score}. Stop.')
                    logging.debug(f'{size_kbytes} downloaded.')
                    keep_going = False
                process_running = False
            else:
                logging.debug(f'process still running.')
        end = dt.datetime.now()
        delta = end - start
        ts = delta.total_seconds()
        logging.info(f'Done. Elapsed: {format_interval(ts)} Downloaded: {format_storage(size_kybtes)}')

def format_interval(delta):
    days, remainder = divmod(delta, 86400 )
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{days}d:{hours}h:{minutes}m:{seconds}s'    

def format_storage(total_kbytes):
    '''
    assuming 1k block counts from du -ks
    '''
    tb, remainder = divmod(total_kbytes, 1073741824 )
    gb, remainder = divmod(remainder, 1048576)
    mb, kbytes = divmod(remainder, 1024)
    return f'{tb}T:{gb}G:{mb}M:{kbytes}k'
        

def score_dest(destroot):
    size_bytes = space_check(destroot)
    file_count = filecount_check(destroot)
    dest_score = size_bytes * file_count
    logging.debug(f'got dest_score={dest_score} for {file_count} entries consuming {size_bytes} bytes')
    return (dest_score, size_bytes)
    

def space_check(destroot):
    spacecheck_cmd = ['du',
       '-ks', destroot 
       ]
    try:
        scp = run_command(spacecheck_cmd)
        size_bytes = int(scp.stdout.split('\t')[0])
        return size_bytes
        
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))        
        raise 

def filecount_check(destroot):    
    filecount_cmd = ['find',
       '-L', destroot 
       ]
    try:
        fcp = run_command(filecount_cmd)
        file_count = len(fcp.stdout.split('\n'))
        return file_count
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))        
        raise 



def wrap_command(cmd):
    
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))        
        raise 

       

def handle_urlfile(urlfile, dest):
    logging.debug(f'downloading urls in  {urlfile} to {dest}')


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
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

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')

 
    parser.add_argument('source' ,
                        metavar='source', 
                        type=str,
                        help='root or file URL to copy')

    parser.add_argument('dest' ,
                        metavar='dest', 
                        type=str,
                        help='destination path ')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)
    
    logging.debug(f'urlroot={args.source} dest={args.dest}')
    
    handle_urlroot(args.source, args.dest)
        


