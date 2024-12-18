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


  
def handle_urlroot(urlroot, dest, delay=0.1):
    dest=os.path.abspath(dest)
    delay=float(delay)
    start = dt.datetime.now()
    
    logging.debug(f'downloading {urlroot} to {dest}')
    o = urlparse(urlroot, scheme='', allow_fragments=True)
    scheme = o.scheme
    hostname = o.hostname
    port = o.port
    urlpath = o.path
    logging.debug(f'scheme={scheme} host={hostname} path={urlpath}')
    
    # size is checked very size_check_interval * sleeptime seconds. 
    # be careful if adjusted that du can finish. 
    sleeptime = 10
    size_check_interval = 12 
    
    keep_going = True    
    progress_score = 0
    prerun_kbytes = 0
    runcount = 0
    loopcount = 0
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
        logging.debug(f'created outroot={outroot}')
    else:
        logging.info(f'outroot={outroot} exists. Check output size...')
        (progress_score, size_kbytes) = score_dest(outroot)
        prerun_kbytes = size_kbytes   
    
    while keep_going:
        runcount += 1
        logging.info(f'runcount={runcount} {urlpath} -> {outroot}')
        
        copy_cmd = ['wget',
               '--recursive',
               '--continue',
               f'--cut-dirs={n_cut}',  
               '--no-parent',
               '-e robots=off',
               '-l', '12',
               '-w', f'{delay}',   
               '--no-host-directories',
               '-R', 'index.html*',
               f'--directory-prefix={outroot}', 
               urlroot ]

        try:
            logging.debug(f"wget run:{runcount} command: {' '.join(copy_cmd)}")
            proc = subprocess.Popen(copy_cmd,    
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    close_fds=True)                                    
        
        except KeyboardInterrupt:
            print('Interrupted by Ctrl-C.')
            sys.exit(0)
        
        except NonZeroReturnException as nzre:
            logging.error(f'problem with {infile}')
            logging.error(traceback.format_exc(None))        

        process_running = True
        logging.info(f'Subprocess poll each {sleeptime}s. Progress check each {sleeptime * size_check_interval}s. ')
        try:
            while process_running:
                loopcount += 1
                #logging.debug(f'sleeping {sleeptime} seconds... ')
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
                    #logging.debug(f'process still running.')
                    if loopcount % size_check_interval == 0:
                        (lps, lsk) = score_dest(outroot)
                        if lsk > 0:
                            lsm = lsk / 1024
                            logging.info(f'Process running. {lsm:.2f}MB so far...')
                        else:
                            logging.info(f'Process running, but no download so far...')
        
        except KeyboardInterrupt:
            print('Interrupted by Ctrl-C.')
            sys.exit(0)
             
        # Finalize...
        end = dt.datetime.now()
        delta = end - start
              
        logging.info(f'Elapsed: {format_interval(delta)} Total: {format_storage(size_kbytes)} ')
        logging.info(f'Downloaded: {format_storage(size_kbytes - prerun_kbytes)} Speed: {format_downspeed(size_kbytes - prerun_kbytes, delta)} MB/s   ')

    

def format_downspeed(kbytes, delta):
    '''
    takes kbytes, timedelta and gives MB/s
    '''
    mbs = 0.0
    if kbytes > 0.0 and delta.total_seconds() > 0:
        mbytes = kbytes / 1024
        mbs = mbytes / delta.total_seconds()
    return f'{mbs} MB/s'
    

def format_interval(delta):
    '''
    delta is standard datetime delta object.
    '''
    days, remainder = divmod(delta.total_seconds(), 86400 )
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{days}d:{hours}h:{minutes}m:{seconds:.1f}s'    

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

    parser.add_argument('-i','--interval', 
                        metavar='delay_interval',
                        required=False,
                        default=0.1,
                        type=float, 
                        help='Delay interval between web requests. Seconds.')

 
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
    
    logging.debug(f'urlroot={args.source} dest={args.dest} delay={args.interval}')
    
    handle_urlroot(args.source, args.dest, args.interval)
        


