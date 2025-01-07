#!/usr/bin/env python
#
#   wrap aria2c providing <src>  <dst> syntax.  
#
#   aria2c -i URIFILE.txt --continue true  --auto-file-renaming false
#   these args allow safe interruption and re-runs without issues. 
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
import requests
from bs4 import BeautifulSoup

def collect_tree(url, depth=0, max_depth=10, urilist=None, dest=None):
    '''
    Recursively list the directory structure of a website.
    Record uris of every file, and relative path to dest for use by aria2c
    return list uf (uri, localdest) tuples

    '''
    logging.debug(f'handling url={url}')
    if urilist is None:
        urilist = []
    if dest is None:
        dest = './'
        
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes

        soup = BeautifulSoup(response.content, 'html.parser')

        links = soup.find_all('a', href=True)
        logging.debug(f'found {len(links)} links.')
        for link in links:
            href = link['href']
            if href.startswith('/'):
                # upward parent link. 
                pass
            elif href.startswith('?'):
                # ignore dross...
                pass
            else:
                # href not parent. directory or file. 
                if href.endswith('/'):
                    # href is sub-directory
                    logging.debug(' ' * depth + href)
                    if depth < max_depth:
                        urilist  = collect_tree(url + href, depth + 1, max_depth, urilist, dest + href )
                else:
                    # href is file in this directory 
                    #print('  ' * depth + href)
                    uri = url + href
                    urilist.append(( uri, dest ) )
                    logging.debug(f'uri={uri} dest={dest}')
        return urilist

    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")


def handle_transfer(source, dest, urifile, force=False):
    '''
    collect full web tree from source root
    create temporary URI file for aria2c
    run aria2c on URI file
    
    '''
    logging.debug(f'source={source} dest={dest} urifile={urifile} force={force}')
    if os.path.exists(urifile) and not force: 
        logging.warning(f'{urifile} exists, skipping scanning and using.')
    else:   
        urilist = collect_tree(source, dest=dest )
        logging.debug(f'got urilist with {len(urilist)} elements. writing to {urifile}')
        with open(urifile, 'w') as fh:
            for (uri, dest) in urilist:
                fh.write(f'{uri}\n   dir={dest}\n')
    # run aria2c -i uris.txt --continue true --auto-file-renaming false
    logging.info(f'running aria2c with urifile {urifile}...')
    
    cmd = [ 'aria2c',
           '-i', urifile, 
           '--continue','true',
           '--auto-file-renaming','false'
           ]
    try:
        logging.debug(f"command: {' '.join(cmd)}")
        #proc = subprocess.Popen(cmd,    
        #                        stdout=subprocess.DEVNULL,
        #                        stderr=subprocess.DEVNULL,
        #                        close_fds=True)                                    
        run_command_shell(cmd)
    
    except KeyboardInterrupt:
        print('Interrupted by Ctrl-C.')
        sys.exit(0)
    
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))      



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

    parser.add_argument('-f','--force',
                    action="store_true",
                    default=False,
                    help='Overwrite URIs even if it exists (rescanning source)') 

    parser.add_argument('-u','--urifile', 
                    metavar='urifile',
                    required=False,
                    default=os.path.abspath('./uris.txt'), 
                    type=str, 
                    help='File to read/write source URL URIs to.')
 
    parser.add_argument('source' ,
                        metavar='source', 
                        type=str,
                        help='root or file URL to copy')

    parser.add_argument('dest' ,
                        metavar='dest',
                        nargs='?',
                        default='./', 
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
    
    dest = os.path.abspath(os.path.expanduser(args.dest))
    dest = dest + '/'
    os.makedirs(dest, exist_ok=True)
      
    if args.source.endswith('/'):
        source = args.source
    else:
        source = args.source + '/'
            
    handle_transfer(source, dest, urifile=args.urifile, force=args.force )
    #urilist = collect_tree(source, dest=dest ) 
    #for (uri, dest) in urilist:
    #    print(f'{uri}\n   dir={dest} ')


