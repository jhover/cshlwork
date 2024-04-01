#!/usr/bin/env python
#
#  generic library to handle pickle-based caching...
#
import logging
import os
import pickle
import sys

CACHEDIR = os.path.expanduser('~/tmp')
   
def ensure_cachedir():
    logging.debug(f'checking for cachedir {CACHEDIR}')
    if not os.path.exists(CACHEDIR):
        logging.debug(f'creating cachedir {CACHEDIR}')
        os.makedirs(CACHEDIR)
    else:
        logging.debug('cachedir exists.')


def cache_object_exists(key):
    ensure_cachedir()
    savepath = f'{CACHEDIR}/{key}.pickle'
    if os.path.exists(savepath):
        return True
    else:
        return False


def store_cache_object(obj, key):
    ensure_cachedir()
    savepath = f'{CACHEDIR}/{key}.pickle'
    logging.debug(f"saving object key='{key}' ...")
    with open(savepath, 'wb') as fh:
        pickle.dump(obj, fh)
    logging.debug(f"saved object key='{key}' to {savepath}")
    
        
def get_cache_object(key):
    ensure_cachedir()
    savepath = f'{CACHEDIR}/{key}.pickle'
    logging.debug(f"loading key='{key}' from {savepath}")
    with open(savepath, 'rb') as fh:
        obj = pickle.load(fh)
    logging.debug(f'got object from {savepath}')
    return obj

