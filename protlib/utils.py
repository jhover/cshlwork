
import os
import logging
import sys
import traceback

from configparser import ConfigParser
import io
import pandas as pd


def read_identifier_file(filepath, flatten=True):
    """
    read and parse several formats for protein file
        
    1 per line
    multi per line: 
        if flatten=True just add to overall list. 
        otherwise each item in return list represents one line. multi-tiems in sub-list
    comma-separated
    ignore empty line
    ignore comment lines/extensions
    remove duplicates
    
    return list of items. 
    """
    idlist = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                idx = line.find('#')
                # get rid of comments
                if idx > -1:
                    line = line[:idx].strip()  
                if len(line) > 0:
                    if ',' in line:
                        line = line.replace(',',' ')
                    fields = line.split()
                    if flatten:
                        for f in fields:
                            f = f.strip()
                            if len(f) >0:
                                idlist.append(f)
                    else:
                        idlist.append(fields)
        idlist = list(set(idlist))
        logging.debug(f'got list with {len(idlist)} items.')
        return idlist
    except:
        return []    