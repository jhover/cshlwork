#!/usr/bin/env python
# 
#  Utilities for handling expression data. 
# 
#
import os, sys
gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import argparse
import datetime
import io
import logging
import traceback

from configparser import ConfigParser

import numpy as np
import pandas as pd
import h5py

# for clustering
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from dynamicTreeCut import cutreeHybrid


def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp

def load_coexp(hdf5file):
    f = h5py.File(hdf5file, 'r')
    dset = f['agg']
    npdset = np.array(dset)
    collist = list(np.array(f['col']).astype(str))
    rowlist = list(np.array(f['row']).astype(str))
    df = pd.DataFrame(npdset, columns = collist, index=rowlist )
    return df


def load_df(filepath):
    """
    Convenience method to load DF consistently accross modules. 
    """
    filepath = os.path.expanduser(filepath)
    df = pd.read_csv(filepath, sep='\t',index_col=0, keep_default_na=False, dtype =str, comment="#")
    df.fillna(value='', inplace=True)
    df = df.astype('str', copy=False)
    return df


def merge_write_df(newdf, filepath,  mode=0o644):
    """
    Reads existing, merges new, drops duplicates, writes to temp, renames temp. 
    """
    log = logging.getLogger('utils')
    log.debug(f'inbound new df:\n{newdf}')
    filepath = os.path.expanduser(filepath)
    if os.path.isfile(filepath):
        df = load_df(filepath)
        log.debug(f'read df:\n{df}')
        df = df.append(newdf, ignore_index=True)
        df.fillna(value='', inplace=True)
        df = df.astype('str', copy=False)
        log.debug(f'appended df:\n{df}')
    else:
        df = newdf
        df.fillna(value='', inplace=True)
        df = df.astype('str', copy=False)
    logging.debug(f"df length before dropping dupes is {len(df)}")
    df.drop_duplicates(inplace=True, ignore_index=True, keep='first')
    logging.debug(f"df length after dropping dupes is {len(df)}")
    df = df.reset_index(drop=True)
    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        df.to_csv(tfname, sep='\t')
        os.rename(tfname, filepath)
        os.chmod(filepath, mode)
        logging.info(f"wrote df to {filepath}")

    except Exception as ex:
        logging.error(traceback.format_exc(None))
        raise ex


def check_range(df, threshold):
    tdf=threshold_boolean(df, threshold)
    return sort_node_degree(tdf)


def threshold_boolean(df,  thresh):
    newdf = df.copy()
    newdf[newdf > thresh] = 1
    newdf[newdf <= thresh] = 0 
    return newdf

def threshold_value(df,  thresh, sense='above'):
    '''
    Only keep values relative to threshold. 
    'above' means higher, with others set to 0
    'below' means lower, with others set to 0
    
    '''
    thresh = float(thresh)
    newdf = df.copy()
    if sense == 'above':
        newdf[newdf > thresh] = newdf
        newdf[newdf <= thresh] = 0 
    if sense == 'below':
        newdf[newdf < thresh] = newdf
        newdf[newdf >= thresh] = 0 
    return newdf

def sort_node_degree(df):
    ts = df.sum(axis=1)
    ts.sort_values(inplace=True, ascending=False)
    return ts


def cluster_coexp(exphd5='~/data/cococonet/yeast_AggNet.hdf5', threshold=0.95, test=False, sense='above'):
    exphd5=os.path.expanduser(exphd5)
    logging.debug(f"clustering {exphd5} threshold={threshold} test={test} sense='{sense}'")
    edf = load_coexp(exphd5)
    logging.debug(edf)
    
    tdf = threshold_value(edf, threshold, sense=sense)
    np.fill_diagonal(tdf.values, 0)
    sdf = tdf.copy()
    sdf.drop(sdf.loc[sdf.sum(axis=1)==0].index, inplace=True)
    sdf.drop(columns=sdf.columns[sdf.sum()==0], inplace=True)
    logging.debug(sdf)
    
    # calc 1 - coexpression = distance. 
    dedf = 1 - sdf
    # pull a subset of the distance matrix
    if test:
        logging.info('test is true. subsetting to 200x200')
        dedf = dedf.iloc[0:200, 0:200]

    distances = pdist(dedf, metric="euclidean")
    logging.debug(distances)
    link = linkage(distances, "average")
    s = datetime.datetime.now()
    clusters = cutreeHybrid(link, distances)
    e = datetime.datetime.now()
    elapsed = e - s
    logging.info(f'Clustering took {elapsed.seconds} seconds.')
    logging.debug(len(clusters['labels']))
    logging.debug(clusters["labels"])
    larray = clusters["labels"]
    ldf = pd.DataFrame(larray, columns=['label'])
    ldf.index = dedf.index
    
    clusterdict = {}
    # extract gene names using assigned labels:
    outdf = pd.DataFrame(columns=['cluster','locus'])
    for lidx in range(1, larray.max() + 1 ):
        glist = list( dedf.loc[ldf['label'] == lidx].index )
        logging.debug(f'label={lidx} len={len(glist)} -> {glist} \n')
        tdf = pd.DataFrame(data=glist, columns=['locus'])
        tdf['cluster'] = lidx
        outdf = pd.concat([outdf, tdf], ignore_index=True  )

    return outdf





if __name__ == '__main__':

    gitpath = os.path.expanduser("~/git/alphaexp")
    sys.path.append(gitpath)

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
                            default='~/git/cshlwork/etc/exputils.conf',
                            help='Config file path [~/git/cshlwork/etc/exputils.conf]')

    parser.add_argument('-i', '--infile',
                        metavar='infile',
                        type=str,
                        default=None,
                        help='infile. ')

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.conffile is not None:
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile)) 
    else:
        pass
        
    cs = get_configstr(cp)
    logging.debug(f"got config: {cs}")
    logging.debug(f"args: {args}")    
    
    df = load_coexp(args.infile)
    print(df)
    