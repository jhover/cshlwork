#!/usr/bin/env python
#
import mkl
mkl.set_num_threads(16)
​
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scanpy as sc
import bottleneck
from scipy import stats
import gc
​
import matplotlib.pyplot as plt
import seaborn as sns
​
sns.set(style='white', font_scale=1.25)
plt.rc("axes.spines", top=False, right=False)
plt.rc('xtick', bottom=True)
plt.rc('ytick', left=True)
​
from itertools import combinations
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
​
import networkx as nx
​
import sys
sys.path.append('../scripts/')
sys.path.append('/home/bharris/Correlation_Coexpression/scripts/')
sys.path.append('/home/bharris/vshape/scripts/')
​
from rank import rank
from processify import processify​
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description='Decompose it!')
    parser.add_argument('--dataset', type=str, help='Which dataset to run on')
    return parser.parse_args()
​
​
def rank(data):
    """Rank normalize data
    
    Rank standardize data to make nonparametric
    
    Arguments:
        data {np.array} -- 2-D coexpression network
    
    Returns:
        np.array -- Rank normalized between 0 and 1 array
    """
    orig_shape = data.shape
    data = bottleneck.nanrankdata(data) - 1
    return (data / np.sum(~np.isnan(data))).reshape(orig_shape)
​
def create_nw(data, replace_nans):
    nw = np.corrcoef(data)
    np.fill_diagonal(data, 1)
    nw = rank(nw)
    if replace_nans:
        nw[np.isnan(nw)] = bottleneck.nanmean(nw)
    return nw
​
​
dataset_dict = pd.read_csv(
    '/home/bharris/biccn_paper/data/dataset_dict_biccn_sets_7.csv',
    index_col=0).to_dict()
​
genes = np.genfromtxt(
    '/home/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv',
    dtype=str)
​
args = parse_args()
dataset = args.dataset
​
logging.info(dataset)
andata = sc.read_h5ad(dataset_dict[dataset]['andata'])
sc.pp.normalize_total(andata, target_sum=1e6)
andata2 = andata[:, genes]
del andata
gc.collect()
​
metacell_assignment = pd.read_csv(
    f'/home/bharris/metacells_9_datasets/{dataset}/{dataset}meta_cell_assignment.csv',
    index_col=0)
andata = andata2[metacell_assignment.index]
del andata2
gc.collect()
​
expression = andata.to_df()
agg_nw = np.zeros([genes.shape[0], genes.shape[0]])
metacell_values = np.unique(metacell_assignment['x'].values)
logging.info(np.max(metacell_values))
for metacell in metacell_values:
    logging.info(metacell)
    #Generate Mask for slicing
    
    mask = metacell_assignment['x'] == metacell
    
    if mask.sum() < 20:
        logging.info(f'{metacell} too small')
        del mask
        gc.collect()
        continue
    
    #Expression is DataFrame of cells x genes 
    data = expression[mask].values.T
    nw = create_nw(data, True)
    agg_nw += nw
    del nw, mask , data
    gc.collect()
​
dataset_nw = pd.DataFrame(rank(agg_nw), index=genes, columns=genes)
dataset_nw.to_hdf(
    f'/home/bharris/biccn_paper/data/networks/metacells/metacell_agg_nw_{dataset}.hdf5',
    'nw')