#!/usr/bin/env python
#
# Various statistics, plotting utils
#
#
import logging
import pandas as pd
import numpy as np
import seaborn as sns
import math
import matplotlib.pyplot as plt
import seaborn as sns

def rank_standardize(series):
    '''
     takes standard Pandas column array/Series and returns a rank-standardized series in 
     the original order, such that the values represent the fraction of all values that 
     fall below that rank. 
        
    '''
    myser = series.copy()
    n_elem = len(myser)
    rnorm = (n_elem - myser.rank(ascending = False)) / n_elem
    return rnorm

#
#  Statistical functions...
#

def gini_coefficient_fast(X):
    """ 
        expects a CSR sparse matrix
        Sorting is O(n log n ) (here n is at most number of genes)
        loops over cells (m) instead of gene pairs. 
        Overall, at most O( m n log n)  but realistically, 
        density of 10% -> `effective n` is 0.1 * n
        
        only looks at nonzero elements
        Cells with no expression get a gini score of 0       
    """    
    # x = np.asarray(x)
    g = np.zeros(X.shape[0])    # ncells
    n = X.shape[1]          # ngenes
    for i in range(X.shape[0]): # loops for all cells
        # take the nonzero elements of the ith cell
        x = X[i,:]  
        x= x[:, x.indices].A.flatten()

        sorted_x = np.sort(x)   
        cumx = np.cumsum(sorted_x, dtype=float)

        if len(cumx) == 0 : # cell with zero expression - perfect equilibrium
            g[i] = 0
        else :
            g[i] =(n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return g


def sparse_pairwise_corr(A, B=None):
    """
    Compute pairwise correlation for sparse matrices. 
    Currently only implements pearson correlation.
    If A is N x P
       B is M x P
    Result is N+M x N+M symmetric matrix
    with off diagonal blocks as correlations between 
        elements in A with elements in B
    and main diagonal blocks as correlations between
        elements in A (or B) with elements in A (or B)
    """
    logging.debug(f'A.shape={A.shape} ')
    if B is None:
        logging.debug(f'B is none. copying A.')
        B = A.copy()

    n = A.shape[1]
    m = B.shape[1]
    logging.debug(f'B.shape={B.shape} ')

    assert n == m

    numer = np.dot(A,B.T).todense() 
    asum = A.sum(1)
    bsum = B.sum(1) 
    numer = n*numer - np.dot(asum,bsum.T) 

    sa =  np.sqrt(n*A.multiply(A).sum(1) - np.multiply(asum, asum))
    sb =  np.sqrt(n*B.multiply(B).sum(1) - np.multiply(bsum, bsum))

    denom = np.dot(sa, sb.T)
    return(np.asarray(numer/denom))


def pairwise_minmax_corr(X,chunksize = 5000 ):
    # X should be a cell x gene csr matrix 
    # be sure to onlycalculate the upper tri
    max_corr = np.ones(X.shape[0]) *-100
    min_corr = np.ones(X.shape[0]) * 100
    
    if chunksize == None:
        chunksize = min(X.shape[0] , 5000 ) 

    nchunks = int(np.ceil(X.shape[0] /  chunksize))

    logging.debug(f'nchunks={nchunks}, chunksize={chunksize} ')

    for i in range(nchunks ): 
        
        A = X[i*chunksize : (i+1)*chunksize,:] 
        for j in range(i,nchunks):
            B = X[j*chunksize : (j+1)*chunksize ,:] 
            logging.debug(f'working on: {i+1}/{j+1} of {nchunks}/{nchunks}' )
            current_corr = sparse_pairwise_corr(A,B )

            if i == j :
                # A ~= B distinct groups
                np.fill_diagonal(current_corr , np.nan)
                

            # np.argpartition(current_corr, n_neighbors  )
            cur_min = np.nanmin(current_corr,axis = 1)
            cur_max = np.nanmax(current_corr,axis = 1)
            
            # fmin/fmax ignores nan value where correlation is with itself. 
            max_corr[i*chunksize : (i+1)*chunksize] = np.fmax(cur_max ,max_corr[i*chunksize : (i+1)*chunksize] )
            min_corr[i*chunksize : (i+1)*chunksize] = np.fmin(cur_min ,min_corr[i*chunksize : (i+1)*chunksize] )

    return( min_corr, max_corr)


# EGAD functions compliments of Ben
def rank(data, nan_val=.5):
    """Rank normalize data
    
    Rank standardize inplace 
    Ignores Nans and replace with .5
    
    Does not return 
    Arguments:
        data {np.array} -- Array of data
    
    """
    finite = np.isfinite(data)
    ranks = bottleneck.rankdata(data[finite]).astype(data.dtype)

    ranks -= 1
    top = np.max(ranks)
    ranks /= top
    data[...] = nan_val
    data[np.where(finite)] = ranks

    return(data)


    



