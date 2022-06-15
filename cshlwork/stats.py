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

