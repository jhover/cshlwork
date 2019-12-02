#!/usr/bin/env python
#
#  Peform processing steps, output prediction file(s)
#  One file per species. 
#
#  cafa [OPTIONS}  <targetfile>
#      -c <configfile>
#      -o <predictionfile> [gillislab_1_
#      -profile <profilename>
#
#     
#
#   HOVER GILLIS_LAB
#   MODEL 1
#   KEYWORDS 
#   [ACCURACY] 
#   T100900000001  
#

#importing libraries
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import cafalib

#%matplotlib inline