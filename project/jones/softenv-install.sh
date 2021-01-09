#!/bin/bash -l
#  For jones project. pytorch
#
set -x

# Conda
conda create -y -n pytorch python=3.8
conda activate pytorch

# Conda available packages
conda install -c conda-forge -c bioconda snakemake pytorch pandas numpy scipy seaborn matplotlib plotly h5py ipython

