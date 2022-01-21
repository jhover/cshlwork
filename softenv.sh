#!/usr/bin/env -S bash -l
set -e
# install conda
if ! command -v conda&> /dev/null; then
	echo "installing miniconda..."
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
	bash Miniconda3-py39_4.9.2-Linux-x86_64.sh -b
	rm -f Miniconda3-py39_4.9.2-Linux-x86_64.sh
	~/miniconda3/bin/conda init 
	echo "miniconda installed. restart terminal."
	exit 0
else
	echo "miniconda installed already."
fi
# be sure to restart terminal to allow conda to start

conda create -y -n cshlwork python=3.9
sleep 5
conda activate cshlwork

# make sure we have conda-forge
conda config --add channels conda-forge

# update conda if necessary
conda update -y -n base -c defaults conda
conda install -y pandas ipython
conda install -y -c bioconda hmmer
