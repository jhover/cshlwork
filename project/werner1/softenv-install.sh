#!/bin/bash -l
#  For werner1 project.   labelled gatk. 
#
#  samtools-1.11
#  bedtools-2.29.2
#  gatk-4.1.9.0
#  igvtools-2.8.9
#  STAR 2.7.2a  (tied to genome dir version)
#  genome 2.4.2a?
set -x

# Conda
conda create -y -n gatk python=3.8
conda activate gatk

# Conda available packages
conda install -y pandas numpy scipy seaborn matplotlib plotly h5py 
conda install -c conda-forge -c bioconda snakemake

# Tarball installs
#  Samtools 1.11
#  samtools 1.9 for compatiblity?
#
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar -xvjf samtools-1.11.tar.bz2
cd samtools-1.11
./configure --prefix=$CONDA_PREFIX
make
make install
cd ..
rm -rf samtools-1.11.tar.bz2 samtools-1.11

#  bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
cp bin/* $CONDA_PREFIX/bin/
cd ..
rm -rf bedtools2 bedtools-2.29.1.tar.gz

#  gatk-4.1.9.0 
# ? 4.1.4.1 for compatibility?
wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
unzip gatk-4.1.9.0.zip
mv gatk-4.1.9.0 $CONDA_PREFIX/
cd $CONDA_PREFIX/bin
ln -s ../gatk-4.1.9.0/gatk ./
cd -
rm gatk-4.1.9.0.zip

# igvtools-2.8.9
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.13.zip
mv IGV_2.8.13 $CONDA_PREFIX/
cd $CONDA_PREFIX/bin/
ln -s  $CONDA_PREFIX/IGV_2.8.13/igvtools ./
cd -
rm IGV_2.8.13.zip

# STAR
# may need to brew install gcc llvm libomp on mac. 
wget https://github.com/alexdobin/STAR/archive/2.7.2a.tar.gz
tar -xvzf 2.7.2a.tar.gz
cd STAR-2.7.2a/source
make STARforMacStatic CXX=/usr/local/Cellar/gcc/10.2.0/bin/g++-10 
cp STAR $CONDA_PREFIX/bin/
cd ..
rm -rf STAR-2.7.2a 2.7.2a.tar.gz




	




