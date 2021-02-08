#
#
#
#
conda create -n fastqtl

# For chapter 2  in eQTL analysis 
conda install -c conda-forge -c bioconda r-base r-genomictools r-tidyverse r-reticulate

# Python/Pandas based...
conda install -c bioconda pyranges

#
#  gr = pr.read_gtf("ensembl.gtf")
#  # as DataFrame
#df = gr.df
#
#

cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh


# Install RMATH
cd ~/CONDA_PREFIX/ ; mkdir src ; cd src
wget http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz
tar xzvf R-3.2.0.tar.gz
cd R-3.2.0
./configure --prefix=$CONDA_PREFIX --with-x=no
cd src/nmath/standalone


sudo yum install boost-devel gsl-devel

# install fastqtl
wget http://fastqtl.sourceforge.net/files/FastQTL-2.184.linux.tgz
tar -xvzf FastQTL-2.184.linux.tgz
cd FastQTL-2.184
# edit Makefile to set RMATH to dir established above...
make cleanall && make
mkdir $CONDA_PREFIX/bin
cp bin/* $CONDA_PREFIX/bin/


#
# Only needed when RPM not installable...
#

#install gsl
#wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
#./configure --prefix=$CONDA_PREFIX
# make 
# make install


#install boost
#https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz
#./bootstrap.sh
#b2


# install zlib
wget http://www.zlib.net/zlib-1.2.11.tar.gz
 ./configure --prefix=$CONDA_PREFIX
make
make install

#rmath?
# conda install -c r r-essentials         
#	(1) Download R source code:             wget http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz
#        (2) Unzip R source code:                        tar xzvf R-3.2.0.tar.gz
#        (3) Go to R source code folder:         cd R-3.2.0
#        (4) Configure Makefile:                         ./configure
#        (5) Go to R math library folder:        cd src/nmath/standalone
#        (6) Compile the code:                           make
#        (7) Go 2 folder backward:                       cd ../..
#        (8) Save the current path:                      RMATH=$(pwd)
wget http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz


