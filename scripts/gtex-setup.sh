# werner/gatk/samtools/bedtools

# GATK=/home/werner/tools/gatk4.1.3.0/gatk/gatk
#    https://github.com/broadinstitute/gatk/releases
#   

# igvtools=/home/werner/tools/IGVTools/igvtools.jar 
#    http://www.broadinstitute.org/software/igv/igvtools_commandline
#    unpack at ~/.conda/envs/gatk/
#    cd ~/.conda/envs/gatk/bin ; ln -s ~/.conda/envs/gatk/igvtools-x.y.z/bin/igvtools
#

# samtools=/home/werner/tools/samtools-1.9/samtools 1.9
#     http://www.htslib.org/
#      https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
#
# bedtools
#  wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
#  make
#  cp bin/* ~/. 
#

module load Anaconda3
module load Java

conda create -n gatk
source activate gatk

conda install -c bioconda gatk      # 3.8
#conda install -c bioconda igvtools  # 2.5.3
#conda install -c bioconda bedtools  # 2.29.2
#conda install -c bioconda samtools  # 1.7

# bedtools
# samtools
#
