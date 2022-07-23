#!/usr/bin/env python
#
#   Vangone,A. and Bonvin,A.M. (2015) Contacts-based prediction of binding affinity in protein–protein 
#   complexes. eLife, 4, e07454. https://doi.org/10.7554/eLife.07454.
#
#   Requires prodigy (and its dep freesasa) in PATH, e.g.
#
#    conda create -n prodigy
#    conda activate prodigy
#    conda install -c conda-forge -c bioconda python=3 biopython json-c libxml2 libcxx compilers autoconf
#    git clone https://github.com/mittinatten/freesasa.git
#    cd freesasa/ ;  git submodule init ; git submodule update
#    autoreconf -i 
#    ./configure --prefix=$CONDA_PREFIX
#    make 
#    make install
#    git clone http://github.com/haddocking/prodigy
#    cd prodigy
#    pip install . 
#

import argparse
import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

def run_prodigy_on_dir(indir, outdir=None):
    '''
    runs prodigy on every .pdb file in directory. 
    optionally outputs .prd (prodigy) output files in different outdir. otherwise in indir. 
    
    usage: prodigy [-h] -q] [-V]
               [--selection A B [A,B C ...]]
               structf
  --distance-cutoff DISTANCE_CUTOFF
                        Distance cutoff to calculate ICs
  --acc-threshold ACC_THRESHOLD
                        Accessibility threshold for BSA analysis
  --temperature TEMPERATURE
                        Temperature (C) for Kd prediction
  --contact_list        Output a list of contacts
  --pymol_selection     Output a script to highlight the interface (pymol)        

  --selection A B [A,B C ...]
 
      By default, all intermolecular contacts are taken into consideration,
      a molecule being defined as an isolated group of amino acids sharing
      a common chain identifier. In specific cases, for example
      antibody-antigen complexes, some chains should be considered as a
      single molecule.
  
      Use the --selection option to provide collections of chains that should
      be considered for the calculation. Separate by a space the chains that
      are to be considered _different_ molecules. Use commas to include multiple
      chains as part of a single group:
  
      --selection A B => Contacts calculated (only) between chains A and B.
      --selection A,B C => Contacts calculated (only) between chains A and C; and B and C.
      --selection A B C => Contacts calculated (only) between chains A and B; B and C; and A and C.
          
    '''
    pdbfiles = glob.glob( f'{indir}/*.pdb')

    for fname in pdbfiles:        
        cmd = ['prodigy',
               fname
               ]
        try:
            run_command(cmd)
        except NonZeroReturnException as nzre:
            logging.error(f'problem with {infile}')
            logging.error(traceback.format_exc(None))
            raise  
    
    

if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
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
     
    parser.add_argument('-o','--outdir', 
                        metavar='outdir',
                        required=False,
                        type=str,
                        help='top output dir ')

    parser.add_argument('-i','--indir', 
                        metavar='indir',
                        required=True,
                        type=str,
                        help='input dir. one file per infile ')

    args= parser.parse_args()
    
    if args.debug:  
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
        

