#!/usr/bin/env python
#
#  Splits of first N bases of infile to bcfile, remaining sequence written to seqfile. 
#  Used to fix input files for STARsolo processing
#
#  Usage example:  
#     fastq_split.py -d SRR15686117_1.fastq SRR15686117_1_bc.fastq SRR15686117_1_seq.fastq
#
 
import argparse
import logging
import os
import sys

def process_files(infile, bcfile, seqfile, n_bases):
    infile  = os.path.abspath(infile)
    bcfile  = os.path.abspath(bcfile)
    seqfile = os.path.abspath(seqfile)
    logging.debug(f'infile={infile} bcfile={bcfile} seqfile={seqfile} n_bases={n_bases}')
    with open(infile) as infh:
        with open(bcfile, 'w') as bcfh:
            with open(seqfile, 'w') as sfh:
                for line in infh:
                    logging.debug(line)
                    if line.startswith('@') or line.startswith('+'):
                        bcfh.write(line)
                        sfh.write(line)                  
                        logging.debug('wrote header line without change...')
                    else:
                        b = line[:n_bases]
                        s = line[n_bases:]
                        bcfh.write(f'{b}\n')
                        sfh.write(s)
                        logging.debug(f'wrote {len(b)} bases to barcode file. {len(s) -1} bases to sequence file.')
                            

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

    parser.add_argument('-b','--bases', 
                        metavar='bases',
                        required=False,  
                        type=int,
                        default=28,
                        help='Number of bases to strip to bcfile')

    parser.add_argument('infile' ,
                        metavar='infile', 
                        type=str,
                        default=None, 
                        help='Fastq to split')    

    parser.add_argument('bcfile', 
                        metavar='bcfile', 
                        type=str, 
                        help='Fastq file for barcodes.')    

    parser.add_argument('seqfile', 
                        metavar='seqfile', 
                        type=str, 
                        help='Fastq file for remaining sequence.')    

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
        
    process_files(os.path.expanduser(args.infile), 
                  os.path.expanduser(args.bcfile), 
                  os.path.expanduser(args.seqfile),
                  args.bases  )
    