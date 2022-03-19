#!/usr/bin/env python
# takes in genome.fa and annotation.gtf
#   laptop times:
#   read in mouse genome (2.5GB):      17s 
#   process triplets in chrY (91Mbp)   55s
#   
#{'gff_id': {('chr1',): 116585,
#            ('chr10',): 86823,
# .
# .
#            ('chrY',): 13694},
# 'gff_source': {('ENSEMBL',): 139569, 
#                ('HAVANA',): 1729538},
# 'gff_source_type': {('ENSEMBL', 'CDS'): 45873,
#                     ('ENSEMBL', 'UTR'): 10566,
#                     ('ENSEMBL', 'exon'): 55778,
#                     ('ENSEMBL', 'gene'): 6373,
#                     ('ENSEMBL', 'start_codon'): 4665,
#                     ('ENSEMBL', 'stop_codon'): 4656,
#                     ('ENSEMBL', 'transcript'): 11658,
#                     ('HAVANA', 'CDS'): 481982,
#                     ('HAVANA', 'Selenocysteine'): 65,
#                     ('HAVANA', 'UTR'): 175369,
#                     ('HAVANA', 'exon'): 786174,
#                     ('HAVANA', 'gene'): 48986,
#                     ('HAVANA', 'start_codon'): 55246,
#                     ('HAVANA', 'stop_codon'): 50999,
#                     ('HAVANA', 'transcript'): 130717},
# 'gff_type': {('CDS',): 527855,
#              ('Selenocysteine',): 65,
#              ('UTR',): 185935,
#              ('exon',): 841952,
#              ('gene',): 55359,
#              ('start_codon',): 59911,
#              ('stop_codon',): 55655,
#              ('transcript',): 142375}}
#<class 'map'>

#  Useful ipython setup lines...
#  genome_file = os.path.expanduser('~/data/genomes/Mus_musculus/genome.fa')
#  annot_file = os.path.expanduser('~/data/genomes/Mus_musculus/annotation.gtf')
#  logging.getLogger().setLevel(logging.DEBUG)
#  genome_dict, annot_dict, gene_idx = scan_genome(genome_file, annot_file)

import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from genomelib.genome import read_genome, read_annot, scan_genome

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
    
    parser.add_argument('-g','--genome', 
                        metavar='genome',
                        required=False,  
                        type=str,
                        default=os.path.expanduser('~/data/genomes/Mus_musculus/genome.fa'), 
                        help='Genome .fasta file') 

    parser.add_argument('-a','--annot', 
                        metavar='annot',
                        required=False,  
                        type=str,
                        default=os.path.expanduser('~/data/genomes/Mus_musculus/annotation.gtf'), 
                        help='Annotation .gtf file') 

    parser.add_argument('-o','--outfile', 
                        metavar='outfile',
                        required=False, 
                        type=str, 
                        help='out file.')    

    parser.add_argument('-s','--species', 
                        metavar='species',
                        required=False, 
                        type=str,
                        default='Mus_musculus', 
                        help='Linnean species w/ underscore.')
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
    
    genome_dict, annot_dict, chr_idx = scan_genome(args.genome, args.annot, args.species)
    
    
    