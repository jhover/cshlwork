import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import pprint
from collections import defaultdict

from Bio import SeqIO
# https://github.com/chapmanb/bcbb/tree/master/gff
from BCBio import GFF
from BCBio.GFF import *

from genomelib.biomart import get_mapping_by_gene


#   http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/
#   http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/ 
#            genome           annot
#   chr17   95294699         95185749
#   
#


def read_genome(filepath, chr=None):
    logging.info(f'opening file {filepath}')
    seq_recs = {}
    for seq_record in SeqIO.parse(filepath, 'fasta'):
        logging.debug(f'reading {seq_record.id} length={len(seq_record)}...')
        seq_recs[seq_record.id] = seq_record
    
    logging.info(f'collected {len(seq_recs)} chromosomes/contigs from {filepath}')
    return seq_recs
        
    
def read_annot(filepath, chr=None, gff_type=None):
    #for f in rec.features:
    #    features.append(f)
    #print(rec.features[0])
    #    print(rec)   
    #limits = examiner.available_limits(f)
    #pprint.pprint(limits)
    #logging.debug('building map...')
    #map = examiner.parent_child_map(f)

    #strmap = str(map)
    #logging.debug(f'done with map, print: {strmap}')
    #pprint.pprint(map)
    
    #limit_info = dict(gff_id=["chr1"], gff_source=['ENSEMBL'])
    #limit_info = dict(gff_type=['exon'])

    #for rec in GFF.parse(f, target_lines=1000):
    #recs = []
    #features = []
    #for rec in GFF.parse(fh, limit_info=limit_info):
    #for rec in GFF.parse(fh, limit_info=limit_info):
    
    logging.info(f'opening file {filepath}')
    try:
        fh = open(filepath)    
        seq_recs = {}
        for rec in GFF.parse(fh):
            logging.debug(f'handling record, id={rec.id}')
            #pprint.pprint(rec.features)
            seq_recs[rec.id] = rec
            logging.debug(f'rec {rec.id} has {len(rec.features)} features.')

    except Exception as ex:
        logging.error(traceback.format_exc(None))
    finally:
        fh.close()
    #return seq_recs mapped by chromosome id. 
    logging.info(f'collected annotations for {len(seq_recs)} chromosomes/contigs from {filepath}')
    return seq_recs


def load_combined(genome_file, annot_file):
    logging.info(f'opening file {genome_file}')
    sfh = open(genome_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(sfh, "fasta"))
    sfh.close()
    
    seq_recs = []
    afh = open(annot_file)
    for rec in GFF.parse(afh, base_dict=seq_dict):
        seq_recs.append(rec)
   
    afh.close()
    logging.info(f'collected {len(seq_recs)} chromosomes/contigs from {genome_file} and {annot_file}')
    return seq_recs


def index_annotations(annot_dict, tg_map ):
    '''
        @arg annot_dict        Dictionary of annotations from read_annot
        @arg transcript-id -> gene_id map    Dict of transcript_ids to gene_ids
    
        @return
    ''' 
    # take dict of annotation SeqRecords and index by chromosome.
    # within chromosomes index by Gene (MUSG......XXXXX)
    # each gene has a .sub_features list of transcripts (splice variants)
    #
    #  srec.qualifiers['gene_id'] identifier, stripped of version number. 
    #
    #  gene_idx =  { 'chr1' : { 'ENSMUSG00000104991' : 
    #                                  [ SeqFeature; type=gene .sub_features= 
    #                                            [ SeqFeature; type = inferred_parent, .sub_features= [ List of SeqFeatures..] 
    #                                             SeqFeature; type= inferred_parent, ... ]  
    #
    #   Example:   Rfc3 ENSMUSG00000033970     7 splice variants. 
    #                         ENSMUST00000038131.10

    #                         ENSMUST00000140067.8
    #                         ENSMUST00000132709.2
    #                            ...

    chr_dict = {}
    for chr in annot_dict.keys():
        logging.info(f'indexing {chr}...')
        gene_dict = {}
        num_transcripts = 0
        num_genes = 0
        chr_sr = annot_dict[chr]
        logging.debug(f'handling {len(chr_sr.features)} ')
        logging.debug(f'handling genes only... ')
        for sf in chr_sr.features:
            #logging.debug(f' {transcript_sf} of type {type(transcript_sf)}')
            #  type <class 'Bio.SeqFeature.SeqFeature'>
            try:
            # handle GENES
                if sf.type == 'gene':
                    num_genes +=1
                    gene_id = sf.qualifiers['gene_id'][0].split('.')[0]
                    gene_type = sf.qualifiers['gene_type']
                    gene_name = sf.qualifiers['gene_name']
                    sf.id = gene_id
                    # id empty
                    # sub_features empty
                    logging.debug(f'gene_id is = {gene_id} getting dict entry... ')
                    gene_dict[gene_id] = sf
            except Exception as ex:
                pass 

        # Now we assume all genes on chromosome are present in gene_dict
        # run through all seqfeatures again, only dealing with transcripts. 
        logging.debug(f'added {len(gene_dict.keys())} genes, handling transcripts... ')
        for sf in chr_sr.features:
            try:
                # handle TRANSCRIPTS
                if sf.type == 'inferred_parent':
                    trans_id = sf.id.split('.')[0]  # id is transcript id
                    subflist = sf.sub_features  # subfeatures   type=transcript, exon, CDS, etc.
                    gene_id = tg_map[trans_id]
                    gsf = gene_dict[gene_id]
                    gsf.sub_features.append(sf)                    
                    logging.debug(f'adding {trans_id} to {gene_id}')
                    num_transcripts += 1
            except Exception as ex:
                pass     
        
        logging.debug(f'handled {num_transcripts} transcripts, {num_genes} genes for chromosome/contig {chr}')
        logging.debug(f'final gene_dict has {len(gene_dict)} genes')
        chr_dict[chr] = gene_dict
    logging.debug(f'created chr_dict len={len(chr_dict)}')
    return chr_dict

def reclist_index_transcripts(reclist, tg_map ):
    '''
        @arg annot_dict        Dictionary of annotations from read_annot
        @arg transcript-id -> gene_id map    Dict of transcript_ids to gene_ids
    
        @return
    ''' 
    # take dict of annotation SeqRecords and index by chromosome.
    # within chromosomes index by Gene (MUSG......XXXXX)
    # each gene has a .sub_features list of transcripts (splice variants)
    #
    #  srec.qualifiers['gene_id'] identifier, stripped of version number. 
    #
    #  gene_idx =  { 'chr1' : { 'ENSMUSG00000104991' : 
    #                                  [ SeqFeature; type=gene .sub_features= 
    #                                            [ SeqFeature; type = inferred_parent, .sub_features= [ List of SeqFeatures..] 
    #                                             SeqFeature; type= inferred_parent, ... ]  
    #
    #   Example:   Rfc3 ENSMUSG00000033970     7 splice variants. 
    #                         ENSMUST00000038131.10

    #                         ENSMUST00000140067.8
    #                         ENSMUST00000132709.2
    #                            ...
    for rec in reclist:
        logging.info(f'indexing {rec.id}...')
        gene_dict = {}
        num_transcripts = 0
        num_genes = 0
        logging.debug(f'handling {len(rec.features)} features..')
        logging.debug(f'handling genes only... ')
        for sf in rec.features:
            try:
            # handle GENES
                if sf.type == 'gene':
                    num_genes +=1
                    gene_id = sf.qualifiers['gene_id'][0].split('.')[0]
                    #gene_type = sf.qualifiers['gene_type']
                    gene_name = sf.qualifiers['gene_name']
                    sf.id = gene_id
                    # id empty
                    # sub_features empty
                    logging.debug(f'gene_id is = {gene_id} setting dict entry... ')
                    gene_dict[gene_id] = sf
            except Exception as ex:
                pass 
        logging.debug(f'added {len(gene_dict.keys())} genes, handling transcripts... ')
        # Now we assume all genes on chromosome are present in gene_dict
        # run through all seqfeatures again, only dealing with transcripts. 

        for sf in rec.features:
            try:
                # handle TRANSCRIPTS
                if sf.type == 'inferred_parent':
                    trans_id = sf.id.split('.')[0]  # id is transcript id
                    subflist = sf.sub_features  # subfeatures   type=transcript, exon, CDS, etc.
                    gene_id = tg_map[trans_id]
                    gsf = gene_dict[gene_id]
                    logging.debug(f'handling gene {gsf}')
                    gsf.sub_features.append(sf)                    
                    logging.debug(f'adding {trans_id} to {gene_id}')
                    num_transcripts += 1
            except Exception as ex:
                pass     
        
        logging.debug(f'handled {num_transcripts} transcripts, {num_genes} genes for chromosome/contig {rec}')
        logging.debug(f'final gene_dict has {len(gene_dict)} genes')
        
    #logging.debug(f'created chr_dict len={len(rec)}')
    return reclist



def add_transcripts(annot_dict, tg_map ):
    '''
        @arg annot_dict        Dictionary of annotations from read_annot
        @arg transcript-id -> gene_id map    Dict of transcript_ids to gene_ids
    
        @return
    ''' 
    for chr in annot_dict.keys():
        logging.info(f'handling {chr}...')
        chr_sr = annot_dict[chr]
        # build gene object index
        gene_idx = {}
        num_genes = 0
        num_transcripts = 0
        for sf in chr_sr.features:
            #logging.debug(f' {transcript_sf} of type {type(transcript_sf)}')
            #  type <class 'Bio.SeqFeature.SeqFeature'>
            try:
            # handle GENES
                if sf.type == 'gene':
                    num_genes +=1
                    gene_id = sf.qualifiers['gene_id'][0].split('.')[0]
                    gene_type = sf.qualifiers['gene_type']
                    gene_name = sf.qualifiers['gene_name']
                    sf.id = gene_id
                    # id empty
                    # sub_features empty
                    logging.debug(f'gene_id is = {gene_id} getting dict entry... ')
                    gene_idx[gene_id] = sf
                    num_genes +=1
            except Exception as ex:
                pass 

        # Now we assume all genes on chromosome are present in gene_dict
        # run through all seqfeatures again, only dealing with transcripts. 
        logging.debug(f'added {len(gene_idx.keys())} genes, handling transcripts... ')
        for sf in chr_sr.features:
            try:
                # handle TRANSCRIPTS
                if sf.type == 'inferred_parent':
                    trans_id = sf.id.split('.')[0]  # id is transcript id
                    subflist = sf.sub_features  # subfeatures   type=transcript, exon, CDS, etc.
                    gene_id = tg_map[trans_id]
                    gsf = gene_idx[gene_id]
                    gsf.sub_features.append(sf)                    
                    logging.debug(f'adding {trans_id} to {gene_id}')
                    num_transcripts +=1
            except Exception as ex:
                pass     
        
        logging.debug(f'handled {num_transcripts} transcripts, {num_genes} genes for chromosome/contig {chr}')
        logging.debug(f'done handling annot. returning ')
    return annot_dict

def merge_sequence(genome_dict, annot_dict):
    '''
    Replaces UnknownSeq in annotation dict with actual genomic sequence in genome_dict.  
    '''
    for chr in annot_dict.keys():
        sr = annot_dict[chr]
        gseq = genome_dict[chr].seq
        logging.debug(f'replacing unknown seq of len {len(sr.seq)} with genome seq of len {len(gseq)}')
        #sr.seq = gseq
    

def assess_dict_index(chr_index):
    one = 0
    two = 0
    three = 0
    four = 0
    five = 0
    
    
    chrs = chr_index.keys()
    logging.debug(f'index has {len(chrs)} chromosomes.')    
    for chr in chrs:
        genes = chr_index[chr]
        logging.debug(f'chromsome {chr} has {len(genes)} genes.')
        for g in genes.keys():
            tlist = genes[g]
            logging.debug(f'gene {g} has {len(tlist)} transcripts.')
            if len(tlist) ==1:
                one += 1
            if len(tlist) ==2:
                two += 1            
            if len(tlist) ==3:
                three += 1
            if len(tlist) ==4:
                four += 1  
            if len(tlist) ==5:
                five += 1
         
    logging.debug(f'one={one} two={two} three={three} four={four} five={five}')


def get_sequence_stats(seq):
    '''
        overall length
        first non-N position
        last non-N position
        final all-N length. 
    '''
    pass
    




def scan_genome(genome_file, annot_file, species='Mus_musculus'):
    #
    #  returns 
    #      genome dict of SeqRecords, one per chromosome/contig
    #      annot_dict of annotated SeqRecord objects, one per transcript
    #      gene_idx dict  indexed map of chr -> gene_id -> list of transcript SeqRecord objects
    # 
    #(recs, features) = read_annot(args.annot, 'chr1')
    logging.info(f'getting genome from {genome_file} ...')    
    genome_dict = read_genome(genome_file)
    logging.info(f'getting annotations from {annot_file} ...')    
    annot_dict = read_annot(annot_file)    
    logging.info(f'getting transcript gene mappings from ensembl ...')  
    transmap = get_mapping_by_gene(species=species)
    logging.info(f'building gene->transcript index...')    
    index_annotations(annot_dict, transmap)
    logging.info(f'built index of {len(chr_idx.keys())} chromsomes/contigs.')
    full_dict =merge_sequence(genome_dict, annot_dict)
    
    #logging.debug(chr_idx)
    return(genome_dict, annot_dict)


def save_gff3(annot_dict):
    srlist = []
    for k in annot_dict.keys():
        srlist.append(annot_dict[k])
        
    
    
    
    
