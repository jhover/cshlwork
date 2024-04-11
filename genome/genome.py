#!/usr/bin/env python
# 
#  Module for generic genomic processing
#  

import logging
import operator
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.caching import *
from cshlwork.utils import *


import pandas as pd

import pprint
from collections import defaultdict

from Bio import SeqIO
# https://github.com/chapmanb/bcbb/tree/master/gff
from BCBio import GFF
from BCBio.GFF import *

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation



class NonZeroReturnException(Exception):
    """
    Thrown when a subprocess command has non-zero return code. 
    """
    

#   http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/
#   http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/ 
#            genome           annot
#   chr17   95294699         95185749

UTR_TYPES = ['five_prime_utr', 'three_prime_utr']


# generic genomic processing. stop assuming particular structure (vs. cellreader). 
#
#   seqrecords have base sequences. 
#   annotrecords have structure. 
#
#   SeqRecord  seq id
#            .features ->   [ SeqFeature
#                                        .subfeatures
#
#

def metric_gc_content(sequence):
    """Returns GC content"""
    return round((sequence.count("G") + sequence.count("C"))/(len(sequence)),3)

def load_seqrecords(genome_file):
    logging.info(f'opening file {genome_file}')
    with open(genome_file) as sfh:
        seq_dict = SeqIO.to_dict(SeqIO.parse(sfh, "fasta"))
    logging.debug(f'handled genome file {genome_file} w/ {len(seq_dict)} records.')
    return seq_dict


def get_parent_dict(annot_file):
    logging.info(f'opening file {annot_file}')
    examiner = GFFExaminer()
    with open(annot_file) as afh:
        pcmap = examiner.parent_child_map(afh)
    return pcmap

def get_limits_dict(annot_file):
    logging.info(f'opening file {annot_file}')
    examiner = GFF.GFFExaminer()
    with open(annot_file) as afh:
        limits_dict = examiner.available_limits(afh)
    return limits_dict

def load_annotrecords(annot_file):
    '''
    limits_info = dict( gff_type=['mRNA','CDS']) ?
    
    '''
    logging.info(f'opening file {annot_file}')
    arecord_list = list( GFF.parse(annot_file))
    logging.info(f'handled annotation file {annot_file} w/ {len(arecord_list)} sets of features.')
    return arecord_list


def get_feature_type_list(sub_features_list):
    tset = set()
    for sf in sub_features_list:
        tset.add(sf.type)
    return list(tset)


def get_genes_generic_gtf(arecord_list, 
              srecord_dict=None, 
              feature= 'gene',  # top level object below record 
              sub_feature='transcript',  # mid-level to handle.        
              seq_feature='exon', #  lowest-level sequences to be extracted.
              ):
    '''
    arecord_list:  list of annotation/feature object trees. 
    srecord_dict:  dictionary of seqrecords that are parent seqeunce of annotation records. 
       
    Generic GTF 3 LAYERS 
    Seqrecord(seq, id).features
        SeqFeature(type=gene ).sub_features
            SeqFeature(type=mRNA| type=inferred_parent).sub_features
                SeqFeature( type=exon )  | SeqFeature( type=CDS )

    '''
    transcript_list = []
    # annotation records
    arecs_handled = 0
    
    for arec in arecord_list:
        # arec is list of features on a sinlge parent SeqRecord
        rec_id = arec.id  # id of parent SeqRecord object. 
        if srecord_dict is not None:
            parent_record = srecord_dict[rec_id]
        else:
            parent_record = arec
        logging.debug(f'handling annotation record {rec_id}')
        feats_handled = 0
        sub_feats_handled = 0
        for feat in arec.features:        
            if feat.type == feature:
                logging.debug(f'handling feature {feat.id} type={feat.type}')
                feat_id = feat.id.replace('"', ':')
                feat_strand = feat.location.strand
                
                logging.debug(f'handling sub_features')                
                
                for sf in feat.sub_features:
                    if sf.type == sub_feature:
                        seqflist = []
                        for ssf in sf.sub_features:
                            if ssf.type == seq_feature:
                                seqflist.append(ssf)
                        if (len(seqflist) > 0):
                            # make new SeqFeature
                            mloc = get_merged_locations(seqflist)
                            exon1 = seqflist[0]
                            seqf = SeqFeature(location=mloc)
                            seqf_length = len(seqf)
                            if srecord_dict is not None:
                                seqf = seqf.extract(parent_record)
                            try:
                                gene_name = str( sf.qualifiers['gene_name'][0] )
                                gene_name = gene_name.upper()
                            except KeyError:
                                gene_name = 'undefined'
                            try:
                                transcript_id = str( sf.qualifiers['transcript_id'][0] )
                                transcript_id = transcript_id.replace('"', ':')
                            except KeyError:
                                transcript_id = 'undefined'                            
                            
                            seqf.type = sub_feature
                            seqf.id = f'{gene_name}:{transcript_id}'
                            seqf.name = f'{parent_record.id}:{feat_id}'
                            seqf.qualifiers = exon1.qualifiers
                            seqf.description = f"{parent_record.id}:{feat_id}:gene={gene_name}:type={seqf.type}:len={seqf_length}:n_{seq_feature}={len(seqflist)}:strand={feat_strand}"
                            if srecord_dict is not None:
                                logging.debug(f'made merged, extracted SeqRecord: {seqf}')
                            else:
                                logging.debug(f'made un-merged SeqFeature: {seqf}')
                            transcript_list.append(seqf)
                        else:
                            logging.debug(f'type {subfeat.type} not in {sub_features}')
                            
                    sub_feats_handled +=1
            else:
                logging.warning(f'type {feat.type} !=  {feature}')
            feats_handled += 1
        logging.debug(f'{feats_handled} for {rec_id}')
    logging.info(f'made list of {len(transcript_list)} objs')    
    return transcript_list



def get_genes_refseq_gff(arecord_list, 
              srecord_dict=None, 
              feature= 'gene',  # top level object below record 
              sub_feature='mRNA',  # mid-level to handle.        
              seq_feature='CDS', #  lowest-level sequences to be extracted.
              ):
    '''
    arecord_list:  list of annotation/feature object trees. 
    srecord_dict:  dictionary of seqrecords that are parent seqeunce of annotation records. 
       
    Refseq GFF   4 LAYERS
    Seqrecord(seq, id).features
        SeqFeature(type= gene | biological_region, pseudogene, enhancer |protein_binding_site ).sub_features   
            SeqFeature(type= mRNA | transcript | lnc_RNA ).sub_features
                SeqFeature( type=exon ) | SeqFeature( type=CDS )

    '''
    transcript_list = []
    # annotation records
    arecs_handled = 0
    
    for arec in arecord_list:
        # arec is list of features on a sinlge parent SeqRecord
        rec_id = arec.id  # id of parent SeqRecord object. 
        if srecord_dict is not None:
            parent_record = srecord_dict[rec_id]
        else:
            parent_record = arec
        logging.debug(f'handling annotation record {rec_id}')
        feats_handled = 0
        sub_feats_handled = 0
        for feat in arec.features:        
            if feat.type == feature:
                logging.debug(f'handling feature {feat.id} type={feat.type}')
                feat_id = feat.id.replace('"', ':')
                feat_strand = feat.location.strand
                
                logging.debug(f'handling sub_features')                
                
                for sf in feat.sub_features:
                    if sf.type == sub_feature:
                        seqflist = []
                        for ssf in sf.sub_features:
                            if ssf.type == seq_feature:
                                seqflist.append(ssf)
                        if (len(seqflist) > 0):
                            # make new SeqFeature
                            mloc = get_merged_locations(seqflist)
                            exon1 = seqflist[0]
                            seqf = SeqFeature(location=mloc)
                            seqf_length = len(seqf)
                            if srecord_dict is not None:
                                seqf = seqf.extract(parent_record)
                            try:
                                gene_name = str( feat.qualifiers['gene'][0] )
                                gene_name = gene_name.upper()
                            except KeyError:
                                gene_name = 'undefined'
                            try:
                                transcript_id = str( sf.qualifiers['transcript_id'][0] )
                            except KeyError:
                                transcript_id = 'undefined'                            
                            
                            seqf.type = sub_feature
                            seqf.id = f'{gene_name}:{transcript_id}'
                            seqf.name = f'{parent_record.id}:{feat_id}'
                            seqf.qualifiers = exon1.qualifiers
                            seqf.description = f"{parent_record.id}:{feat_id}:gene={gene_name}:type={seqf.type}:len={seqf_length}:n_{seq_feature}={len(seqflist)}:strand={feat_strand}"
                            if srecord_dict is not None:
                                logging.debug(f'made merged, extracted SeqRecord: {seqf}')
                            else:
                                logging.debug(f'made un-merged SeqFeature: {seqf}')
                            transcript_list.append(seqf)
                        else:
                            logging.debug(f'type {subfeat.type} not in {sub_features}')
                            
                    sub_feats_handled +=1
            else:
                logging.warning(f'type {feat.type} !=  {feature}')
            feats_handled += 1
        logging.debug(f'{feats_handled} for {rec_id}')
    logging.info(f'made list of {len(transcript_list)} objs')    
    return transcript_list


def get_genes(arecord_list, srecord_dict, annot_type = 'refseq-gff'):
    '''
    annot_types = refseq-gff   refseq-gtf  ensembl-gff ensembl-gtf generic-gtf
    
    '''
    if annot_type == 'refseq-gff':
        tlist = get_genes_refseq_gff(arecord_list, srecord_dict)
    elif annot_type == 'generic-gtf':
        tlist = get_genes_generic_gtf(arecord_list, srecord_dict)
    
    return tlist
    
    
def get_genes_combo(arecord_list, 
              srecord_dict=None, 
              features=['gene','inferred_parent'],  # top level object 
              sub_features=['mRNA'],  # mid-level to handle.        
              seq_feature='exon', #  lowest-level sequences to be extracted.
              seq_flag_type = 'transcript' # custom flag indicating exon level (for special case)
              ):
    '''
    arecord_list:  list of annotation/feature object trees. 
    srecord_dict:  dictionary of seqrecords that are parent seqeunce of annotation records. 
    
    Generic GTF 3 LAYERS 
    Seqrecord(seq, id).features
        SeqFeature(type=gene | inferred_parent).sub_features
            SeqFeature(type=mRNA| type=inferred_parent).sub_features
                SeqFeature( type=exon )  | SeqFeature( type=CDS )
    
        SeqRecord -> gene -> mRNA  ->   exon|CDS            
    
    Refseq GFF   3 LAYERS
    Seqrecord(seq, id).features
        SeqFeature(type=gene | biological_region, pseudogene, enhancer |protein_binding_site ).sub_features   
            SeqFeature(type= mRNA| transcript | lnc_RNA ).sub_features
                SeqFeature( type=exon ) | SeqFeature( type=CDS )

    Refseq GTF   3 LAYERS
    Seqrecord(seq, id).features
        SeqFeature(type=gene | inferred_parent).sub_features   
                SeqFeature( type=exon ) | SeqFeature( type=CDS )
    
    
    Custom GFF/GTF (Scotomys_teguina)
    SeqRecord(seq,id).features
        SeqFeatures(type='inferred_parent).sub_features
            SeqFeature(type=transcript | exon | CDS )
     [ no sub-sub-features ]
    
        SeqRecord -> inferred_parent ->  transcript + exon(s)
                     feature               sub-feature
    '''
    transcript_list = []
    # annotation records
    arecs_handled = 0
    
    for arec in arecord_list:
        # arec is list of features on a sinlge parent SeqRecord
        rec_id = arec.id  # id of parent SeqRecord object. 
        if srecord_dict is not None:
            parent_record = srecord_dict[rec_id]
        else:
            parent_record = arec
        # SeqRecord(seq=Seq(None, length=884039), id='scf7180000321890.891056'
        logging.debug(f'handling annotation record {rec_id}')
        feats_handled = 0
        sub_feats_handled = 0
        for feat in arec.features:        
            if feat.type in features:
                logging.debug(f'handling feature {feat.id} type={feat.type}')
                feat_id = feat.id.replace('"', ':')
                feat_strand = feat.location.strand
                logging.debug(f'handling sub_features. checking for custom...')                
    
                sf_types = get_feature_type_list(feat.sub_features)
                if seq_flag_type in sf_types:
                    logging.info(f'found custom hierarchy, with {seq_flag_type}...')
                    # this is a custom, 2-level annotation set. handle special. 
                    # treat this level as a sub-sub-feature
                    seqflist = []
                    for sf in feat.sub_features:
                        if sf.type == seq_feature:
                            logging.debug(f'adding {sf.type} to sequence SeqFeature list. ')
                            seqflist.append(sf)

                    if (len(seqflist) > 0):
                        # make new SeqFeature
                        mloc = get_merged_locations(seqflist)
                        exon1 = seqflist[0]
                        seqf = SeqFeature(location=mloc)
                        seqf_length = len(seqf)
                        if srecord_dict is not None:
                            seqf = seqf.extract(parent_record)
                        try:
                            gene_name = str( exon1.qualifiers['gene_name'][0] )
                        except KeyError:
                            gene_name = 'undefined'
                        seqf.type = 'mRNA'
                        seqf.id = f'{gene_name}'
                        seqf.name = f'{parent_record.id}:{feat_id}'
                        seqf.qualifiers = exon1.qualifiers
                        seqf.description = f"{parent_record.id}:{feat_id}:gene={gene_name}:type={seqf.type}:len={seqf_length}:n_{seq_feature}={len(seqflist)}:strand={feat_strand}"
                        if srecord_dict is not None:
                            logging.debug(f'made merged, extracted SeqRecord: {seqf}')
                        else:
                            logging.debug(f'made un-merged SeqFeature: {seqf}')
                        transcript_list.append(seqf)
                    
                else:
                    # 
                    # This is a conventional, properly hierarchical feature
                    #
                    logging.info(f'found standard hierarchy, with feat_type={feat.type}')
                    sf_type = None
                    for subfeat in feat.sub_features:    
                        if subfeat.type in sub_features:
                            sf_type = subfeat.type        
                            sf_id = subfeat.id
                            ssflist = []            
                            for ssf in subfeat.sub_features:
                                if ssf.type == seq_feature:
                                    logging.debug(f'adding {ssf.type} to sub-sub-feature list. ')
                                    ssflist.append(ssf)                    
                            
                            if ( sf_type is not None ) and (len(ssflist) > 0):
                                # make new SeqFeature
                                mloc = get_merged_locations(ssflist)
                                seqf = SeqFeature(location=mloc)
                                seqf_length = len(seqf)
                                if srecord_dict is not None:
                                    seqf = seqf.extract(parent_record)
                                seqf.id = f'{feat_id}:{sf_id}'
                                seqf.type = sf_type
                                seqf.name = f'{feat_id}:{feat_strand}'
                                seqf.description = f"{parent_record.id}:{feat_id}:type={seqf.type}:len={seqf_length}:n_{seq_feature}={len(ssflist)}:strand={feat_strand}"
                                if srecord_dict is not None:
                                    logging.debug(f'made merged, extracted SeqRecord: {seqf}')
                                else:
                                    logging.debug(f'made un-merged SeqFeature: {seqf}')
                                transcript_list.append(seqf)
                        else:
                            logging.debug(f'type {subfeat.type} not in {sub_features}')
                            
                    sub_feats_handled +=1
            else:
                logging.warning(f'type {feat.type} not in {features}')
            feats_handled += 1
            
        logging.debug(f'{feats_handled} for {rec_id}')
    logging.info(f'made list of {len(transcript_list)} objs')    
    return transcript_list


def write_seqfeature_fasta(seqfeature_list, seqrecord_dict, outfile, field='description', sep=":", translate=False):
    '''
    assumes seqrecord key is first field of description
    '''
    logging.debug(f'len(seqfeature_list) = {len(seqfeature_list)}, len(seqrecord_dict) = {len(seqrecord_dict)}, outfile={outfile}, translate={translate}')
    outfile = os.path.abspath(outfile)
    seqtype = 'dna'
    if translate:
        seqtype = 'protein'
    logging.info(f'writing {len(seqfeature_list)} {seqtype} sequences to {outfile} ')
    with open(outfile, 'w') as ofh:
        for sf in seqfeature_list:
            parent_key = sf.description.split(sep)[0]
            parent_record = seqrecord_dict[parent_key]
            sr = sf.extract(parent_record)
            if translate:
                sr = sr.translate()
            sr.id = sf.id
            sr.type = sf.type
            sr.name = sf.name
            sr.description = sf.description
            #sf.description = f"{parent_record.id}:{feat_id}:type={sf.type}:n_{seq_feature}={len(ssflist)}:strand={feat_strand}"      
            SeqIO.write(sr , ofh, 'fasta')
    

def gtf_to_df(annot_file):
    
    GTF_COLS = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    df = pd.read_csv(annot_file, sep='\t', header=None, names = GTF_COLS)    
    return df

def parse_gtf_attribute(pdser, sep=';'):
    '''
    takes GTF attributes column and returns multi-column DF with key=colname, val = cellval. 
    '''
    




def make_gene_transcript_df(genome_file, annot_file, species, cache, records):
    '''
    take raw .fasta, .gtf/gff3, and produce
    dataframe. 
    sequence    gene_name    strand    
    
    '''
    #rec_list = build_genome(genome_file, annot_file, species, cache, records)
    afh = open(annot_file)
    sfh = open(genome_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(sfh, "fasta"))
    gffgen = GFF.parse(afh, base_dict=seq_dict)
    gfflist = list(gffgen)





def add_transcripts(reclist, tg_map ):
    '''
        Adds correct transcripts to .sub_features of every gene feature in SeqRecord list.
        
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


def index_by_id(record_list):
    '''
    BioPython chromosome info is in lists of SeqRecords, we often want to pull them out by id. 
    '''
    rec_idx = {}
    for sr in record_list:
        rec_idx[sr.id] = sr
    return rec_idx



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
    

def build_genome(genome_file, annot_file, species='Mus_musculus', cache=True, records=None):
    '''
      returns 
          genome dict of SeqRecords, one per chromosome/contig
          annot_dict of annotated SeqRecord objects, one per transcript
          gene_idx dict  indexed map of chr -> gene_ids -> list of transcript SeqRecord objects
    ''' 
   
    outlist = []
    
    if records is None:
        key = f'reclist.{species}.all'
        if cache_object_exists(key) and cache:
            outlist = get_cache_object(key)
        else:
            logging.info(f'loading from {genome_file} and {annot_file} for species {species}')
            reclist = load_combined(genome_file, annot_file)
            #logging.info(f'getting transcript gene mappings from ensembl ...')  
            #transmap = get_mapping_by_gene(species=species)
            #logging.info(f'adding transcripts to genes.')
            #add_transcripts(reclist, transmap)
            #logging.info(f'transcripts added.')
            #logging.debug(f'created genome with genes->transcripts: {reclist}')
            store_cache_object(reclist, key)
            for r in reclist:
                key = f'reclist.{species}.{r.id}'
                store_cache_object(r, key )
            outlist = reclist

    else:
        needed = []               
        for record in records:
            key = f'reclist.{species}.{record}'
            if cache_object_exists(key) and cache:
                obj = get_cache_object(key)
                outlist.append(obj)
            else:
                needed.append(record)
        logging.debug(f'records uncached={needed}')
        if len(needed) > 0:
            # there are some records not already cached, we might as well parse and cache 
            logging.info(f'loading from {genome_file} and {annot_file} for species {species}')
            reclist = load_combined(genome_file, annot_file, records)
            #logging.info(f'getting transcript gene mappings from ensembl ...')  
            #transmap = get_mapping_by_gene(species=species)
            #logging.info(f'adding transcripts to genes.')
            #add_transcripts(reclist, transmap)
            #logging.info(f'transcripts added.')
            #logging.debug(f'created genome with genes->transcripts: {reclist}')
            for r in reclist:
                key = f'reclist.{species}.{r.id}'
                store_cache_object(r, key )
            outlist = reclist
    return(outlist)


def save_gff3(annot_dict):
    srlist = []
    for k in annot_dict.keys():
        srlist.append(annot_dict[k])

       

def get_seqfeatures(record_list, feature_type='transcript', min_len=75, top_n = 3):
    '''
    Creates a nested dict (genes, transcripts) with lists of ONLY requested feature types with structure 
    
    Exons to be spliced are merged into lists.     
            
    [exon|intron|UTR|transcript|CDS] 
    exon: 
        get most shared (top_n) exons for each *gene*. 
    
    intron:
        get list of introns for each *gene*

    UTR:
        get all 5UTRs and 3UTRs for each *gene*
               
    transcript:
       merge all exons for each *transcript*
    
    CDS:
       merge all CDS for each *transcript*

  
    takes list of chromosome records
    finds all genes
       finds all transcripts for all genes
           finds seqfeatures of given type
    
    Returns processed dict structure of transcript SeqFeatures, indexed by chromosome id, and gene_id.
    Seqfeatures ready to be handed to concatenating routine (so preserve order). 
    record_list index necessary to extract sequence from FeatureLocations
                                                               rec_dict 
        { '1' : { 'ENSMUSG00000021743' :                       gene_dict
                    { 'ENSMUST00000224714' :                   trans_dict 
                          [ SeqFeature, SeqFeature, ... ],     individual sf_list for   
                      'ENSMUST00000022262' : 
                          [ [SeqFeature, SeqFeature], .. ],    e.g. merged exon list for transcript. 
                    },
                  '<gene2> ...
                },
                
                  'ENSMUSG00000103609': [ SeqFeature ],
                },
          '2' : ...
        }
    '''
    logging.info(f'getting all feature_type={feature_type} ')
    rec_dict = {}
    for rec in record_list:
        rec_id = rec.id
        gene_dict = {}
        for gene in rec.features:
            # only process genes at this level...
            if gene.type == 'gene':
                # will hold entries for each transcript. 
                trans_dict = {}
                logging.debug(f' handling gene = {gene.id}')
                              
                if feature_type == 'exon':
                    logging.debug(f'type is exon. finding most shared...')
                    sf_list = most_shared_exon(gene, feature_type='exon', min_len=min_len, top_n=top_n) 
                    trans_dict['exon'] = sf_list

                elif feature_type == 'intron':
                    logging.debug(f'type is intron. handling all transcripts...')
                    for trans in gene.sub_features:
                        logging.debug(f'handling transcript {trans.id}')
                        sf_list = get_introns(trans, min_len, top_n)
                        trans_dict[trans.id] = sf_list
                
                elif feature_type == 'five_prime_utr':
                    logging.debug(f'type is five_prime_utr. handling all transcripts...')
                    for trans in gene.sub_features:
                        logging.debug(f'handling transcript {trans.id}')
                        sf_list = []
                        fputr_list = []
                        for sf in trans.sub_features:
                            if sf.type == 'five_prime_utr':
                                fputr_list.append(sf)
                        sf_list.append(fputr_list)
                        trans_dict[trans.id] = sf_list

                elif feature_type == 'three_prime_utr':
                    logging.debug(f'type is three_prime_utr. handling all transcripts...')
                    for trans in gene.sub_features:
                        logging.debug(f'handling transcript {trans.id}')
                        sf_list = []
                        tputr_list = []
                        for sf in trans.sub_features:
                            if sf.type == 'three_prime_utr':
                                tputr_list.append(sf)
                        sf_list.append(tputr_list)
                        trans_dict[trans.id] = sf_list


                elif feature_type == 'CDS':
                    logging.debug(f'type is CDS. merging all CDSs...')
                    for trans in gene.sub_features:
                        logging.debug(f'handling transcript {trans.id}')
                        sf_list = []
                        cds_list = []
                        for sf in trans.sub_features:
                            if sf.type == 'CDS':
                                cds_list.append(sf)
                        sf_list.append(cds_list)
                        trans_dict[trans.id] = sf_list 

                    
                elif feature_type == 'transcript':
                    logging.debug(f'type is transcript. merging all exons...')
                    for trans in gene.sub_features:
                        logging.debug(f'handling transcript {trans.id}')
                        sf_list = []
                        ex_list = []
                        for sf in trans.sub_features:
                            if sf.type == 'exon':
                                ex_list.append(sf)
                        sf_list.append(ex_list)
                        trans_dict[trans.id] = sf_list             
                
                gene_dict[gene.id] = trans_dict
           
            else:
                logging.debug(f'gene.type not gene ({gene.type}). Ignoring.')
        
        rec_dict[rec.id] = gene_dict
    return rec_dict  

def get_merged_locations(seqfeat_list):
    '''
    
    '''
    loclist=[]
    for sf in seqfeat_list:
        loclist.append(sf.location)
    if len(loclist) > 1:
        mergedloc = CompoundLocation(loclist)
    else:
        mergedloc = loclist[0]
    return mergedloc
    


def sort_loclist(feat_location_list):
    '''
        sorts list of exact FeatureLocations by numeric position. of the start (leftmost value) 
        e.g. FeatureLocation(ExactPosition(45620850), ExactPosition(45620913), strand=-1)
        
        Uses in-place .sort() with lambda. returns none. 

        Needed to construct introns between CDS/exons. 
    '''
    feat_location_list.sort(key = lambda x: int(x.start) )
    return feat_location_list



def get_introns( transcript_sf, min_len=75, top_n=3):
    '''
        Take SeqFeature (transcript) and process sub_features.
    
        Return list of intronic SeqFeatures with valid FeatureLocations and same .qualifiers as nearby CDS. 
    
    SeqFeature(location=None, type='', location_operator='', strand=None, id='<unknown id>', 
        qualifiers=None, sub_features=None, ref=None, ref_db=None)
    
    '''
    old_cds = None
    intron_list = []    
    for sf in transcript_sf.sub_features:
        if sf.type == 'CDS':
            if old_cds is not None:
                # we have an intron
                loclist = sort_loclist( [ old_cds.location, sf.location ] )
                # get ExactPositions
                start = loclist[0].end
                end = loclist[1].start
                fl = FeatureLocation(start, end, strand=sf.strand, ref=None, ref_db=None)
                sf = SeqFeature(location=fl, type='intron', strand=sf.strand, qualifiers=old_cds.qualifiers)
                if len(sf) > min_len:
                    intron_list.append(sf)
                old_cds = sf
            else:
                old_cds = sf
    return intron_list
        

    
def most_shared_exon(gene_sf, feature_type='exon', min_len=75, top_n=3):
    '''
    Take SeqFeature (gene)
            SeqFeature.sub_features -> List of SeqFeatures  (transcripts/splice-variants)
                ->SeqFeature.sub_features -> List of SeqFeatures (exon|UTR|stop_codon|CDS|transcript)
    
    Extract exons. 
    Find exon most common to all transcripts longer than min_len             
    
    @arg SeqFeature (gene w/ transcript SeqFeature sub_features.   
    
    @return ordered list of "best" exons
         
    Preference:  most-shared feature between variants that meets minimum length. 
    Note:  Assumes all versions of a feature at a given location start are the same between transcripts.   
    
    ses_lengths >  [75, 102 , 150, 201]
        
    '''
    logging.debug(f'gene seqfeature: {gene_sf}')
    feature_dict = defaultdict(list)
    len_list = []
   
    for tsf in gene_sf.sub_features:
        if tsf.type=='inferred_parent':
            for sf in tsf.sub_features:
                if sf.type == 'exon':
                    logging.debug(f'handling exon: {sf}')
                    sf_loc = sf.location.start.position
                    sf_len = len(sf)
                    if len(sf) > int(min_len):
                        feature_dict[sf_loc].append(sf)
                        len_list.append( (sf_loc, sf_len) ) 
    
    logging.debug(f'feature_dict: {pprint.pformat(dict(feature_dict))} ' )
    
    # basically just making list, removing dupes, and ordering by size...
    logging.debug(len_list)
    len_list = list(set(len_list))
    len_list.sort(key=operator.itemgetter(1), reverse=True)
    logging.debug(len_list)                       
    len_list = len_list[:top_n]
    out_sf_list = []
    for (sf_loc, sf_len) in len_list:
        #logging.debug(f'looking up loc {sf_loc} with len {sf_len}')            
        out_sf_list.append(feature_dict[sf_loc][0])                
    logging.debug(f'found list of {len(out_sf_list)} exons > {min_len} bp. Top {top_n}. ')
    return out_sf_list  
   
MAPDICT = { 'GLT14'  :  'GALNT14',
            'RGRF2'  :  'RASGRF2',
            'GLTL6'  :  'GALNTL6',
            'CADH9'  :  'CDH9',
            'RSLAA'  :  'RASL10A',
            'SORC3'  :  'SORCS3',
            'TSH2'   :  'TSHZ2',
            'HS3S4'  :  'HS3ST4',
            'INP4B'  :  'INPP4B',
            'CAH3'   :  'CAR3',
            'LRRT4'  :  'RRTM4',
            'TEN3'   :  'TENM3',
            'IGS21'  :  'IGSF21',
            'BRNP3'  :  'BRINP3',
            'CNTP4'  :  'CNTNAP4',
            'SEM5A'  :  'SEMA5A'
            }

def apply_mapdict(row, column='r_gene', mapdict=MAPDICT):
    '''
    get alternate items based on row contents. 

    '''
    r = None
    v = row[column]
    try:
        r = mapdict[ v ]
    except KeyError:
        r = v
    return r
    


def prepare_genome_ensembl(genomefile, annotfile, outdir):
    '''
    for ensemble, no assembly report is needed. 
    make <chrnum>label.txt files for each chromosome, with <chrnum> as value inside. 
    make <chrnum>.fa and <chrnum>.fa.fai for each chromosome using index_region
        
     
    '''
    logging.debug(f'creating symlinks')
    if os.path.exists(f'{outdir}/genome.fa'):
        os.remove(f'{outdir}/genome.fa')
    os.symlink(genomefile, f'{outdir}/genome.fa')
    if os.path.exists(f'{outdir}/annotation.gtf'):
        os.remove(f'{outdir}/annotation.gtf')
    os.symlink(annotfile, f'{outdir}/annotation.gtf')

    sfh = open(genomefile)
    seq_dict = SeqIO.to_dict(SeqIO.parse(sfh, "fasta"))
    sfh.close()
    logging.debug(f'read genome file {genomefile} w/ {len(seq_dict)} records. keys={seq_dict.keys()}')
    chrlabels = list(seq_dict.keys())
    for chrlabel in chrlabels:
        fn = f'{outdir}/{chrlabel}label.txt'
        logging.debug(f'writing {fn}...')
        with open(fn, 'w') as f:
            f.write(f'{chrlabel}\n')
    
        make_chr_index( infile=genomefile,
                        genomedir = outdir, 
                        chr = chrlabel, 
                        outfile = f'{outdir}/{chrlabel}.fa')
        logging.debug(f'handled chromsome {chrlabel}')
    
    star_genome(outdir, "6" ,f'{outdir}/annotation.gtf',  f'{outdir}/genome.fa')
    samtools_dict(f'{outdir}/genome.fa',f'{outdir}/genome.dict')
    
    logging.info(f'done.')


def make_chr_index(infile, genomedir, chr, outfile):
    region = get_chr_label(genomedir, chr)
    samtools_faidx_region(infile, outfile, region)  


def get_chr_label(genomedir, chr='chrX'):
    labelfile = f"{genomedir}/{chr}label.txt"
    f = open(labelfile, 'r')
    label = f.read().strip()
    logging.debug(f"retrieved label {label} for {chr} in {labelfile}")
    return label

def make_chr_label(reportfile, outfile, chr='chrX' ):
    """
    reads NCBI assembly report and extracts selected scaffold/assembly label. 
    """
    if os.path.exists(reportfile):
        colnames = ['Sequence-Name','Sequence-Role','Assigned-Molecule',
                    'Assigned-Molecule-Location/Type','GenBank-Accn',
                    'Relationship','RefSeq-Accn','Assembly-Unit',
                    'Sequence-Length','UCSC-style-name']   
        df = pd.read_csv(reportfile, comment="#", sep='\t')
        df.columns = colnames
        if chr.startswith('chr'):
            chrnum = chr[3:]
        else:
            chrnum = chr
        tagval = chrnum
        label = df[ df['Assigned-Molecule'] == tagval]['RefSeq-Accn'].values[0]
        logging.debug(f'extracted label {label} for {tagval} in {reportfile}')
    else:
        label = chr

    f = open(outfile, 'w')
    f.write(f'{label}\n')
    f.close()

def parse_assembly_report(reportfile):
    '''
    return list of tuples  ( chrnum,   contigid )
    
    '''
    colnames = ['Sequence-Name','Sequence-Role','Assigned-Molecule',
                'Assigned-Molecule-Location/Type','GenBank-Accn',
                'Relationship','RefSeq-Accn','Assembly-Unit',
                'Sequence-Length','UCSC-style-name']   
    df = pd.read_csv(reportfile, comment="#", sep='\t', header=None)
    df.columns = colnames
    amdf = df[df['Sequence-Role'] == 'assembled-molecule']
    # idiosyncratic between some refseq genomes and others...
    ndf = amdf[amdf['Assembly-Unit'] == 'Primary Assembly' ]
    if len(ndf) > 2:
        amdf = ndf  
    mdf = amdf[['Assigned-Molecule', 'RefSeq-Accn']]
    tlist = list(mdf.itertuples(index=False, name=None))
    logging.debug(f'extracted {len(tlist)} maps in {reportfile}')
    return tlist

def prepare_genome_refseq(genomefile, annotfile, reportfile, outdir):
    '''
    sets up genome for per-chromosome work, and STAR alignments. 
        make <chrnum>label.txt files for each chromosome, with <chrnum> as value inside. 
    make <chrnum>.fa and <chrnum>.fa.fai for each chromosome using index_region
    
    '''
    logging.debug(f'creating symlinks')
    if os.path.exists(f'{outdir}/genome.fa'):
        os.remove(f'{outdir}/genome.fa')
    os.symlink(genomefile, f'{outdir}/genome.fa')
    
    if os.path.exists(f'{outdir}/annotation.gtf'):
        os.remove(f'{outdir}/annotation.gtf')
    os.symlink(annotfile, f'{outdir}/annotation.gtf')
    
    tlist = parse_assembly_report(reportfile)
    logging.debug(tlist)
    for (chrlabel, contig) in tlist:
        fn = f'{outdir}/{chrlabel}label.txt'
        logging.debug(f'writing {fn}...')
        with open(fn, 'w') as f:
            f.write(f'{contig}\n')

        make_chr_index( infile=genomefile,
                        genomedir = outdir, 
                        chr = chrlabel, 
                        outfile = f'{outdir}/{chrlabel}.fa')
        logging.debug(f'handled chromosome {chrlabel}')
        
    star_genome(outdir,"6" ,f'{outdir}/annotation.gtf',  f'{outdir}/genome.fa')
    samtools_dict(f'{outdir}/genome.fa',f'{outdir}/genome.dict')
    samtools_faidx(f'{outdir}/genome.fa',f'{outdir}/genome.fa.fai')
    logging.info(f'done.')

def prepare_genome_genbank():
    pass


def samtools_faidx_region(infile, outfile, region):
    cmd = ['samtools',
           'faidx',
           '-o', outfile,
           infile,
           region 
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise   



 