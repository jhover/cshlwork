import copy
import itertools
import logging
import os 
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def srlist_print(srlist):
    s = '['
    for sr in srlist:
        s += f"{sr.id} "
    s += ']'
    return s


def process_cross_files(filelist, outdir):
    '''
       
    '''
    filesets = []
    for filepath in filelist:
        filepath = os.path.abspath(filepath)    
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)
        
        srlist = []
        for seq_record in SeqIO.parse(filepath, 'fasta'):
            srlist.append(seq_record)
            
            #seqids.append(seq_record.id)
            #seqdescs.append( seq_record.description)
            #seqs.append(str(seq_record.seq))
        filesets.append(srlist)
        logging.debug(f'file {filename} had {len(srlist)} sequences')
    
    logging.debug(f'created {len(filesets)} sets of records. ')
    for idx in range(0,len(filesets)):
        logging.debug(f'len fileset {idx} = {len(filesets[idx])} first entry: {filesets[idx][0]}')
            
    if outdir is None:
        outdir = os.getcwd()
    else:
        outdir = os.path.abspath(outdir)    
    
    combos = []
    sourceset = filesets[0]
    for idx in range(0, len(filesets) - 1):
        logging.debug(f'loop idx={idx} len(sourceset) = {len(sourceset)}')
        targetset = filesets[idx + 1]
        logging.debug(f'handling targetset {idx + 1} len={len(targetset)}')
        combos = []
        for sr in sourceset:
            logging.debug(f'type source sr= {type(sr)}')
            for t in targetset:
                ssr = copy.copy(sr)
                logging.debug(f'target t = {t.id}')        
                #combo = list(itertools.chain.from_iterable([s,t]))
                if isinstance(ssr, SeqRecord ):
                    #logging.debug(f'type of sr is {type(sr)}')
                    logging.debug(f'combo is SeqRecord ssr={ssr.id} t={t.id}')
                    combos.append([ssr,t])
                elif isinstance(sr, list):
                    logging.debug(f'ssr is list: {srlist_print(ssr)}, appending...')
                    ssr.append(t)
                    logging.debug(f'len(ssr) after append: {len(ssr)}')
                    combos.append(ssr)
                else:
                    logging.warning(f'ssr is {ssr}')
        
        sourceset = combos    
        logging.debug(f'len(sourceset) = {len(sourceset)}')
    allcombos = sourceset                
    logging.debug(f'created set of {len(allcombos)} outsets.. Example {srlist_print(allcombos[0])} ')            
        
    outfile = f'{outdir}/{base}.fa'
    logging.debug(f'filepath={filepath} filename={filename} outdir={outdir} base={base} ext={ext} outfile={outfile}')        
    

    for combo in allcombos:
        logging.debug(srlist_print(combo))
    
        seqids = []
        seqdescs = []
        seqs = []
        for seq_record in combo:
            seqids.append(seq_record.id)
            #seqdescs.append( seq_record.description)
            #seqs.append(str(seq_record.seq))
        
        newid = 'x'.join(seqids)
        # avoid backslash in filename confuses shell...
        newid = newid.replace('/','--')
        #newdesc = ':'.join(seqdescs)
        #ss = ':'.join(seqs)
        outfile = f'{outdir}/{newid}.fa'
      
        with open(outfile, 'w' ) as fh:
            logging.debug(f'writing to {outfile}')
            SeqIO.write(combo,fh, 'fasta' )
        logging.info(f'wrote {len(combo)} sequences to {outfile}, id={newid}') 