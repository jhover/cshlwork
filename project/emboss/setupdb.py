#!/usr/bin/env python
#
# Downloads .dat.gz files for all taxa to subdir.
# Extracts .dat files. 
# Combines sprot and trembl into all
# Splits .dat files too big for indexing. 
# Indexes DB
#
#

# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
# uniprot_sprot_mammals.dat.gz
# uniprot_trembl_mammals.dat.gz 

TAXA_URL='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions'
TAXA = ['archaea',
        #'bacteria',  # way too big. 
        'fungi',
        'human',
        'invertebrates',
        'mammals',
        'plants',
        'rodents',
        'vertebrates',
        'viruses'        
        ]
REPOS=['sprot', 'trembl']
DBROOT=os.path.expanduser('~/data/embossdb')


for t in TAXA:
    for r in REPOS:
        DLURL=f'{TAXA_URL}/uniprot_{r}_{t}.dat.gz'
        
    
