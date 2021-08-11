#!/bin/bash
#
# Force downloads all .dat files from uniprot. 
#
#
#



CURRENT="https://ftp.uniprot.org/pub/databases/uniprot/current_release"
CDATS="uniprot_sprot uniprot_trembl"
TDATS="archaea bacteria fungi human invertebrates mammals plants rodents vertebrates viruses" 

#echo "Handling complete databases..."
#for dat in $CDATS; do 
#    echo time wget -O ${dat}.dat.gz  ${CURRENT}/knowledgebase/complete/${dat}.dat.gz
#    time wget -O ${dat}.dat.gz  ${CURRENT}/knowledgebase/complete/${dat}.dat.gz
#    echo gunzip -f ${dat}.dat.gz 
#    gunzip -f ${dat}.dat.gz 
#done

#echo "Merging sprot/trembl dbs..."
#echo "cat uniprot_trembl.dat uniprot_sprot.dat > uniprot_all.dat" 
#cat uniprot_trembl.dat uniprot_sprot.dat > uniprot_all.dat 


echo "Handling taxonomic division databases..."
for tdat in $TDATS; do
    echo "Handling ${tdat} taxa/taxon..."
    for cat in $CDATS; do
    	echo time wget -O ${cat}_${tdat}.dat.gz  ${CURRENT}/knowledgebase/taxonomic_divisions/${cat}_${tdat}.dat.gz
    	time wget -O ${cat}_${tdat}.dat.gz  ${CURRENT}/knowledgebase/taxonomic_divisions/${cat}_${tdat}.dat.gz
    	echo gunzip -f ${cat}_${tdat}.dat.gz 
    	gunzip -f ${cat}_${tdat}.dat.gz 
    done
done

echo "Merging sprot/trembl dbs..."
for tdat in $TDATS; do
    echo "cat uniprot_trembl_${tdat}.dat uniprot_sprot_${tdat}.dat > uniprot_all_${tdat}.dat" 
    cat uniprot_trembl_${tdat}.dat uniprot_sprot_${tdat}.dat > uniprot_all_${tdat}.dat 
done


