#!/usr/bin/env python
#
#
#  https://www.biostars.org/p/224572/
#  preferred breakdown:    animals,   fungi,  plants,   protists,   bacteria,  archaea
#  ncbi taxon category codes:  chordata     7711   vertebrata   7742  Insecta  50557 
#                                # in DB
#   1 => 'Eukaryota',
#        2 => 'Animals',         1653
#          Vertebrata            877
#            7 => 'Mammals', 
#            3 => 'Birds',
#            9 => 'Reptiles',
#            2 => 'Amphibians', 
#            4 => 'Fishes', 
#          Invertebrates  
#            5 => 'Flatworms',    38
#            10 => 'Roundworms'   119
#            6 => 'Insects',      412
#            8 => 'Other Animals', },
#        3 => 'Fungi', 
#        4 => 'Other', 
#        5 => 'Plants',         459
#        6 => 'Protists'},DB
#   2 => 'Bacteria' -> all, 
#   3 => 'Archaea',
#   4 => 'Viroids', 
#   5 => 'Viruses');
#   
#   areas?:   eukaryote->animals
#                          amph, birds, fishes, mammals, reptiles   |  flatworms , roundworms 
#            eukaryote->plants               
#            eukaryotes->fungi
#            eukaryotes->protists
#   
#  Convert newick tree to distance matrix in R:  https://www.biostars.org/p/312148/
#    treeText <- readLines(tree.phy)
#        treeText <- paste0(treeText, collapse="")
#        library(treeio)
#        tree <- read.tree(text = treeText) ## load tree 
#    distMat <- cophenetic(tree) ## generate dist matrix
#
#
#
