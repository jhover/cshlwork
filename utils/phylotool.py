#!/usr/bin/env python
#  
#  phylogenetic processing. 
#  ncbi_get_genome_tree  
#                                # in DB
#  1 => 'Eukaryota',
#        2 => 'Animals',
#            7 => 'Mammals', 
#            3 => 'Birds',
#            9 => 'Reptiles',
#            2 => 'Amphibians', 
#            4 => 'Fishes', 
#            5 => 'Flatworms',    38
#            10 => 'Roundworms'   119
#            6 => 'Insects', 
#            8 => 'Other Animals', },
#        3 => 'Fungi', 
#        4 => 'Other', 
#        5 => 'Plants', 
#        6 => 'Protists'},
#         
#   2 => 'Bacteria' -> all, 
#   3 => 'Archaea',
#   4 => 'Viroids', 
#   5 => 'Viruses');
#   
#   areas:   eukaryote->animals
#                          amph, birds, fishes, mammals, reptiles   |  flatworms , roundworms 
#            eukaryote->plants
#                        
#            eukaryotes->fungi
#            eukaryotes->protists
#
#
#   codes:  chordata     7711   vertebrata   7742
#   
#
