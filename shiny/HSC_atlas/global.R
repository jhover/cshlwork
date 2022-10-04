# no need to change directory. Just copy global.R, ui.R and server.R to 'demo-app' directry.

.libPaths("/usr/lib64/R/library")
# .libPaths("/home/johlee/R/x86_64-redhat-linux-gnu-library/4.0")


DATADIR = "/var/ftp/data/HSC_atlas/"
GO_LIST = rhdf5::h5ls(paste0(DATADIR , 'go2gene.hdf5') )$name


TAB1_LOOM_FILE_PATH = paste0(DATADIR,"erythroid_and_monocyte_lineage_adata_no_gaps.loom")
TAB1_H5AD_FILE_PATH = paste0(DATADIR,"erythroid_and_monocyte_lineage_adata_no_gaps.h5ad")

TAB1_MARKER_PATH = paste0(DATADIR, "droplet_scNym_metamarkers.csv.gz")
TAB1_GO_ENRICH = paste0(DATADIR,"cluster_GO_enrichment.csv")

GENELIST = rhdf5::h5read(TAB1_LOOM_FILE_PATH, 'row_attrs/var_names')


TAB2_LOOM_FILE_PATH = paste0(DATADIR,"mouse_hsc_labeled.loom")
TAB2_H5AD_FILE_PATH = paste0(DATADIR,"mouse_hsc_labeled.h5ad")


# TAB2_GENELIST = rhdf5::h5read(TAB1_LOOM_FILE_PATH, 'row_attrs/var_names')
TAB2_MARKER_PATH = paste0(DATADIR, "droplet_FACS_labels_metamarkers.csv.gz")
TAB2_GO_ENRICH = paste0(DATADIR,"cell_state_GO_enrichment.csv")



# Load in Rdata  
# load(paste0(dir,'www/attr.Rdata'))          # known gene names and various IDs for each species
# load(paste0(dir,'www/gotree.Rdata'))        # Gene ontology tree (first few layers)
# load("/ftp/data/gene2go/voc.Rdata")           # descrptions of GO terms
# load(paste0(dir,'www/goa.Rdata'))          # gene x  function binaries for each species
# load(paste0(direc,'www/EGADlite.Rdata'))
# load('www/GOlist.Rdata')        # cc / bp / mf trees


# enableBookmarking(store = "server")

##### config file #######
# /etc/shiny-server/shiny-server.conf

##### To install new packages ...
#  sudo su - -c "R -e \"install.packages('shiny', repos= 'http://cran.rstudio.com/' )\""
# repos='http://cran.rstudio.com/'
#####  
#  (GRCh38) , (GRCm38) (TAIR10) (UMD3.1) (WBcel235) (GRCz11) (BDGP6) (Gallus_gallus-5.0) (Glycine_max_v2.0) (Rnor_6.0) (R64-1-1) (Sscrofa11.1) (rIRGSP-1.0) (maiAGPv4)




# speciesList[3:length(speciesList)] =  speciesList[ind+2]

# library(SingleCellExperiment)
# library(tidyverse)
# library(scattermore)
# library(visNetwork,quietly = T)
# library(Matrix,quietly = T)
# library(tools,quietly=T)
# library(ggplot2,quietly=T)
# library(taRifx,quietly=T)
# library(dplyr,quietly=T)
### library(plotly,quietly=T) #* 
# library(GO.db,quietly=T)
# library(data.table,quietly=T)
library(shiny,quietly=T)
library(shinyWidgets,quietly=T)
library(shinydashboard,quietly=T)
# library(shinyTree,quietly=T)
library(shinydashboardPlus,quietly=T)
library(shinyBS,quietly=T)  # bootstrap
# library(shinyalert,quietly=T)
# library(rhdf5,quietly=T)
# library(matrixStats,quietly=T)
library(RColorBrewer)
# library(DT,quietly=T)
# library(rARPACK) 
# library(sna)
# library(shinythemes,quietly=T)
# library(EGAD)
# library(shinyEventLogger)
# library(igraph,quietly=T)
# library(Matrix,quietly = T)
# library(tools,quietly=T)
library(ggplot2,quietly=T)
# library(SummarizedExperiment)
library(tidyverse)
