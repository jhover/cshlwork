.libPaths("/usr/lib64/R/library")
DATADIR = "/home/ftp/data/HiC"
# t install new R packages - use
# sudo su - -c "R -e \"install.packages('shiny', repos= 'http://cran.rstudio.com/' )\""
# sudo su - -c "R -e \"BiocManager::install('rhdf5')\""
# sudo su - -c "R -e \"devtools::install_github('gillislab/MetaMarkers')\""


# don't load too many packages in global.R - slows initial load time for the app
# Instead load them in server.R
library(shiny,quietly=T)
library(shinyWidgets,quietly=T)
library(shinydashboard,quietly=T)
library(shinydashboardPlus,quietly=T)
library(shinyBS,quietly=T)  
library(shinyjs,quietly=T)

SPECIESLIST = c("human", "mouse", "drosophila")
METRICS <- c("pearson", "coexpression", "shared_contacts", "contacts")
names(METRICS) <- c("Pearson Correlation", "Co-Expression",
                    "Shared Contacts", "Contact Frequency")
RESORDER <- 1:6
names(RESORDER) <- c('1kbp', '25kbp', '40kbp', '100kbp', '250kbp', '500kbp')
