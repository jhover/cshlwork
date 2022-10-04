library(shiny, quietly=T)
library(shinyWidgets, quietly=T)
library(shinydashboard, quietly=T)
library(shinydashboardPlus, quietly=T)
library(shinyBS, quietly=T)
library(shinycssloaders, quietly=T)
# require(GenomicRanges, quietly=T)
# require(GenomicFeatures, quietly=T)
# require(AnnotationDbi, quietly=T)
# require(Gviz, quietly=T)
# require(stringr, quietly=T)
# require(ComplexHeatmap, quietly=T)
require(viridis, quietly=T)
require(DT, quietly=T)
options(ucscChromosomeNames=FALSE)
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)


# dataDirec= "/var/ftp/data/Meta_scATAC/"
.libPaths("/usr/lib64/R/library")

dataset.list = c('BICCN', 'Cusanovich 2018', 'Lareau 2019', 'Chen 2019', 'Spektor 2019', 'Zhu 2019', 'Preissl 2018', 'Predicted by DL (Kawaguchi RK, et al.)')
dataid.list = c('BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990', 'GSE100033', 'DeepLearning')
names(dataid.list) = dataset.list
celltype.list = c('EX_L2.3.IT', 'EX_L5.6.NP', 'EX_L5.ET', 'EX_L5.IT', 'EX_L6b', 'EX_L6.CT', 'EX_L6.IT.Car3', 'IN_Lamp5', 'IN_Pvalb', 'IN_Sncg', 'IN_Sst', 'IN_Vip')

tss <<- NULL
txdb <<- NULL
gtxdb <<- NULL
gdf <<- NULL
gene.dict <<- list()
gtrack <<- NULL
grtrack <<- NULL
gr <<- list()
initialized <<- FALSE
data.label <<- list()

cgr <<- list()
cdata.label <<- list()
cluster.label.list <<- list()
fname = "./www/cluster_names.rds"
if (file.exists(fname)) {
    cluster.label.list <<- readRDS(fname)
}
