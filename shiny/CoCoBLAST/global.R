# no need to change directory. Just copy global.R, ui.R and server.R to 'demo-app' directry.
direc = "/srv/shiny-server/sample-apps/CoCoBLAST/" 
dataDirec= "/var/ftp/data/supplementData/"
consDirec= "/var/ftp/data/CoCoBLAST/coexp_cons/"

.libPaths("/usr/lib64/R/library")

# Load in Rdata  
# load(paste0(dir,'www/attr.Rdata'))          # known gene names and various IDs for each species
# load(paste0(dir,'www/gotree.Rdata'))        # Gene ontology tree (first few layers)
# load("/ftp/data/gene2go/voc.Rdata")           # descrptions of GO terms
# load(paste0(dir,'www/goa.Rdata'))          # gene x function binaries for each species
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
speciesList = list( 
    "Homo sapiens - Human"                  ='human',
    'Mus musculus - Mouse'                  ='mouse',
    "Arabidopsis thaliana - Thale cress"    ="arabidopsis",                   
    "Bos taurus - Cow"                      ="cow", 		   
    "Caenorhabditis elegans - Roundworm"    ="roundworm",
    "Danio rerio - Zebrafish"               ="zebrafish",
    "Drosophila melanogaster - Fruitfly"    ="fruitfly", 
    "Gallus gallus - Chicken"               ="chicken",
    "Glycine max - Soybean"                 ="soybean",
    "Rattus norvegicus - Rat"               ="rat",
    "Saccharomyces cerevisiae - Budding yeast"  ="yeast",
    "Sus scrofa - Boar"                     ="boar",
    "Oryza sativa - Rice"                   ="rice",
    "Zea mays - Maize"                      ="maize",
    "Salmo salar - Salmon"                  ='atlanticsalmon',
    'Apis mellifera - Honey Bee'            ='bee',
    'Bombyx mori - Silkworm'                ='bombyx',
    'Brachypodium distachyon - Brome'       ='brome',
    'Pan troglodytes - Chimpanzee'          ='chimp',
    'Chlamydomonas reinhardtii - Chlamydomonas'  ='chlam',
    'Macaca fascicularis - Crab-eating macaque' ='crabm',
    'Canis lupus familiaris - Dog'          ='dog',
    # 'Schizosaccharomyces pombe - Fission yeast' ='fissionyeast',
    'Capra hircus - Goat'                   ='goat',
    'Vitis vinifera - Grape'                ='grape',
    'Equus caballus - Horse'                ='horse',
    'Medicago truncatula - Medicago'        ='medicago',
    'Aedes aegypti - Mosquito'              ='mosquito',
    'Solanum tuberosum - Potato'            ='potato',
    'Oryctolagus cuniculus - Rabbit'        ='rabbit',
    'Macaca mulatta - Rhesus macaque'       ='rhesusm',
    'Ovis aries - Sheep'                    ='sheep',
    'Sorghum bicolor - Sorghum'             ='sorghum',
    'Solanum lycopersicum - Tomato'         ='tomato',
    # 'Xenopus tropicalis - Western clawed frog'  ='wfrog',
    # 'Xenopus tropicalis - African clawed frog'  = 'afrog',
    # 'Neurospora crassa - Neurospora'        = 'neurospora',
    'Brassica rapa - Mustard'               = 'mustard',
    'Malus domestica - Apple'               = 'apple',
    'Nicotiana tabacum - Tobacco'           = 'tobacco',
    'Oncorhynchus mykiss - Rainbow Trout'   = 'rainbowtrout'
)

# should sort by phylogeny - alphabetical for now
speciesList=speciesList[c(1:2, 2+order(names(speciesList[3:length(speciesList)])))]


# speciesList[3:length(speciesList)] =  speciesList[ind+2]

# library(SingleCellExperiment)
# library(tidyverse)
# library(scattermore)
library(visNetwork,quietly = T)
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
library(shinyalert,quietly=T)
# library(rhdf5,quietly=T)
library(matrixStats,quietly=T)
library(RColorBrewer)
# library(DT,quietly=T)
# library(rARPACK) 
# library(sna)
# library(shinythemes,quietly=T)
# library(EGAD)
# library(shinyEventLogger)
library(igraph,quietly=T)
library(Matrix,quietly = T)
library(tools,quietly=T)
library(ggplot2,quietly=T)

# load('www/attr/gene_annotations_mouse.Rdata')
# assign('mouse_attr',attr)
# load('www/attr/gene_annotations_human.Rdata')
# assign('human_attr',attr)
# load('www/attr/gene_annotations_arabidopsis.Rdata')
# assign('arabidopsis_attr',attr)
# load('www/attr/gene_annotations_rat.Rdata')
# assign('rat_attr',attr)
# load('www/attr/gene_annotations_maize.Rdata')
# assign('maize_attr',attr)
# load('www/attr/gene_annotations_chicken.Rdata')
# assign('chicken_attr',attr)
# load('www/attr/gene_annotations_fruitfly.Rdata')
# assign('fruitfly_attr',attr)
# load('www/attr/gene_annotations_soybean.Rdata')
# assign('soybean_attr',attr)
# load('www/attr/gene_annotations_cow.Rdata')
# assign('cow_attr',attr)
# load('www/attr/gene_annotations_yeast.Rdata')
# assign('yeast_attr',attr)
# load('www/attr/gene_annotations_boar.Rdata')
# assign('boar_attr',attr)
# load('www/attr/gene_annotations_roundworm.Rdata')
# assign('roundworm_attr',attr)
# load('www/attr/gene_annotations_rice.Rdata')
# assign('rice_attr',attr)
# load('www/attr/gene_annotations_zebrafish.Rdata')
# assign('zebrafish_attr',attr)
# rm(attr)


# Location of network data
# loc = 'www/rankNormVals/' 






########################## Graph, description box footer and static slider ##############################
# tabPanel("Result3",icon = icon("bezier-curve"),
#          boxPlus(title = "boxPlus with footer and static slider",width = '',height='470px',status = "success",
#                  closable=F,collapsible = T,footer_padding = F,
#                  fluidRow(column(width=9,plotOutput('tab3plot')),
#                           column(width=3,
#                                  checkboxInput("somevalue2", "Apply some condition", FALSE),
#                                  verbatimTextOutput("value2"),
#                                  sliderInput("tab3slider","Number of observations:",min = 0,max = 1000,value = 200),
#                                  actionButton('tab3SliderButton','Update changes',icon=icon('magic'),width='auto')
#                           )),
#                  footer = fluidRow(
#                      column(width=3,
#                             descriptionBlock(
#                                 number = "17%",
#                                 number_color = "green",
#                                 number_icon = "fa fa-caret-up",
#                                 header = "header1",
#                                 text = "text1")),
#                      column(width=3,
#                             descriptionBlock(
#                                 number = "18%",
#                                 number_color = "yellow",
#                                 number_icon = "fa fa-caret-down",
#                                 header = "header4",
#                                 text = "text4",
#                                 right_border = F))
#                  )
#          )
# )


########################## Graph, description box and dynamic slider ##############################
# tabPanel("Result1",icon = icon("chart-bar"),
#          boxPlus(title = "boxPlus with dynamic sidebar",width='',height='520px',status = "warning",
#                  collapsible = T,closable = T,enable_sidebar = TRUE,sidebar_width = 25,
#                  sidebar_start_open = T,footer_padding = F,sidebar_background = "#001F3F",
#                  sidebar_content = tagList(
#                      checkboxInput("somevalue", "Apply some condition", FALSE),
#                      verbatimTextOutput("value"),
#                      sliderInput("slider_boxsidebar","Number of observations:",
#                                  min = 0,max = 1000,value = 700),
#                      actionButton('tab1SliderButton','Update changes',icon=icon('magic'),width='auto')
#                  ),
#                  plotOutput("boxSidebarPlot",width='auto',height='480px'))),



####################################### Gradient Box ####################################################
# boxPlus(title = "boxPlus with gradient box and boxPad",width = '',height='520px',status = "danger",
#         closable=F,collapsible = T,footer_padding = T,
#         fluidRow(column(width=9,
#                         gradientBox(plotOutput('tab2plot',height='300px'),
#                                     title='Gradient box with boxPad',icon=icon('heart'),
#                                     gradientColor='red',boxToolSize='xs',closable=T,
#                                     width=12,height='450px',collapsible =T,footer_padding=F,
#                                     footer=sliderInput("slider_gradient","Number of observations:",
#                                                        min = 0,max = 1000,value = 700))),
#                  column(width=3,height='480px',
#                         boxPad(color = "green", height='480px',
#                                descriptionBlock(header = "8390",text = "VISITS",right_border = F),
#                                descriptionBlock(header = "30%",text = "REFERRALS",right_border = F),
#                                descriptionBlock(header = "70%",text = "ORGANIC",right_border = F)
#                         )
#                  ))
# )



#######################################Different Tab views#######################################
# #tabItem(tabName = "ortho",
# navbarPage(title = "1-1 Ortholog Finder",
#            tabPanel('Input',
#                     verticalLayout(
#                         selectInput(inputId="species",label=h3("Select species:"),choices=geneList,selected=NULL),
#                         selectInput(inputId="userGenes",label=h3("Select genes:"),multiple = T,choices=speciesList,selected=NULL),
#                         textInput("geneList", h3("Enter gene"),value = "",placeholder ="Multiple genes should be separated by ','"),
#                         fileInput("geneFile", h3("Upload gene list"), buttonLabel = 'Upload'),
#                         actionButton("action", "Generate Results",icon('cogs'))
#                     )
#            ),
#            tabPanel("Result1",value='tab1',icon=icon("chart-bar"),"Result 1 displayed here..."),
#            tabPanel("Result2",value='tab2',icon=icon("hubspot"),"Result 2 displayed here..."),
#            tabPanel("Result3",value='tab3',icon=icon("bezier-curve"),"Result 3 displayed here...")
# )
# )
#
# # net3
# tabItem(tabName = "net3",
#         fluidPage(
#             titlePanel("Co-expression network generator"),
#             fillRow(flex=2,
#                     verticalLayout(
#                         selectInput(inputId="species",label=h3("Select species:"),choices=speciesList,selected=NULL),
#                         selectInput(inputId="userGenes",label=h3("Select genes:"),multiple = T,choices=geneList,selected=NULL),
#                         textInput("geneList", h3("Enter gene"),
#                                   value = "",
#                                   placeholder ="Multiple genes should be separated by ','"),
#                         fileInput("geneFile", h3("Upload gene list"), buttonLabel = 'Upload'),
#                         actionButton("action", "Generate Results",icon('cogs'))),
#                     tabsetPanel(
#                         #side = "left", height = "450px",id = 'netTabResult',
#                         #selected = "Result1", width = '',
#                         tabPanel("Result1",icon = icon("chart-bar"),"Result 1 displayed here..."),
#                         tabPanel("Result2",icon = icon("hubspot"),"Result 2 displayed here..."),
#                         tabPanel("Result3",icon = icon("bezier-curve"),"Result 3 displayed here...")
#                         )
#                     )
#             )
#         ),

################################----individual histogram functions------################################################
################################----GO group nestedlist generator----#######################################
# treeGenerator <- function(children,ori){
#     for(i in seq(length(children))){
#         child = children[[i]]
#         if(!is.null(ori[[child]])){
#             children[[child]] = treeGenerator(ori[[child]],ori)}}
#     return(children)
# }
#
# lister <- function(vec){
#     for(i in seq(length(vec))){
#         vec[[i]]=sapply(as.vector(vec[[i]]),function(x) list(x))}
#     return(vec)
# }
#
# cc <- as.list(GOCCCHILDREN)
# cc <- cc[!is.na(cc)]
# cc <- lister(cc)
# cc_list = treeGenerator(cc[['GO:0005575']],cc)
# bp <- as.list(GOBPCHILDREN)
# bp <- bp[!is.na(bp)]
# bp <- lister(bp)
# bp_list = treeGenerator(bp[['GO:0008150']],bp)
# mf <- as.list(GOMFCHILDREN)
# mf <- mf[!is.na(mf)]
# mf <- lister(mf)
# mf_list = treeGenerator(mf[['GO:0003674']],mf)

# b1 = list()
# for(i in seq(nrow(res5))[-c(144,153:154)]){
#     #a = res5[i,][!is.na(res5[i,])]
#     a = z[i,][!is.na(z[i,])]
#     print(i)
#     switch (length(a),
#         '1' = {b1[[a[1]]]=list()},
#         '2' = {
#             if(is.null(b1[[a[1]]])){
#                 b1[[a[1]]]=list()
#                 b1[[a[1]]][[a[2]]] = a[2]
#             } else {
#                 b1[[a[1]]][[a[2]]] = a[2]
#             }
#         },
#         '3' = {
#             if(is.null(b1[[a[1]]])){
#                 b1[[a[1]]]=list()
#                 b1[[a[1]]][[a[2]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]] = a[3]
#             } else if(is.null(b1[[a[1]]][[a[2]]])) {
#                 b1[[a[1]]][[a[2]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]] = a[3]
#             } else {
#                 b1[[a[1]]][[a[2]]][[a[3]]] = a[3]
#             }
#         },
#         '4' = {
#             if(is.null(b1[[a[1]]])){
#                 b1[[a[1]]]=list()
#                 b1[[a[1]]][[a[2]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]][[a[4]]] = a[4]
#             } else if(is.null(b1[[a[1]]][[a[2]]])) {
#                 b1[[a[1]]][[a[2]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]][[a[4]]] = a[4]
#             } else if(is.null(b1[[a[1]]][[a[2]]][[a[3]]])) {
#                 b1[[a[1]]][[a[2]]][[a[3]]] = list()
#                 b1[[a[1]]][[a[2]]][[a[3]]][[a[4]]] = a[4]
#             } else {
#                 b1[[a[1]]][[a[2]]][[a[3]]][[a[4]]] = a[4]
#             }
#         })
# }

