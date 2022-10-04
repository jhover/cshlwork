tweaks <- 
  list(tags$head(tags$style(HTML("
                                 .multicol { 
                                   height: 200px;
                                   font: Arial;
                                   -webkit-column-count: 5; /* Chrome, Safari, Opera */ 
                                   -moz-column-count: 5;    /* Firefox */ 
                                   column-count: 5; 
                                   -moz-column-fill: auto;
                                   -column-fill: auto;
                                   overflow-x: scroll;
                                 } 
                                 ")) 
                                 ))



body <- dashboardBody(
    tabItems(
        tabItem(
            tabName = 'main',
            fluidPage(
                titlePanel(
                    h3(div('Cell-type specific scATAC-seq profiles at a glance'), align = "center")
                ),
                fluidRow(
                    box(
                        status='primary',width = NULL,
                        solidHeader=F,
                        column(width=8,
                            p(
                            h4(strong("Background: Heterogeneity and sparsity of scATAC-seq datasets")),
                            HTML(
                            "Single-cell Assay for Transposase Accessible Chromatin using sequencing (scATAC-seq) is a valuable resource to learn cis-regulatory elements such as cell-type specific enhancers and transcription factor binding sites. However, cell-type identification of scATAC-seq data is known to be challenging due to the heterogeneity derived from different protocols and the high dropout rate."),
                            h4(strong("Key idea: Meta-analytic marker genes for brain cell types")),
                            HTML(
                            "A key idea of our study is using meta-analytic marker genes, called SF marker set, derived from multiple scRNA-seq datasets from Brain Initiative Cell Census Network (BICCN). We developed a method, MetaMarkers, and used it to define a new meta-analytic marker gene set, which is the expanded marker set that includes co-expressed genes of marker genes. The gene set is made by selecting the top 100 genes to predict 85 cortical cell types accurately by their aggregated signals."),
                            h4(strong("Integrative server: Cell-type specific scATAC-seq profiles from six studies")),
                            HTML(
                            "The heterogeneity of scATAC-seq datasets causes the disparity in terms of the number of cells and clusters, suggesting cluster-level annotation is inadequate and not comparable. Instead, we predicted the cells whose aggregated module activity is ranked within the top 500 using SF marker gene sets. Then we produced the meta scATAC-seq profiles based on those cell-type annotation.
                            In this server, the average read count of each genomic location can be visualized in a genome browser for certain cell types or datasets."),
                            h4(strong("Deep learning prediction: Sequence-dependent and -independent chromatin accessibility regulation")),
                            HTML("Furthermore, to assess the potential of genomic sequences to regulate the cell-type specific cis-regulatory programs, we integrated predicted chromatin accessibility data using a sequence-based deep CNN trained on the BICCN scATAC-seq data.")))
                    ),
                    box(
                        status='primary',width = NULL,
                        solidHeader=F,

                        column( width=8,
                            h4(strong("Example")),
                            HTML("1. Cross dataset analysis"),
                            br(),
                            div(img(src="Server_Fig3.png", height='100%',width='100%'), align='center'),
                            br(),
                            HTML("2. Cross cell-type analysis"),
                            br(),
                            div(img(src="Server_Fig4.png", height='450',width='100%'),
                                align='center')
                        )
                    ),                      
                    box(
                        status='primary',width = NULL,
                        solidHeader=F,
                        column(width=8,
                            p(h4(strong("Reference")),
                                HTML('<ul><li>Kawaguchi RK, et al. Exploiting marker genes for robust classification and characterization of single-cell chromatin accessibility. bioRxiv (2021)</li>
                                      <li>BRAIN Initiative Cell Census Network (BICCN): Yao, Z., et al. A transcriptomic and epigenomic cell atlas of the mouse primary motor cortex. Nature 598, 103–110 (2021)</li>
                                      <li>Preissl, S., et al. Single-nucleus analysis of accessible chromatin in developing mouse forebrain reveals cell-type-specific transcriptional regulation. Nature neuroscience, 21(3):432-439 2018.</li>
                                      <li>Cusanovich DA., et al. A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. Cell, 23;174(5):1309-1324.e18 2018.</li>
                                      <li>Lareau, CA., et al. Droplet-based combinatorial indexing for massive-scale single-cell chromatin accessibility. Nature Biotechnology 37(8):916-924 2019.</li>
                                      <li>Chen, S., et al. High-throughput sequencing of the transcriptome and chromatin accessibility in the same cell. Nature biotechnology, 37(12):1452-1457 2019.</li>
                                      <li>Zhu, C., et al. An ultra high-throughput method for single-cell joint analysis of open chromatin and transcriptome. Nature Structural and Molecular Biology, 2019.</li>
                                      <li>Spektor, R., et al. Single cell atac-seq identifies broad changes in neuronal abundance and chromatin accessibility in down syndrome. bioRxiv, 2019.</li>
                                      </li></ul>')
                            )
                        )
                    ),
                      
                    box(
                        status='primary', width = NULL,
                        solidHeader=F, 
                        column(width=8,
                            p(
                                h4(strong("History")),
                                HTML('2022.05.31 First release.'),
                                br(),
                                HTML('2021.10.21 New data released.'),
                                br(),
                                HTML('2021.10.18 Pre-release.')
                            )
                        )
                    )
                )
            )
        ),
        tabItem(
            tabName = 'celltype',
            fluidPage(
                tags$head(
                    tags$style(
                        HTML("
                        p {
                            line-height: 170%;
                        }"),
                    )
                ),
                titlePanel(
                    div('Cross dataset analysis for each cell type', align="center")
                ),
                box(
                    status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE, collapsed=FALSE, height=320,
                    column(width=8, #style="h6 {line-height: 5000%;}",
                        p( h4(strong("Integrative annotation using meta-analytic marker genes")),
                           HTML("<ul><li>Due to the limitation of technology and resource, there is a <b>variation in the granularity of each scATAC-seq data granularity (9 to 36 clusters)</b>.</li>
                                <li>The annotated cell types include major cell types and more detailed cell types, such as Pvalb or Vip within inhibitory neuronal cells.</li>
                                <li>To overcome the heterogeneity of scATAC-seq datasets, we <b>reannotated all scATAC-seq data for detailed cell types </b>at each cell level rather than clusters.</li>
                                <li>We compute <b>module activity of meta-analytic marker gene sets </b> to capture cell-type specific signal at single-cell level.</li>
                                <li>The top 500 cells are extracted from each dataset according to the enrichment of module activity for each cell type to compute meta scATAC-seq profiles.</li></ul>")
                        )
                    )
                ),
                fluidRow(
                    column(width = 8,
                        textInput("gene", 
                            label = "Gene symbol", value = "Gad2"
                        ),
                        textInput("region", 
                            label = "Genomic region (activated only when a gene name is blank)", value = "chr2:22622205-22693874"
                        ),
                        sliderInput("offset", "Flanking region (+/-)", value = 10000, min = 0, max = 100000, step=5000)
                    ),
                    column(width = 8,
                        selectInput(
                            inputId="dtrack.type", selected='line', 
                            label='Plot type: ',
                            choices=c('line', 'gradient', 'horizon')
                        ),
                        selectInput(
                            inputId="ctmarker", selected='EX_L5.ET',
                            label='Cell type: ',
                            choices=celltype.list
                        ),
                        actionButton("go", "Run"),
                        br(),br()
                    ),
                    column(width = 8,
                        withSpinner(
                            plotOutput('ctrack', height='600px', width='100%')
                        , type=4, color = "#0dc5c1")
                    )
                ),
                box(title = 'Surrounding gene information', status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE, collapsed=FALSE, 
                    fluidRow(
                        column(width = 8,
                            DT::dataTableOutput('overlappedGenes')
                        )
                    )                    
                )
            )
        ),
        tabItem(
            tabName = 'meta',
            fluidPage(tweaks,
                titlePanel(
                    div('Cross cell-type analysis within each dataset', align="center")
                ),
                box(
                    status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE, collapsed=FALSE, height=200,
                    column(width=8, 
                        p(h4(strong("Cell-type specific meta profiles within each dataset")),
                           HTML("<ul><li>Cluster-specific pseudo-bulk profiles are computed for each dataset.</li>
                                <li>Moreover, predicted scATAC-seq profiles from DNA sequence only are shown as a representative of sequence-dependent chromatin accessibility.</li>
                                <li>The model for prediction is trained for all clusters found in the BICCN datasets.</li></ul>")
                        )
                    )
                ),
                fluidRow(
                    column(width = 12,
                        textInput("cgene", 
                            label = "Gene symbol", value = "Sst"
                        ),
                        textInput("cregion", 
                            label = "Genomic region (activated only when a gene name is blank)", value = "chr16:23889573-23890958"
                        ),
                        sliderInput("coffset", "Flanking region (+/-)", value = 10000, min = 0, max = 100000, step=5000)
                    ),
                    column(width = 8,
                        selectInput(
                            inputId="cdtrack.type", selected='line', 
                            label='Plot type: ',
                            choices=c('line', 'gradient', 'horizon')
                        ),
                        selectInput(
                            inputId="study", selected=dataset.list[1],
                            label='Dataset: ',
                            choices=dataset.list[-c(5)]
                        ),
                        br(),
                        tags$div(align='left', class = 'multicol', 
                            withSpinner(uiOutput("clusterSelection"), type=1, color = "#0dc5c1", size=0.5)
                        ),
                        checkboxInput("unselectAll", "Unselect all checkboxes", FALSE),
                        checkboxInput("selectAll", "Select all checkboxes", FALSE),
                        actionButton("cgo", "Run"),
                        br(),br()
                    ),
                    column(width = 8,
                        withSpinner(uiOutput("ui"), type=4, color = "#0dc5c1")
                    )
                ),
                box(title = 'Overlapped genes', status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE, collapsed=FALSE, 
                    fluidRow(
                        column(width = 8,
                            DT::dataTableOutput('overlappedGenesCluster')
                        ),
                        column(width=4,
                            uiOutput('downloadflag'),
                            uiOutput('option'),
                            uiOutput('connectionNumber'),
                            uiOutput('downloadButton')  
                        )
                    )                    
                )
            )
        ),
        tabItem(
            tabName = 'single',
            fluidPage(
                titlePanel(
                    div('Reannotation of original clusters', align="center")
                ),
                box( status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE, collapsed=FALSE, height=200,
                    fluidRow(
                        column( width=8,
                            p(
                                h4(strong('Module activity of meta-analytic marker genes for clusters')),
                                h5('Cell-type characterization of each cluster based on overlap of cluster-specific genes and marker gene set.'),
                            ),
                            p(style="margin-left: 40px", 
                                h5(strong('1. Major cell type:'), 'IN (inhibitory), EX (excitatory), and NN (non neuronal)'),
                                h5(strong('2\' Neuronal subtype (IN)'), 'Lamp5, Pvalb, Sncg, Sst, and Vip'),
                                h5(strong('2\'\' Neuronal subtype (EX): '), 'L2.3 IT, L5.IT, L5.ET, L5.6.NP, L6.IT.Car3, L6.CT, L6.b'),
                                h5('Other conditions - ', 'score: normalized Jaccard index, # of cluster-specific genes: 1000, # of marker genes: 100.')
                            )
                        )
                    )
                ),
                fluidRow(
                    column( width=8,
                        selectInput(
                            inputId="hm.data", selected=dataset.list[1], selectize=F,
                            label='Dataset: ',
                            choices=dataset.list[-c(7,8)]
                        ),
                        selectInput(
                            inputId="hm.marker", selected='Major celltype', selectize=F,
                            label='Marker resolution: ',
                            choices=c('Major celltype', 'Neuronal subtype')
                        ),
                    ),
                    column(width = 8,
                        withSpinner(plotOutput('heatmap', height='600px'), type=4, color = "#0dc5c1")
                    )
                )
            )
        ),
        tabItem(
            tabName = 'download',
            fluidRow(
                box(
                    status='primary', width = NULL,
                    solidHeader=F, title='Marker gene sets',
                    column(width=12,
                        p(
                            HTML('<a href="https://github.com/carushi/Catactor/tree/main/marker_genes">https://github.com/carushi/Catactor/tree/main/marker_genes</a>')
                        )
                    ),
                ),
                box(
                    status='primary', width = NULL,
                    solidHeader=F, title='Pseudo-bulk profiles of scATAC-seq data',
                    column(width=12,
                        p(
                            HTML('<a href="http://labshare.cshl.edu/shares/gillislab/resource/Meta_scATAC">http://labshare.cshl.edu/shares/gillislab/resource/Meta_scATAC</a>')
                        )
                    ),
                ),
            )
        ),
        tabItem(
            tabName = 'link',
            fluidPage(
                box(
                    status='primary', width = NULL,
                    solidHeader=F, title='Meta-analysis of scATAC-seq cell-type characterization',
                    column(width=12,
                        p(
                            HTML('<a href="ttps://www.biorxiv.org/content/10.1101/2021.04.01.438068v1">Exploiting marker genes for robust classification and characterization of single-cell chromatin accessibility</a></br>
                            Risa Karakida Kawaguchi, Ziqi Tang, Stephan Fischer, Rohit Tripathy, Peter K. Koo, Jesse Gillis. bioRxiv 2021.04.01.438068; </br>doi: https://doi.org/10.1101/2021.04.01.438068')
                        )
                    ),
                ),
                box(
                    status='primary', width = NULL,
                    solidHeader=F, title='Deep learning model',
                    column(width=12,
                        p(
                            HTML('<a href="https://github.com/crajesh6/biccn">https://github.com/crajesh6/biccn</a>')
                        )
                    )
                ),
                box(
                    status='primary', width = NULL,
                    solidHeader=F, title='Marker set construction',
                    column(width=12,
                        p(
                            HTML('<a href="https://www.biorxiv.org/content/10.1101/2021.04.16.439807v2">How many markers are needed to robustly determine a cell’s type?</a></br>
                            Stephan Fischer, Jesse Gillis. bioRxiv 2021.04.16.439807; </br>doi: https://doi.org/10.1101/2021.04.16.439807')
                        )
                    ),
                ),
                box(
                    status='primary', width = NULL,
                    solidHeader=F, title='Lab website',
                    column(width=12,
                        p(
                            HTML('<a href="http://gillislab.labsites.cshl.edu/">Gillis lab website</a>')
                        )
                    )
                )

            )
        )
    )
)



sidebar <-  dashboardSidebar(
    sidebarMenu(
		menuItem("Top", tabName = "main", icon = icon("compass")),
		menuItem("Cross cell-type analysis", tabName = "meta", icon = icon("brain")),
		menuItem("Cross dataset analysis", tabName = "celltype", icon = icon("solar-panel")),
		menuItem("Reannotation of clusters", tabName = "single", icon = icon("book")),
		menuItem("Download", tabName = "download", icon = icon("download")),
		menuItem("Link", tabName = "link", icon = icon("external-link-alt"))
    )
)

header = dashboardHeaderPlus(
    title = div('Meta scATAC-seq analysis for mouse brain @ Gillis lab'),
    titleWidth=900
)

dashboardPage(
    header, sidebar, body, skin='black'
)
