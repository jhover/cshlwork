
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#   
# http://shiny.rstudio.com
#
# header = dashboardHeaderPlus(title = 'CoCoCoNet',titleWidth=150)

function(request){
    header = dashboardHeaderPlus( 
        title = 'FAIRit',
        titleWidth=150)


    sidebar = dashboardSidebar(
        width = 150,
        sidebarMenu(
            id = 'tabs', 
            menuItem('Query',tabName = 'query',icon = icon('search')),
            menuItem('Project',tabName = 'project',icon = icon('database')),
            menuItem('Reference',tabName = 'tab2',icon = icon('sitemap')),
        
            menuItem("Download",icon = icon("download"),href='http://labshare.cshl.edu/shares/gillislab/resource/MetaQC'),
            uiOutput('bookmarkurl')
        )

    )


    ##################################################################################################*
    ######################################### Start Main body #########################################
    ##################################################################################################*

    body = dashboardBody( 
        # dashboardthemes::shinyDashboardThemes(
        #     theme = "spacelab"
        # ),

        tabItems(
            tabItem(tabName = 'query',
                fluidPage(theme =shinythemes::shinytheme("spacelab"),
                    box(    width=NULL, title = 'FAIRit-BRAIN',
                        fluidRow(
                            column(width = 9,
                                "FAIRit aims to make single-cell RNAseq data more FAIR (Findable Accessible Interoperable Reusable).",br(),
                                "Filter datasets (projects) based on desired constraints and quality control metrics.",br(),
                                "Download metadata for resultant subset.",br(),
                                "Click on project id/UMAP to go to project-specific page with additional statistics and plotting."
                            ),
                            column(width = 3,
                                # p('Download section - metadata'),
                                
                                # uiOutput('downloadflag'),
                                uiOutput('downloadoptions'),
                                # uiOutput('connectionNumber'),
                                uiOutput('downloadButton')  ,
                                br()
                            )
                        ),
                    
                        box(width =NULL,
                            column(width = 2,
                                h4('Technology:'),
                                checkboxInput('is10x' , 
                                    label = 'Chromium 10x',
                                    value = TRUE
                                ),
                                checkboxInput('isSS' , 
                                    label = 'SmartSeq',
                                    value = TRUE
                                ),
                                h4('Data Source:'),
                                
                                checkboxInput('isSRA' , 
                                    label = 'Sequence Read Archive (SRA)',
                                    value = TRUE
                                ),
                                
                                checkboxInput('isNemo' , 
                                    label = 'NeMO',
                                    value = TRUE
                                ),

                                h4('Most prevalent class:'),
                                pickerInput('primary_class',
                                    # label ="Predominant class label:",
                                    choices = unique(SEARCH_DF$primary_class),
                                    selected = unique(SEARCH_DF$primary_class),
                                    multiple=TRUE,
                                    options = pickerOptions(
                                        virtualScroll=10,
                                        liveSearch=TRUE,
                                        size=10,
                                        actionsBox = TRUE
                                    ),
                                    width='100%',
                                    inline=FALSE

                                ),
                                

                                
                                h4('Most prevalent subclass:'),
                                pickerInput('primary_subclass',
                                    # label ="Predominant class label:",
                                    choices = unique(SEARCH_DF$primary_subclass),
                                    selected = unique(SEARCH_DF$primary_subclass),
                                    multiple=TRUE,
                                    options = pickerOptions(
                                        virtualScroll=10,
                                        liveSearch=TRUE,
                                        size=10,
                                        actionsBox = TRUE
                                    ),
                                    width='100%',
                                    inline=FALSE

                                ),

                                chooseSliderSkin("Square", color = 'grey'),
                                h4('Publish Date:'),
                                shinyWidgets::sliderTextInput('date_slider',label=NULL,
                                    choices =zoo::as.yearmon(seq(from = min(SEARCH_DF$pdate -5,na.rm =TRUE),to =max(SEARCH_DF$pdate+28,na.rm=TRUE),by='mon')),
                                    selected=c(min(zoo::as.yearmon(SEARCH_DF$pdate),na.rm=TRUE),  max(zoo::as.yearmon(SEARCH_DF$pdate),na.rm=TRUE) ),
                                    # grid=TRUE,
                                    width='100%'
                                ),
                                h4('Number of Cells'),

                                shinyWidgets::sliderTextInput('ncell_slider',label=NULL,
                                    choices = format(2^seq(from = floor(log2(min(SEARCH_DF$ncell))), to =ceiling(log2(max(SEARCH_DF$ncell))) , by=1) ,big.mark =','),
                                    selected = format(c(2^floor(log2(min(SEARCH_DF$ncell))) , 2^ceiling(log2(max(SEARCH_DF$ncell)))),big.mark =','),
                                    grid = TRUE,
                                    width='100%'
                                ),
                                

                                # sliderInput('ncell_slider',label=NULL,
                                #     min = floor(min(log(SEARCH_DF$ncell))*10 )/10, 
                                #     max = ceiling(max(log(SEARCH_DF$ncell))*10 )/10,
                                #     value=c(floor(min(log(SEARCH_DF$ncell))*10 )/10 , 
                                #             ceiling(max(log(SEARCH_DF$ncell))*10 )/10),
                                #     step = .25,
                                #     width='100%'
                                # ),
                                
                                                                
                                h4('Gini Coefficient'),
                                sliderInput('gini_ds_slider',label=NULL,
                                    min = floor(min(SEARCH_DF$gini_ds)*10 )/10, 
                                    max = ceiling(max(SEARCH_DF$gini_ds)*10 )/10,
                                    value=c(floor(min(SEARCH_DF$gini_ds)*10 )/10 , 
                                            ceiling(max(SEARCH_DF$gini_ds)*10 )/10),
                                    step = .01,
                                    width='100%'
                                ),
                                
                                shinyBS::bsTooltip(                      
                                    id="gini_ds_slider", 
                                    title= "The Gini coeffiecient is a measure of disparity in total counts for each dataset.",
                                    options = list(container = "body"),
                                    placement='left'
                                ),  
                                
                                h4('Class Precision Recall'),
                                sliderInput('class_pr_slider',label=NULL,
                                    min = floor(min(SEARCH_DF$class_pr)*10 )/10, 
                                    max = ceiling(max(SEARCH_DF$class_pr)*10 )/10,
                                    value=c(floor(min(SEARCH_DF$class_pr)*10 )/10 , 
                                            ceiling(max(SEARCH_DF$class_pr)*10 )/10),
                                    step = .01,
                                    width='100%'
                                ),
                                
                                shinyBS::bsTooltip(                      
                                    id="class_pr_slider", 
                                    title= "The Precision Recall is a measure of the quality of Class labels predicted by MetaMarkers.",
                                    options = list(container = "body"),
                                    placement='left'
                                ),  
                                
                                h4('Subclass Precision Recall'),
                                sliderInput('subclass_pr_slider',label=NULL,
                                    min = floor(min(SEARCH_DF$subclass_pr_hier)*10 )/10, 
                                    max = ceiling(max(SEARCH_DF$subclass_pr_hier)*10 )/10,
                                    value=c(floor(min(SEARCH_DF$subclass_pr_hier)*10 )/10 , 
                                            ceiling(max(SEARCH_DF$subclass_pr_hier)*10 )/10),
                                    step = .01,
                                    width='100%'
                                ),
                                shinyBS::bsTooltip(                      
                                    id="subclass_pr_slider", 
                                    title= "The Precision Recall is a measure of the quality of Subclass labels predicted by MetaMarkers.",
                                    options = list(container = "body"),
                                    placement='left'
                                ),  
                                h4('Percent Dropout'),
                                sliderInput('percent_dropout_slider',label=NULL,
                                    min = floor(min(SEARCH_DF$percent_dropout) ), 
                                    max = ceiling(max(SEARCH_DF$percent_dropout)),
                                    value=c(floor(min(SEARCH_DF$percent_dropout)), 
                                            ceiling(max(SEARCH_DF$percent_dropout)) ),
                                    step = 1,
                                    width='100%'
                                ),                   
                                
                                shinyBS::bsTooltip(                      
                                    id="percent_dropout_slider", 
                                    title= "Dropouts are defined as genes with zero expression within the dataset.",
                                    options = list(container = "body"),
                                    placement='left'
                                ),               
                                h4('Mitochondria AUROC'),
                                sliderInput('mt_slider',label=NULL,
                                    min = floor(min(SEARCH_DF$mt)*10 )/10, 
                                    max = ceiling(max(SEARCH_DF$mt)*10 )/10,
                                    value=c(floor(min(SEARCH_DF$mt)*10 )/10 , 
                                            ceiling(max(SEARCH_DF$mt)*10 )/10),
                                    step = .01,
                                    width='100%'
                                ),                   

                                shinyBS::bsTooltip(                      
                                    id="mt_slider", 
                                    title= "Mitochondria AUROCs are a measure of dissimilarity between the mitochondria counts in the dataset to the global distribution.",
                                    options = list(container = "body"),
                                    placement='left'
                                ),                          
                            ),
                                
                            
                            column(width = 10,
                                div(style = 'height:1100px; overflow-y: scroll;',
                                    DT::dataTableOutput('dataset_search')
                                )       
                            )

                        )
                    )
                )
            ),

            tabItem( tabName = "project",
                # Helper text goes here
                fluidPage(   theme =shinythemes::shinytheme("spacelab"),
                    # titlePanel(),
                
                
                    # main tab with figures
                    ### Main info tab 
                    box(title = textOutput('title') ,width = NULL,
                        column(width=7,
                            
                                box( width = NULL, 
                                    div(style = "height:125px;",
                                    fluidRow(
                                        column( width=2,
                                            h4('Identifiers:')  
                                        ),
                                        
                                        column( width=3,
                                            p('SRA:') 
                                        ),
                                        column(width=4,
                                            pickerInput(                   
                                                inputId="projid",
                                                label=NULL, multiple=FALSE,
                                                choices=AVAILABLE_PROJECTS,
                                                selected= 'SRP179101' ,
                                                width='fit',
                                                options = pickerOptions(
                                                    virtualScroll=10,
                                                    liveSearch=TRUE,
                                                    size=10
                                                ),
                                                inline=FALSE
                                            )  
                                        ),
                                        column(width=3,
                                                                                                
                                        )
                                    ),
                                    
                                    fluidRow(
                                        column(width=2),
                                        column(width=3,
                                            htmlOutput('ext_ids_sub')
                                            # p('GEO:')
                                        ),
                                        column(width=6, align ='left',
                                            htmlOutput('ext_ids')
                                            # p('GEO accession')
                                        )
                                    ),
                                    fluidRow(
                                        column(width=2),
                                        column(width=10,
                                            br(),
                                            htmlOutput('labshare_url')
                                            
                                        )
                                        
                                    )                                    

                                )
                            ),
                            box(width=NULL,
                                fluidRow(
                                    column( width=2,

                                        h4('Abstract:')  
                                    ),
                                    column(width =10,
                                        div(style = "height:403px;overflow-y: scroll",
                                        textOutput('abstract')
                                        )

                                        # p('Some text for the abstract')
                                    )
                                
                                )
                            )
                        ),
                        column( width = 5,
                            box(width = NULL, 
                                fluidRow(
                                    column(width=12,
                                        
                                        h4('Information:'),
                                        div(style = 'height:150px; overflow-y: scroll', 
                                            tableOutput('basic_info') 
                                        )
                                    )
                                )
                            ),
                            box(width=NULL,
                                fluidRow(
                                    column(width=12,
                                        h4('Sample Attributes:'),
                                        div(style = 'height:300px; overflow-y: scroll',
                                            DT::dataTableOutput('sample_attr') #overflow-x: scroll;
                                        )   
                                    )
                                )
                                # DT::dataTableOutput("",
                                #     style = ";overflow-x: scroll;")
                                # p('Table of sample attributes.')
                            )
                        )

                    ),

                    ### Main figures tab 
                    # idea - allow users to flag cells based on figures 
                    box( title = 'Quality Control statistics',width = NULL,
                        fluidRow( 
                            column(width = 8, 
                                pickerInput(                    # select cell specific stats 
                                    inputId="obsname",
                                    label="Statistic:", multiple=FALSE,
                                    choices=CTS_LABS,
                                    width='fit',
                                    options = pickerOptions(
                                        virtualScroll=20,
                                        liveSearch=FALSE,
                                        size=10
                                    ),
                                    inline=TRUE
                                ),br(),                            
                            # ),
                            # column(width = 2, 

                                # deprecated
                                # pickerInput(                    # select 
                                #     inputId="obsunit",
                                #     label="Metric:", multiple=FALSE,
                                #     choices=OBSUNITS, 
                                #     selected = OBSUNITS[1],
                                #     width='fit',
                                #     options = pickerOptions(
                                #         liveSearch=FALSE
                                #     ),
                                #     inline=TRUE
                                # )                            
                            )
                        ),
                        fluidRow(
                            column(width = 5,
                                plotOutput('umap',
                                    dblclick = "umap_dblclick",
                                    brush = brushOpts(
                                        id = "umap_brush",
                                        resetOnNew = TRUE
                                    )                                
                                ),
                                
                                shinyBS::bsTooltip(                      # info box example
                                    id="umap", options = list(container = "body"),placement='right',
                                    title=paste("Click and hold to brush and select a region. Double click inside the region to zoom in, or double click to zoom out.")
                                ),
                            ),
                            column(width = 5,
                                plotOutput('hist')
                            ),
                            column(width = 2,
                                htmlOutput( 'stat1_info')

                                    # plotOutput('umap_zoom',
                                    #     brush = brushOpts(
                                    #         id = "plot1_brush",
                                    #         resetOnNew = TRUE
                                    #     )
                                    # )
                            )                            
                        ),
                        fluidRow(
                            column(width = 8, 
                                pickerInput(                    # select cell specific stats 
                                    inputId="obsname2",
                                    label="Statistic 2:", multiple=FALSE,
                                    choices=CTS_LABS,
                                    width='fit',
                                    options = pickerOptions(
                                        virtualScroll=20,
                                        liveSearch=FALSE,
                                        size=10
                                    ),
                                    inline=TRUE
                                ),
                                br()
                            ),
                        ),
                        fluidRow(
                            column(width=5,
                                plotOutput('hist2' )
                            ),
                            column(width=5,
                                plotOutput('stat_comp')
                            ),
                            column(width=2,
                                
                            )
                        ),
                        br(),
                        fluidRow(
                            column(width = 1),
                            column(width = 8,
                            
                                plotOutput('pca_expvar',
                                    # dblclick = "pca_expvar_dblclick",
                                    # brush = brushOpts(
                                        # id = "pca_expvar_brush",
                                        # resetOnNew = TRUE,
                                    # )  ,
                                    height = '200px'
                                ),
                                helpText('Note: PCA without zero-centering was used: the explained variance does not correspond to the exact statistical defintion.
                                    The first component, e.g., might be heavily influenced by different means.
                                    The following components often resemble the exact PCA very closely.')
                            ),
                            column(width = 1)
                        )
                    ),


                    box( title = 'MetaMarkers: cell type annotation',width = NULL,
                        fluidRow( 
                            column(width = 4, 
                                
                                uiOutput('class_label_box'), br(),
                                uiOutput('class_mm_metric_box')
                        
                            ),
                            column(width = 4, 
                                uiOutput('subclass_label_box'), br(),
                                uiOutput('subclass_mm_metric_box')
                            )
                        ),
                        fluidRow(
                            column(width = 4,
                                plotOutput('class_umap_mm',
                                    dblclick = "class_umap_mm_dblclick",
                                    brush = brushOpts(
                                        id = "class_umap_mm_brush",
                                        resetOnNew = TRUE
                                    )                                
                                ),
                                
                                shinyBS::bsTooltip(                      # info box example
                                    id="class_umap_mm", options = list(container = "body"),placement='right',
                                    title=paste("Click and hold to brush and select a region. Double click inside the region to zoom in, or double click to zoom out.")
                                ),
                            ),
                            column(width = 4,
                                plotOutput('subclass_umap_mm',
                                    dblclick = "subclass_umap_mm_dblclick",
                                    brush = brushOpts(
                                        id = "subclass_umap_mm_brush",
                                        resetOnNew = TRUE
                                    )                                
                                ),
                                
                                shinyBS::bsTooltip(                      # info box example
                                    id="subclass_umap_mm", options = list(container = "body"),placement='right',
                                    title=paste("Click and hold to brush and select a region. Double click inside the region to zoom in, or double click to zoom out.")
                                ),
                            ),                            
                            column(width = 4,
                                plotOutput('class_mm_pr',
                                    dblclick = "class_umap_mm_pr_dblclick",
                                    brush = brushOpts(
                                        id = "class_umap_mm_pr_brush",
                                        resetOnNew = TRUE
                                    )   
                                )
                            )                         
                        ),
                        # fluidRow(

                        #     column(width = 5,
                        #         plotOutput('subclass_mm_pr',
                        #             dblclick = "subclass_umap_mm_pr_dblclick",
                        #             brush = brushOpts(
                        #                 id = "subclass_umap_mm_pr_brush",
                        #                 resetOnNew = TRUE
                        #             )   
                        #         )
                        #     )                   
                        # ),
                    )
                )
            ), # end tab item (tab_1) 

            tabItem(tabName = 'tab2',
                fluidPage(   theme =shinythemes::shinytheme("spacelab"),
                    # titlePanel(),
                
                
                    # main tab with figures
                    ### Main info tab 
                    box(width = NULL, 
                        h4('MetaNeighbor accesses replicability across datasets within the BICCN'),
                        box(width=NULL,
                            

                            p(
                                'The Brain Initiative - Cell Census Network (BICCN) is a ... 
                                Here we take high quality scRNA seq data from the BICCN, 
                                and assess replicability of cell types across datasets using MetaNeighbor (cite).
                                MetaNeighbor uses a neighbor voting algorithm... etc'
                            ),
                            
                            p(
                                'Something about the hierarchical cell type annotations'
                            ),
                            
                            pickerInput(                   
                                inputId="label_set",
                                label='View:', multiple=FALSE,
                                choices=c('Class labels'='class','Subclass labels'='subclass','Cluster labels'='cluster'),
                                selected='class' ,
                                width='fit',
                                options = pickerOptions(
                                    virtualScroll=10,
                                    liveSearch=FALSE,
                                    size=10
                                ),
                                inline=FALSE
                            
                            ),
                            fluidRow(
                                column(width=6,

                                    div(style = 'height:550px',
                                        plotOutput('biccn_MN_heatmap')
                                    )   
                                    
                                ),
                                column(width=6,
                                    
                                        DT::dataTableOutput('biccn_MN_top_hits')
                                )
                            ),
                            fluidRow(
                                column(width=6,
                                    # div(style = 'height:550px',
                                        visNetwork::visNetworkOutput('biccn_MN_network')
                                    # )
                                ),
                                column(width=6,
                                    div(style = 'height:425px;overflow-y: scroll',
                                        DT::dataTableOutput('biccn_MN_metaclust')
                                    )
                                )
                            )
                        ),

                        h4('MetaMarkers looks for specific and sensitive markers across the BICCN datasets'),
                        fluidRow(
                            
                            box(width=8,
                                p('some text'),
                                fluidRow(
                                    column(
                                        width = 6,
                                        pickerInput(                   
                                            inputId="mm_class_label",
                                            label='Class label:', multiple=FALSE,
                                            choices=sort(CLASS_COLORS$label[1:(length(CLASS_COLORS$label)-1)]),
                                            selected= sort(CLASS_COLORS$label)[1] ,
                                            width='fit',
                                            options = pickerOptions(
                                                virtualScroll=10,
                                                liveSearch=FALSE,
                                                size=10
                                            ),
                                            inline=FALSE
                                        )  ,
                                        plotOutput('pareto_front_class',
                                            dblclick = "class_par_dblclick",
                                            brush = brushOpts(
                                                id = "class_par_dblclick_brush",
                                                resetOnNew = TRUE
                                            ),
                                            hover = hoverOpts(
                                                id='class_par_hover',
                                                nullOutside = FALSE
                                            )                                                                                        
                                        ),
                                        htmlOutput("class_par_hover_near")

                                    ),
                                    column(
                                        width = 6,

                                        pickerInput(                   
                                            inputId="mm_subclass_label",
                                            label='Subclass label:', multiple=FALSE,
                                            choices=sort(SUBCLASS_COLORS$label[1:(length(SUBCLASS_COLORS$label)-2)]),
                                            selected=sort(SUBCLASS_COLORS$label)[1] ,
                                            width='fit',
                                            options = pickerOptions(
                                                virtualScroll=10,
                                                liveSearch=FALSE,
                                                size=10
                                            ),
                                            inline=FALSE
                                        )  ,                                        
                                        plotOutput('pareto_front_subclass',
                                            dblclick = "subclass_par_dblclick",
                                            brush = brushOpts(
                                                id = "subclass_par_dblclick_brush",
                                                resetOnNew = TRUE
                                            ),
                                            hover = hoverOpts(
                                                id='subclass_par_hover' ,
                                                nullOutside = FALSE
                                            )
                                        ),
                                        htmlOutput("subclass_par_hover_near")
                                        

                                    ),
                                    # column(
                                    #     width = 4,
                                    #     plotOutput('pareto_front_cluster')
                                    # )
                                )
                            ),
                            box(width=4,
                                p('some text'), 
                                DT::dataTableOutput('pareto_front_cluster')
                            )
                        )
                    )
                )
            )
                        

        )   # end tabitems
    )   # end dashboard body           

    # build everything
    dashboardPagePlus(skin = "black",header,sidebar,body,enable_preloader = TRUE, sidebar_fullCollapse = T,loading_duration = 0)



}


# between 100-500
# SRP066314
# SRP108034
# SRP228572



# broken projects
# SRP166780 - no markers
# SRP237263 - no markers
# SRP250815 - no markers
# SRP254109 - no markers
# SRP256315 - breaks - still working
# SRP269209 - no markers
# SRP273307 - breaks - still working
# SRP304784 - breaks - still working
