
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#   
# http://shiny.rstudio.com
#
# header = dashboardHeaderPlus(title = 'CoCoCoNet',titleWidth=150)
header = dashboardHeaderPlus(
    title = 'Mouse HSC Atlas',
    titleWidth=150)


sidebar = dashboardSidebar(
    width = 150,
    sidebarMenu(
        id = 'tabs', 
        menuItem("Home",tabName =  "home", icon = icon("home")),
        # menuItem("Comparison", tabName = "coexnet", icon = icon("project-diagram")), 

        menuItem('Cell Clusters',tabName = 'tab_1',icon = icon('sitemap')),
        
        menuItem('Cell States',tabName = 'tab_2',icon = icon('sitemap')),

        menuItem('Pseudotime',tabName = 'tab_3',icon = icon('sitemap')),
        
        menuItem('Cross Species',tabName = 'tab_4',icon = icon('sitemap')),

        menuItem("Data download",icon = icon("download"),href='http://labshare.cshl.edu/shares/gillislab/resource/HSC_atlas/')
        # menuItem("Help",icon = icon("question"),href = "https://github.com/johnlee4/CoCoCoNet")
    )
)


##################################################################################################*
######################################### Start Main body #########################################
##################################################################################################*

body = dashboardBody(   
    # shinyDashboardThemes(
    #     theme = "poor_mans_flatly"
    # ),

    tabItems(
        tabItem( tabName = 'home',
            fluidPage(
                titlePanel(
                    h1(div('Mouse HSC Atlas'), align = "center")
                )
            ),
            fluidRow(
                box(

                    status='primary',width = NULL,
                    solidHeader=F,
                    column(width=12,

                        p( h2("A Meta-Analytic Single-Cell Atlas of Mouse Bone Marrow Hematopoietic Development" ),
                            'Harris, Lee, and Gillis 2021'
                        ), br(),
                        HTML(' '),
                
                

                        br(),
                        fluidRow(
                            
                            column(width=3,
                                h4('Cell Cluster Analysis' ),
                                p('Clusters from Figures 1 and 2')  
                            ),
                            
                            column(width=3,
                                # actionButton(inputId='toConservation',label='Co-expression Conservation',icon=icon('arrow-right'))
                                h4('Cell State Analysis' ),
                                p('Cell states from Figure 3')  

                            ),
                            column(width= 3,
                                h4('Pseudotime Analysis' ) ,
                                p('Pseudotime analysis from Figure 4')  
                            ),
                                h4('Cross Species Analysis' ) ,
                            column(width= 3,
                                p('Cross species analysis from Figure 5')  
                            )
                        
                        ),
                        fluidRow(
                            column(width=3,
                                actionBttn(inputId='to_tab_1',label='Clusters',icon=icon('arrow-right') )
                            ),
                            column(width=3,
                                actionBttn(inputId='to_tab_2',label='Cell State',icon=icon('arrow-right'))
                            ),
                            column(width=3,
                                actionBttn(inputId='to_tab_3',label='Pseudotime',icon=icon('arrow-right'))
                            ),
                            column(width=3,
                                actionBttn(inputId='to_tab_4',label='Cross Species',icon=icon('arrow-right'))
                            ),
                        ),
                        br(),br(),
                    

                        fluidRow(
                        
                            p(

                                #  <a href='https://academic.oup.com/nar/article/48/W1/W566/5835815' target='_blank'> CoCoCoNet</a>
                                # column( width=7,
                                #     tags$h4('Descriptive Header 1'),
                                #         HTML("Some text..." 
                                #         ),br(),

                                #     tags$h4('Descriptive Header 2'),
                                #         HTML("Some text...")

                                    

                                

                                # ),
                                column( width=12,
                                    div(img(src="intro_cartoon.png",height='100%',width='100%'), align='center'   )                       
                                )
                            )

                        ),
                        

                        br(),
                        p(
                        tags$h4('Did you find this HSC Atlas helpful?'),
                            HTML('If so, please consider citing us at: <br/><br/>
                            Benjamin D. Harris, John Lee, and Jesse Gillis,<br/> <b>A Meta-Analytic Single-Cell Atlas of Mouse Bone Marrow Hematopoietic Development</b>, 
                            <i>BioRxiv</i>, August 2021, <br/><a href="https://doi.org/10.1101/2021.08.12.456098" target="_blank">https://doi.org/10.1101/2021.08.12.456098</a>'
                            )

                        ),

                        br(),
                        p(
                        tags$h4('Have a question or issue?'),
                            HTML('Feel free to contact John Lee at <b>johlee@cshl.edu</b> to report an issue, suggest an update, or ask for help. ')
                        

                        ),

                    ),
                    # column(width=4, 
                    #     div(img(src="cococonet_v3tree-01.png", height= "70%",width= "80%"), align='center'   )                       
                    # )
                

                )
            )
                
                
              
            
        ),

        tabItem( tabName = "tab_1",
            # Helper text goes here
            fluidPage(   
                titlePanel('Cell Cluster Analysis')
            ),

            # main tab with figures
            fluidRow(
                box(status='primary',width = NULL, solidHeader=F,
                    h2('Use'),
                    h3('Genes: '),
                    'The cluster analysis allows you to look at the expression of gene and GO terms in labled clusters',
                    'If you have a gene you are interested in you can see the expression of it in the gene on the UMAP',
                    'You can also look at the the differential expresssion summary statistics using the table exploration',
                    'An important statistic to look at is auroc, a value >.9 tells you that this gene is a very good marker for the given cluster',
                    'An AUROC=.5 means that the gene does no better than random. Statistics like fold_change and fold_change_detection serve as the signal-to-noise ratio or effect size.',
                    br(),
                    br(),
                    br(),

                    
                    # slow! render empty list first then update in server.
                    # row of gene input
                    fluidRow(  
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_1_gene",
                                label='Consider this gene...', multiple=FALSE,
                                choices=NULL
                            )
                        )
                    ), 
                    # top panels
                    fluidRow( 
                        column( width = 6, 
                             plotOutput('tab1_umap', height='600px')#, width="600px")
                        ), 
                        column( width = 6,
                            fluidRow(
                                plotOutput('tab1_gene_violin' , height='300px' )#, width="600px") 
                            ),
                            fluidRow(
                                DT::dataTableOutput('tab1_marker_table' , height='300px' )#, width="600px") 
                            )
                        )
                    ),
                    
                    # row of GO term input
                    # slow! render empty list first then update in server.
                    h3('Functional Annotation'),
                    'For looking at functional annotations, the AUROC of the GO term selected will appear as a Star on the violin plots showing the distributino for all GO terms',
                    'The AUROC is showing how well a given cluster can be identified based on expression of the genes in that GO term, using MetaNeighbor',
                    'The Dotplot shows the normalized (across datasets) expression of the genes that are members of the GO term selected.',
                    br(),
                    br(),
                    fluidRow(  
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_1_go",
                                label='Consider this GO term...',
                                choices=NULL      # change this to go terms
                            )   
                        )
                    ), 

                    # bottom panels
                    fluidRow( 
                        column(width = 6, 
                            plotOutput('tab1_go_enrich_violin' , height = '600px')
                        ), 
                        column(width = 6, 
                            # plotly::plotlyOutput('tab1_go_dot' , height = '600px')
                            plotly::plotlyOutput('tab1_go_dot' , height = '600px')

                        )
                    )
                )

            )
            


        ) , # end tab item (tab_1) 

        tabItem( tabName = "tab_2",
            # Helper text goes here
            fluidPage(   
                titlePanel('Cell State Analysis')
            ),

            # main tab with figures
            fluidRow(
                box(status='primary',width = NULL, solidHeader=F,
                   h2('Use'),
                    h3('Genes: '),
                    'The cluster analysis allows you to look at the expression of gene and GO terms in labled clusters',
                    'If you have a gene you are interested in you can see the expression of it in the gene on the UMAP',
                    'You can also look at the the differential expresssion summary statistics using the table exploration',
                    'An important statistic to look at is auroc, a value >.9 tells you that this gene is a very good marker for the given cluster',
                    'An AUROC=.5 means that the gene does no better than random. Statistics like fold_change and fold_change_detection serve as the signal-to-noise ratio or effect size.',
                    br(),
                    br(),
                    br(),
                    
                    # slow! render empty list first then update in server.
                    # row of gene input
                    
                    fluidRow(  
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_2_gene",
                                label='Consider this gene...', multiple=FALSE,
                                choices=NULL
                            )
                        )
                    ), 
                    # top panels
                    fluidRow( 
                        column( width = 6, 
                            plotOutput('tab2_umap', height='600px')#, width="600px")
                        ), 
                        column( width = 6,
                            fluidRow(
                                plotOutput('tab2_gene_violin' , height='300px' )#, width="600px") 
                            ),
                            fluidRow(
                                DT::dataTableOutput('tab2_marker_table' , height='300px' )#, width="600px") 
                            )
                        )
                    ),
                    
                    # row of GO term input
                    # slow! render empty list first then update in server.
                     h3('Functional Annotation'),
                    'For looking at functional annotations, the AUROC of the GO term selected will appear as a Star on the violin plots showing the distributino for all GO terms',
                    'The AUROC is showing how well a given cluster can be identified based on expression of the genes in that GO term, using MetaNeighbor',
                    'The Dotplot shows the normalized (across datasets) expression of the genes that are members of the GO term selected.',
                    br(),
                    br(),
                    fluidRow(  
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_2_go",
                                label='Consider this GO term...',
                                choices=NULL      # change this to go terms in server
                            )   
                        )
                    ), 

                    # bottom panels
                    fluidRow( 
                        column(width = 6, 
                           plotOutput('tab2_go_enrich_violin' , height = '600px')
                        ), 
                        column(width = 6, 
                            plotly::plotlyOutput('tab2_go_dot' , height = '600px')
                            # plotOutput('tab2_go_dot' , height = '600px')

                        )
                    )
                )

            )
        ) ,  # end tab item (tab_2)

        tabItem( tabName = "tab_3",
            # Helper text goes here
            fluidPage(   
                titlePanel('Pseudotime Analysis')
            ),

            # main tab with figures
            fluidRow(
                box(status='primary',width = NULL, solidHeader=F,
                    h2('Use: '),
                    'The left plots show each cell, within each dataset, colored by the pseudotime. On the right side you can select a gene', 
                    'To see the expression of it in each individual dataset select a gene from the dropdown menu',
                    'You can search through the table on the bottom left for your gene of interest to see how well is associated with pseudotime',
                    'fc is the fold change (signal to noise ratio/ effect size)',
                    'p_adj is the P value for testing whether a gene is associated with psuedotime ordering, irrespective of lineage',
                    'while p_adj_eryth and p_adj_mono are p values from testing association with pseudotime + lineage',
                    br(),
                    
                    # slow! render empty list first then update in server.
                    # row of gene input
                    fluidRow(  
                        column(width=6),
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_3_gene",
                                label='Consider this gene...', multiple=FALSE,
                                choices=NULL
                            )
                        )
                    ), 
                    # top panels
                    fluidRow( 
                        column( width = 6, 
                            div(img(src="pseudotime_umap_pseudotime_color.png",height='600px',width="100%"), align='center'   )                       
                            # plotOutput('tab3_umap', height='600px')#, width="600px")
                        ), 
                        column( width = 6,
                            fluidRow(
                                plotOutput('tab3_umaps' , height = '600px')
                            )
                        )
                    ),
                    

                    # bottom panels
                    fluidRow( 
                        column(width = 6, 
                            DT::dataTableOutput('tab3_pt_stats' , height='600px' )#, width="600px") 

                        ), 
                        column(width = 6, 
                            DT::dataTableOutput('tab3_go_enrich' , height = '600px')
                        )
                    )
                )

            )
        ),   # end tab item (tab_3)


        tabItem( tabName = "tab_4",
            # Helper text goes here
            fluidPage(   
                titlePanel('Cross Species')
            ),

            # main tab with figures
            fluidRow(
                box(status='primary',width = NULL, solidHeader=F,
                    h2('Use: '),
                    'Here we have all of the cells from the 3 species we used, Mouse, Human, and Zebrafish.',
                    'You can select a gene symbol from one of the species and it will display the expression of that gene and',
                    ' 1:1 orthologs if available for the other two species. ',
                    br(),
                    br(),
                    
                    # row of gene input
                    fluidRow(  
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_4_mouse_gene",
                                label='Consider this gene in Mouse', multiple=FALSE,
                                choices=NULL
                            )
                        ),
                        column( width=1),
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_4_human_gene",
                                label='with this ortholog in Human', multiple=FALSE,
                                choices=NULL
                            )
                        ),
                        column( width=1),
                        column( width=3,
                            selectizeInput(                    # select species
                                inputId="tab_4_zebrafish_gene",
                                label='and this ortholog in Zebrafish', multiple=FALSE,
                                choices=NULL
                            )
                        )

                    ), 
                    # top panels
                    fluidRow(
                        column( width = 4, 
                            div(img(src="mouse_umap_cell_type.png",height='100%',width='100%'), align='center'   )
                        ), 
                        column( width = 4, 
                            div(img(src="human_umap_facs_labels.png",height='100%',width='100%'), align='center')
                        ), 
                        column( width = 4, 
                            div(img(src="zebrafish_estimated_cell_type.png",height='100%',width='100%'), align='center')
                        ),                                                                         
                    ),
                    # bottom panels
                    fluidRow( 
                        column( width = 4, 
                            plotOutput('tab4_mouse_umap', height='400px')
                        ), 
                        column( width = 4, 
                            plotOutput('tab4_human_umap', height='400px')
                        ), 
                        column( width = 4, 
                            plotOutput('tab4_zebrafish_umap', height='400px')
                        ), 
                    ),
                )
            )  
        )  # end tab item (cons) 



    )   # end tabitems
)   # end dashboard body           


# build everything
dashboardPagePlus(header,sidebar,body,enable_preloader = TRUE, sidebar_fullCollapse = T,loading_duration = 0)



