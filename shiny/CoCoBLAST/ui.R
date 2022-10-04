
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#  
# http://shiny.rstudio.com
#
# header = dashboardHeaderPlus(title = 'CoCoCoNet',titleWidth=150)
header = dashboardHeaderPlus(
    title = div('CoCoBLAST'),
    titleWidth=150)


sidebar = dashboardSidebar(
    width = 150,
    sidebarMenu(
        id = 'tabs',
        menuItem("Home",tabName =  "home", icon = icon("home")),
        # menuItem("Comparison", tabName = "coexnet", icon = icon("project-diagram")), 

        menuItem('Conservation',tabName = 'cons',icon = icon('sitemap')),
        
        menuItem("Data download",icon = icon("download"),href='http://labshare.cshl.edu/shares/gillislab/resource/CoCoBLAST')
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
        tabItem( # home page
            tabName = 'home',
            fluidPage(
                titlePanel(
                    h1(div('CoCoBLAST'), align = "center")

                ),
                fluidRow(
                    box(
 
                        status='primary',width = NULL,
                        solidHeader=F,
                        column(width=8,

                            p( h2("Welcome to CoCoBLAST!" ),
                                'How can we find gene regulatory differences between long-diverged species? The answer is coexpression.'
                            ), br(),
                            HTML('With <b>CoCoBLAST</b>,  you can input a gene of interest within a reference species and search for genes with conserved coexpression partners in a target species.'),
                    
                    

                            br(),
                            fluidRow(
                                
                                column(width=3) ,
                                column(width=6,
                                    p('To explore conserved genes across species, go here!' ),
                                    actionBttn(inputId='toConservation',label='Co-expression Conservation',icon=icon('arrow-right')),
                                    # actionButton(inputId='toConservation',label='Co-expression Conservation',icon=icon('arrow-right'))
                                ) ,
                                column(width= 3)
                            
                            ),
                            br(),br(),
                        
   
                            p(
                            
                                fluidRow(

                                    #  <a href='https://academic.oup.com/nar/article/48/W1/W566/5835815' target='_blank'> CoCoCoNet</a>
                                    column( width=7,
                                        tags$h4('The wide breadth of species surveyed'),
                                            HTML("We use aggregate coexpression networks from CoCoCoNet
                                            <sup><a href='https://academic.oup.com/nar/article/48/W1/W566/5835815' target='_blank'> 1</a></sup> 
                                            for 37 species chosen for their high quality RNA-seq expression data and orthology predictions from OrthoDB
                                            <sup><a href='https://academic.oup.com/nar/article/47/D1/D807/5160989' target='_blank'> 2</a> </sup>. 
                                            This large collection of data permits the evaluation of coexpression conservation across the tree of life: from fungi to animals to plants." 
                                            ),br(),

                                        tags$h4('Our co-expression conservation metric'),
                                            HTML("To quantify gene similarity between each pair of species, 
                                            we compare each gene’s top coexpression partners. We treat this as a supervised 
                                            learning task, using the ranks of the coexpression strengths from one species 
                                            to predict the top coexpression partners from the second species, 
                                            then repeating this task in the opposite direction, 
                                            finally averaging the scores. We refer to this as a measure of 
                                            <b>coexpression conservation</b> and note that it is formally equivalent to the 
                                            average area under the receiver operator characteristic curve (AUROC). ")

                                        

                                    

                                    ),
                                    column( width=5,
                                        div(img(src="method.png",height='90%',width='90%'), align='center'   )                       
                                    )
                                )

                            ),
                            

                            br(),
                            p(
                            tags$h4('Did you find CoCoBLAST Helpful?'),
                                HTML('If so, please consider citing us at: <br/><br/>
                                Megan Crow, Hamsini Suresh, John Lee, Jesse Gillis,<br/> <b>Coexpression reveals conserved mechanisms of transcriptional cell identity</b>, 
                                <i>bioRxiv</i>, 10 Nov, 2020, <a href="https://doi.org/10.1101/2020.11.10.375758" target="_blank">  https://doi.org/10.1101/2020.11.10.375758</a>'
                                )

                            ),

                            br(),
                            p(
                            tags$h4('Have a question or issue?'),
                                HTML('Feel free to contact John Lee at <b>johlee@cshl.edu</b> to report an issue, suggest an update, or ask for help. ')
                            

                            ),

                        ),
                        column(width=4, 
                            div(img(src="cococonet_v3tree-01.png", height= "70%",width= "80%"), align='center'   )                       
                        )
                    

                    )
                )
                    
                
              
            )
        ),

        tabItem(

            tabName = "cons",
            fluidPage(   
                titlePanel('CoCoBLAST'),
                HTML("To query co-expression conservation scores, input a gene from a reference species, select a target species, then hit “Submit”. CoCoBLAST will then do several tasks:"),
                tags$ol(
                    tags$li('It will find the 50 top co-expressed genes in the target species, showing each gene pair’s rank, as well as the reciprocal rank.'),
                    tags$li('Next, it looks for paralogs for each input gene in the reference species, maps these via orthology to the target species, and plots a heatmap and network of coexpression conservation scores.'),
                    tags$li(HTML("CoCoBLAST also showcases the expression of query genes when mapped to either Tabula Muris
                        <sup><a href='https://www.nature.com/articles/s41586-018-0590-4' target='_blank'>3</a> </sup> 
                        or single-cell RNA-seq data from Arabidopsis roots
                        <sup><a href='http://www.plantcell.org/content/31/5/993?ijkey=833f8c899258ff67f83c0dc5917375f1693aa70a&keytype2=tf_ipsecsha' target='_blank'>4,</a> </sup> 
                        <sup><a href='https://www.sciencedirect.com/science/article/pii/S1534580719301455?via%3Dihub' target='_blank'> 5,</a> </sup> 
                        <sup><a href='http://www.plantphysiol.org/content/179/4/1444?ijkey=41a5fbcf207863cc617cc115e9d9b1afee25fd4f&keytype2=tf_ipsecsha' target='_blank'> 6,</a> </sup> 
                        <sup><a href='https://www.sciencedirect.com/science/article/pii/S2211124719305273' target='_blank'> 7</a> </sup>.  
                        ") 
                    )
                ),

                
                # "),
                # maps these genes to the target species and a heatmap of scores, and a network between orthologs. CoCoCoBLAST also showcases the expression of these genes when mapped to either Tabula Muris 
                # <sup><a href='https://www.nature.com/articles/s41586-018-0590-4' target='_blank'>3</a> </sup> 
                # or Arabidopsis Root data 


                fluidRow(
                    box(
                        title='Co-expression Conservation',status='primary',width = NULL,
                        solidHeader=T,

                        #input box (typing)
                        # begin paste gene list
                        fluidRow(
                            # column( width=3),
                            
                            column( width=3.,
                                textInput(  "cons_input", 
                                    label ="Look for the following genes...",value = "Tubb4b"
                                    
                                ),                          
                                bsTooltip(                      
                                    id="cons_input", options = list(container = "body"),placement='right',
                                    title= "Mutiple genes can be queried simultaneously with mutiple genes separated by a comma."
                                ),                                  
                                selectInput(                    # select species
                                    inputId="cons_specA",selected='mouse',selectize=F,
                                    label='In this species...',
                                    choices=speciesList
                                ),
                                                            selectInput(                    # select species
                                    inputId="cons_specB",selected='human',selectize=F,
                                    label='Then compare them to orthologs in... ',
                                    choices=speciesList
                                ),

                                actionBttn(
                                    'cons_submit','Submit',icon=icon('search') 
                                )    ,

                                bsTooltip(id='cons_submit',placement = 'right',options = list(container = "body"),
                                    title = "Warning: Only click submit once: A delay may occur if our servers are busy."
                                )

                                # bsTooltip(                      
                                #     id='cons_submit', options = list(container = "body"),placement='right',
                                #     title= "Mutiple genes can be queried simultaneously. Mutple genes should be separated by a comma."
                                # ),                                  
                                
                            )

                        )
                    ),
                    box(title = 'Tables', status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE,collapsed=FALSE, 
                        fluidRow( 
                            column(width = 1,
                                htmlOutput('viewText11') # view
                            ),
                            column(width = 2,
                                uiOutput('tableSelector')
                            ),
                            column(width = 1,
                                htmlOutput('viewText22')     # for...
                            ),
                            column(width = 2,
                                uiOutput('geneSelector1')      
                            ),
                            column(width=1,
                                # htmlOutput('viewText3')     #  paralogs and orthologs
                            ),
                            column(width=2, 
                                # uiOutput('tnseSpecSelector')
                            ), 
                            column(width=3,
                                # uiOutput('bipart_slider')
                            )
                        ),br(),
                        fluidRow(
                            column(width = 8,
                                DT::dataTableOutput('cons_table')
                            ),
                            column(width=4,
                                uiOutput('downloadflag'),
                                uiOutput('option'),
                                uiOutput('connectionNumber'),
                                uiOutput('downloadButton')  
                            )
                        )
                    ),
                    box(title = 'Figures', status = 'primary', width=NULL, solidHeader=T, collapsible = TRUE,collapsed=FALSE,
                        fluidRow( 
                            column(width = 1,
                                htmlOutput('viewText1') # view
                            ),
                            column(width = 2,
                                uiOutput('figSelector')
                            ),
                            column(width = 1,
                                htmlOutput('viewText2')     # for...
                            ),
                            column(width = 2,
                                uiOutput('geneSelector')      
                            ),
                            column(width=2,
                                htmlOutput('viewText3')     #  paralogs and orthologs
                            ),
                            column(width=2, 
                                uiOutput('tnseSpecSelector')
                            ), 
                            column(width=2,
                                uiOutput('bipart_slider')
                            )
                            
                        ),
                        fluidRow(
                            column(width = 8,
                                plotOutput('mainFig',height='600px') # ,click = 'click_pt'
                            ),

                            column(width=4,
                                uiOutput('scoreThres'),
                                visNetworkOutput('bipartite',height='600px')
                            )
                        )

                    
                        
                    )

                )
            ) # end fluid page
        )  # end tab item (cons)

    )   # end tabitems
)   # end dashboard body           


# build everything
dashboardPagePlus(header,sidebar,body,enable_preloader = TRUE, sidebar_fullCollapse = T,loading_duration = 0)



