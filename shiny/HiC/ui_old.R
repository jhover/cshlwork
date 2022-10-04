useShinyjs()


header = dashboardHeaderPlus(
			         title = 'Hi-C Aggregate Resource',
				     titleWidth=210
				 )



sidebar = dashboardSidebar( 
			       width = 210,
			           sidebarMenu(
					               # pro tip - re order tabs to change landing page for testing.
					               id = 'tabs', 
						               menuItem("Home", tabName =  "home", icon = icon("home")),
						               # menuItem('Gene Comparison', tabName = 'gcomp', icon = icon('sitemap')),
						               # menuItem('Gene Data', tabName = 'gdata', icon = icon('sitemap')),
						               menuItem("Download", tabName = "download", icon = icon("download"))
							               # menuItem("Download", icon = icon("download"), href='http://labshare.cshl.edu/shares/gillislab/resource/HSC_atlas')
							           )
			       )

# NOTES for Nathan: 
### commas are super annoying in ui.R. Fair warning.
### use combinations of fluidRow() and column() to generate grids. 
### use box() to make groups prettier 
### column widths should be integers that sum to 12 
### large select/picker inputs should be rendered empty and then updated in server for speed
### 

 
body = dashboardBody(
		         tabItems( 
				          tabItem(
						              tabName = 'home',
							                  fluidPage(
										                    titlePanel(
													                           h1(
																                             div('Hi-C Aggregate Resource'),
																			                            align = "center")
																                   )
												                ),
						              box(
								                  status='primary', width = NULL, solidHeader=TRUE,
										                  column(
													                     width=12,
															                         # p( h2("A Meta-Analytic Single-Cell Atlas of Mouse Bone Marrow Hematopoietic Development" ),
															                         #     'Harris, Lee, and Gillis 2021'
															                         # ), 
															                         br(),
															                         HTML(paste0('We provide a Hi-C data aggregate resource. The aggregate ',
																			                                     'is constructed of at least 3500 human, 6600 mouse and 400',                                ' drosophila SRA runs.')),
													             
													                     br(),
															                         HTML(paste0('The data is available at various resolutions under the ',
																			                                     '"Download" tab ')),
													                     fluidRow(
																                              column(
																				                                 width=3,
																								                             h4('Full Aggregate for All Resolutions' ),
																								                             p(paste0('Download aggreagted binsXbins contact frequency data'))
																											                             ),
																                              column(
																				                                 width=3,
																								                             h4('Gene Vectors for All Resolutions' ),
																								                             p(paste0('Download contact frequency of a gene with every other bin in bed file format'))
																											                             )
																			                          ),
													                     # fluidRow(
													                     #     column(width=3,
													                     #         actionBttn(inputId='to_gdata',label='to tab1',icon=icon('arrow-right') )
													                     #     ),
													                     #     column(width=3,
													                     #         actionBttn(inputId='to_gcomp',label='to tab2',icon=icon('arrow-right'))
													                     #     )
													                     # ),

													                     br(),
															                         p(
																		                           tags$h4('Did you find this app helpful?'),
																					                           HTML('If so, please consider citing us at: <br/><br/>
																									                        Ruchi Lohia, Nathan Fox, and Jesse Gillis<br/> <b>title of paper</b>, 
																												                        <i>BioRxiv</i>, pub date,<br/> <a href="https://doi.org/10.1101/2021.08.12.456098" target="_blank">link text</a>'
																												                        )
																					                       ),

													                     br(),
															                         p(
																		                           tags$h4('Have a question or issue?'),
																					                           HTML('Feel free to contact John Lee at <b>johlee@cshl.edu</b> to report an issue, suggest an update, or ask for help. ')
																					                       )
															                     )
										              )   
							              ), # end home tab

				          # tabItem(
				          #     tabName = "gdata",
				          #     # Helper text goes here
				          #     fluidPage(   
				          #         titlePanel('Gene Data')
				          #     ),

				          #     box(status='primary', width = NULL, solidHeader=TRUE,
				          #         fluidRow(
				          #             column(
				          #                 width=5,
				          #                 # h2('Title'),
				          #                 # h3('Subtitle'),
				          #                 # p('some text'),

				          #                 # br(),
				          #                 div("NOT FUNCTIONAL", style="font-weight: bold; color: red; font-size: 2em;"),
				          #                 pickerInput(                    # select species
				          #                     inputId="gdata_species",
				          #                     label='Species', multiple=FALSE,
				          #                     choices=SPECIESLIST ,selected=SPECIESLIST[1],
				          #                     width='fit',
				          #                     inline=FALSE
				          #                 ),
				          #                 selectizeInput(                    # select gene for the x axis
				          #                     inputId="gdata_gene",
				          #                     label='Gene 1', multiple=FALSE,
				          #                     choices=NULL, selected=NULL,    # updated in server.R
				          #                     width='16em'
				          #                 ),
				          #                 fluidRow(
				          #                     column(
				          #                         width = 6,
				          #                         textInput(                    
				          #                             inputId="gdata_region1", label='genomic region 1', placeholder='chr1:1-100'
				          #                         ),
				          #                         pickerInput(                    # select gene for the x axis
				          #                             inputId="gdata_binsize",
				          #                             label='Bin Size', multiple=FALSE,
				          #                             choices=c("1kbp", "5kbp", "10kbp", "40kbp",
				          #                                       "100kbp", "250kbp", "500kbp"),
				          #                             selected="10kbp",    # updated in server.R
				          #                             width='fit',
				          #                             inline=FALSE,
				          #                             options = pickerOptions(
				          #                                 virtualScroll=20,       # max 20 items appear
				          #                                 liveSearch=FALSE         # type to search for gene
				          #                             )
				          #                         ),
				          #                     ),
				          #                     column(
				          #                         width = 6,
				          #                         textInput(                    
				          #                             inputId="gdata_region2", label = 'genomic region 2',
				          #                             placeholder= 'chr2:100-200'
				          #                         )                    
				          #                     ),
				          #                 ),
				          #                 fluidRow(
				          #                     box(
				          #                         # specify download options
				          #                         prettyRadioButtons(
				          #                             inputId = "downloadOption",
				          #                             label = "Download:",
				          #                             thick = TRUE,
				          #                             choiceNames = c('Gene x Bins Data', 'Heatmap Data', 'Heatmap Image'),
				          #                             choiceValues = c(1, 2, 3),
				          #                             animation = "pulse",
				          #                             status = "info", outline = T
				          #                         ),

				          #                         # download executer
				          #                         downloadBttn('downloader',label='Download:',color='primary',style='stretch',size='sm',block=T)   
				          #                     )
				          #                 )
				          #             ), 

				          #             # figures
				          #             column(width = 7, 
				          #                 plotOutput(
				          #                     'tab1_heatmap', 
				          #                     height='600px',
				          #                     dblclick = "heatmap_dblclick",
				          #                     brush = brushOpts(
				          #                         id = "heatmap_brush",
				          #                         resetOnNew = TRUE
				          #                     )         
				          #                 )
				          #             ), 
				          #             # column( width = 6,
				          #             #     fluidRow(
				          #             #         DT::dataTableOutput('tab1_marker_table' , height='300px' )
				          #             #     )
				          #             # )
				          #         )
				          #         
				          #     )

				          # ) , # end tab item (gdata) 

				          # tabItem(
				          #     tabName = "gcomp",
				          #     # Helper text goes here
				          #     fluidPage(   
				          #         titlePanel('Gene Comparison')
				          #     ),

				          #     # main tab with figures
				          #     box(status='primary', width = NULL, solidHeader=F,
				          #         fluidRow(
				          #             column(
				          #                 width=5,
				          #                 # h2('Title'),
				          #                 # h3('subtitle'),
				          #                 p('Explore characteristics of genes similar to a target gene.' ),
				          #                 h4("Instructions", style="margin-top: 1.25em;"),
				          #                 tags$ol(
				          #                     tags$li("Select a species and a selection metric."),
				          #                     tags$li("Select a target gene."),
				          #                     tags$li(paste("Choose how many genes will be examined. These will",
				          #                                   "be the genes that are the most similar to the target",
				          #                                   "gene, based on the selection metric")),
				          #                     tags$li("Select an x-axis metric, and a y-axis metric.")
				          #                 ),
				          #                 p(paste("The resulting scatterplot will show each of the N examined genes,",
				          #                         "plotted by their similarity to the target gene on two different metrics")),
				          #                 p(paste("For example, we select 50 genes from Pearson Correlation to target Gene1",
				          #                         "and plot by Contact Frequency and Shared Contacts. This scatterplot will",
				          #                         "show the 50 genes whose Hi-C profiles are most highly correlated with",
				          #                         "Gene1's Hi-C profile. Those 50 genes will be plotted on the number of",
				          #                         "contacts they have with Gene1 against the number of shared contacted",
				          #                         "regions with Gene1.")),
				          #                 # input selection
				          #                 br(),
				          #                 h4("Pearson Correlation is not functional.",
				          #                    style="color: red;"),
				          #                 br(),
				          #                 pickerInput(                    # select species
				          #                     inputId="gcomp_species",
				          #                     label='Species', multiple=FALSE,
				          #                     choices=SPECIESLIST ,selected=SPECIESLIST[1],
				          #                     width='fit',
				          #                     inline=FALSE
				          #                 ),
				          #                 pickerInput(                    # select species
				          #                     inputId="gcomp_resolution",
				          #                     label='Resolution', multiple=FALSE,
				          #                     # choices=c("10kbp", "40kbp"), selected="10kbp",
				          #                     choices=NULL ,selected=NULL,
				          #                     width='fit',
				          #                     inline=FALSE
				          #                 ),
				          #                 pickerInput(                    # select species
				          #                     inputId="gcomp_metric1",
				          #                     label='Compared Genes Selection Metric', multiple=FALSE,
				          #                     choices=names(METRICS), selected='Pearson Correlation',    # updated in server.R
				          #                     width='fit',
				          #                     inline=FALSE
				          #                 ),
				          #                 # pickerInput(                    # select species
				          #                 selectizeInput(                    # select species
				          #                     inputId="gcomp_gene",
				          #                     label='Target Gene', multiple=FALSE,
				          #                     # choices=c(paste("human", "gene", 1:9, sep = "_")), 
				          #                     choices=NULL, selected=NULL,    # updated in server.R
				          #                     width='16em'  # ,
				          #                     # inline=FALSE,
				          #                     # options = pickerOptions(
				          #                     #     virtualScroll=10,       # max 20 items appear
				          #                     #     liveSearch=TRUE         # type to search for gene
				          #                     # )
				          #                 ),
				          #                 numericInput(
				          #                     inputId="gcomp_topngenes",
				          #                     label="Top N Compared Genes",
				          #                     value=50,
				          #                     min=1,
				          #                     max=500,
				          #                     step=1,
				          #                     # width=100,
				          #                 ),
				          #                 pickerInput(                    # select species
				          #                     inputId="gcomp_metric2",
				          #                     label='Metric 1', multiple=FALSE,
				          #                     choices=names(METRICS), selected='Pearson Correlation',    # updated in server.R
				          #                     width='fit',
				          #                     inline=FALSE
				          #                 ),
				          #                 pickerInput(                    # select species
				          #                     inputId="gcomp_metric3",
				          #                     label='Metric 2', multiple=FALSE,
				          #                     choices=names(METRICS), # selected='Pearson Correlation',    # updated in server.R
				          #                     width='fit',
				          #                     inline=FALSE
				          #                 ),
				          #                 actionBttn(inputId="gcomp_submit", label="Submit")
				          #             ),

				          #             # scatter plot
				          #             column(
				          #                 width = 7, 
				          #                 br(),
				          #                 br(),
				          #                 plotOutput('tab2_scatter', height='50em')  #, width="600px")
				          #             )                        
				          #         ),
				          #         br(),
				          #         fluidRow(
				          #             column(
				          #                 width = 2,
				          #                 downloadBttn(outputId = "gcomp_img_download", label = "Download Figure")
				          #             ),
				          #             column(
				          #                 width = 10,
				          #                 downloadBttn(outputId = "gcomp_data_download", label = "Download Data")
				          #             )
				          #         )
				          #     )        
				          # ),   # end tab item (gcomp)
				          tabItem(
						              tabName = 'download',
							                  fluidPage(
										                    titlePanel(
													                           h1(
																                             div('Download Data')
																			                         )
																                   )
												                ),
						              box(
								                  status='primary', width = NULL, solidHeader=TRUE,
										                  column(
													                     width=12,
															                         h4(
																		                           div('Full Aggregate for All Resolutions')
																					                       ),
													                     br(),
															                         div(paste("Download a genome-wide Hi-C aggregate dataset for a",
																			                                 "species/resolution combination. The downloaded file is", "in "),
																		                             a("Hi-C Explorer", href="https://hicexplorer.readthedocs.io/en/latest/"),
																					                             paste(" format. Files range from a few MB to tens of GB, so download
																									                                 times may be long.")),
																													                     br(),
																													                     pickerInput(                    # select species
																																	                         inputId="download_aggspecies",
																																				                         label='Species', multiple=FALSE,
																																				                         choices=SPECIESLIST ,selected=SPECIESLIST[1],
																																							                         width='fit',
																																							                         inline=FALSE
																																										                     ),
																									                       pickerInput(                    # select species
																													                           inputId="download_aggresolution",
																																                           label='Resolution', multiple=FALSE,
																																                           choices=NULL ,selected=NULL,
																																			                           width='fit',
																																			                           inline=FALSE
																																						                       ),
																									                       hr(),
																											                           downloadBttn(outputId = "download_aggdownload", label = "Download")
																											                       )
																					                 ),
													             box(
															                 status='primary', width = NULL, solidHeader=TRUE,
																	                 column(
																				                    width=12,
																						                        h4(
																									                          div('Gene Vectors for All Resolutions')
																												                      ),
																				                    br(),
																						                        pickerInput(                    # select species
																										                            inputId="download_genespecies",
																													                            label='Species', multiple=FALSE,
																													                            choices=SPECIESLIST ,selected=SPECIESLIST[1],
																																                            width='fit',
																																                            inline=FALSE
																																			                        ),
																				                    pickerInput(                    # select species
																								                        inputId="download_generesolution",
																											                        label='Resolution', multiple=FALSE,
																											                        choices=NULL ,selected=NULL,
																														                        width='fit',
																														                        inline=FALSE
																																	                    ),
																				                    selectizeInput(                    # select species
																								                           inputId="download_gene",
																											                           label='Target Gene', multiple=FALSE,
																											                           choices=NULL, selected=NULL,    # updated in server.R
																														                           width='16em'
																														                       ),
																				                    hr(),
																						                        downloadBttn(outputId = "download_genedownload", label = "Download")
																						                    )
																	             ),
													             box(
															                 status='primary', width = NULL, solidHeader=TRUE,
																	                 column(
																				                    width=12,
																						                        # img(src="http://labshare.cshl.edu/shares/gillislab/resource/HiC/link_clipart.jpg", alt="hyperlink symbol",
																						                            # width="100", height = "100"),
																						                        h4(
																									                          tags$a("Link to all aggregates",
																													                               href="http://labshare.cshl.edu/shares/gillislab/resource/HiC/",
																																                                     style="color: #222222;",
																																                                     onmouseover="this.style.color='#222222'; this.style.textDecoration='underline'",
																																				                                   onmouseout="this.style.color='#222222'; this.style.textDecoration='none'")
																												                      )
																						                    )
																	             )
														             ) # end download tab
										      )   # end tabitems
							      )   # end dashboard body


					  # actually builds the page
					  dashboardPagePlus(header,sidebar,body,enable_preloader = TRUE, sidebar_fullCollapse = T,loading_duration = 0)



