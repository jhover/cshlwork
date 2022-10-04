library(shiny,quietly=T)
library(shinyWidgets,quietly=T)
library(shinydashboard,quietly=T)
library(shinydashboardPlus,quietly=T)
library(networkD3)
library(DT)

header = dashboardHeader(title = "",
												 titleWidth=0)

sidebar = dashboardSidebar(
	width = 210,
	collapsed = FALSE,
	sidebarMenu(
		# pro tip - re order tabs to change landing page for testing.
		id = 'tabs',
		menuItem("Home", tabName = "home", icon=icon("house-user") ),
		menuItem("Comparative analyses", tabName = "primate", icon=icon("search") ),		
		menuItem("Data Download", icon=icon("download"), href='http://labshare.cshl.edu/shares/gillislab/resource/Primate_MTG_coexp') 
	) # end sidebarMenu 
)   # end sidebar

body = dashboardBody(
	tabItems(
		tabItem( tabName = 'home',
				 fluidRow(
				 	box(title = "Gene functional conservation across cell types and species", 
				 		status = 'primary', 
				 		width=NULL, 
				 		solidHeader=T, 
				 		collapsible = FALSE, 
				 		collapsed=FALSE,
            p(
              HTML('We aligned single-nucleus atlases of the middle temporal gyrus (MTG) of 5 primates (human, chimp, gorilla, macaque and marmoset) and identified 57 consensus cell types common to all species. We provide this resource for users to: <ol>
                <li>explore the conservation of gene expression across primates at single cell resolution,</li>
                <li>compare with the conservation of gene coexpression across metazoa, and</li>
                <li>identify genes with changes in expression or connectivity that drive the rapid evolution of the human brain</li>
               </ol> Given a query human gene (Gene A), this tool visualizes:')
              ),
             
            div(img(src="web-server-fig-01.png",height='90%',width='90%'), align='center'),
            br(),
            p(
              HTML('If you find this resource helpful, please consider citing us at: <br/>
              Hamsini Suresh, Megan Crow, Nikolas Jorstad, Rebecca Hodge, Ed Lein, Alexander Dobin, Trygve Bakken, Jesse Gillis,<br/> <b>Conserved coexpression at single cell resolution across primate brains</b>, 
                                <i>bioRxiv</i>, Sep, 2022' 
                                #<a href="https://doi.org/10.1101/2020.11.10.375758" target="_blank">  https://doi.org/10.1101/2020.11.10.375758</a>'
              )

            ),
            HTML("If you have any questions or suggestions to improve the tool, please feel free to contact Hamsini Suresh at <b>suresh@cshl.edu</b>")
           )  # end box
	       ) # end fluidRow                 				 	
		),  # end tabItem 'home'						
   
   tabItem( tabName = 'primate',
				fluidRow(
					box(title = "Cross-species coexpression conservation Query", 
						status = 'primary', 
						width=NULL, solidHeader=T, collapsible = FALSE, collapsed=FALSE,
							selectInput(
								inputId = 'gene',
								label = 'Human Gene Name',
                selected = 'GAD2',
								choices = genes,
								selectize = T
							)
					), # end box
             
          box(title = "Gene expression profile in human and NHP", status = 'primary', width=NULL, solidHeader=T, collapsible = FALSE, collapsed=FALSE,                                                     
              fluidRow(
                        column(
                            width = 9,
                            HTML("Expression profile of query gene across 57 consensus cell types in human (red) and non-human primates (grey):"), br(),
                            HTML("Note: Dashed lines separate consensus cell types by cell class (excitatory neurons in the left, inhibitory neurons in the middle, and non-neurons in the right)"),
                            plotOutput(outputId ="human_nhp_exp")                                  
                        ),                           
                        column(
                            width = 3,
                            HTML("Expressolog score: rank-standardized expression profile similarity of 1:1 orthologs relative to all other genes. Expressolog scores over cell types in each cell class:"),
                            DT::dataTableOutput("class_explog_table")                          
                        )
              )                                      
					),
          
          box(title = "Expressolog scores across all primates", status = 'primary', width=NULL, solidHeader=T, collapsible = FALSE, collapsed=FALSE,
              HTML("Expressolog scores calculated over 57 cell types for each primate pair:"), 
              fluidRow(
                        br(),
                        column(
                            width = 5,
                            plotOutput(outputId ="expressolog", height = 300)                              
                        ),                           
                        column(
                            width = 7,
                            DT::dataTableOutput("exp_ortho_table")                          
                        )
              )                            
					),
            
          box(title = "Expressolog scores across cell lineages", status = 'primary', width=NULL, solidHeader=T, collapsible = FALSE, collapsed=FALSE,
              HTML("The consensus cell types are split into 3 classes, 6 subtypes (2 from each major class), and 10 meta-clusters (3 from each inhibitory subtype and 2 from each excitatory subtype), and the expressolog scores calculated over cell types within each group are shown below:"), 
              diagonalNetworkOutput(outputId ="lineage_aurocs", height = 300)              
              #div(img(src="Rplots-zissou1-colorbar.png",height='50%',width='50%'), align='center')  
         ),
            
          box(title = "Expressolog score relative to other genes", status = 'primary', width=NULL, solidHeader=T, collapsible = FALSE, collapsed=FALSE,
          p(
              HTML("Distributions of expressolog scores for 14,131 genes are shown in boxplots, with the expressolog scores for the query gene marked in red. If the query gene lies below the threshold (grey line at AUROC = 0.55) in the Human-NHP boxplots, it is assumed to have diverged expression profiles between human and non-human primates (NHP) in those classes."),
              br()
              ),             
              plotOutput(outputId ="human_nhp_boxplots", height = 250)                             
					),
            
					box(title = "Coexpression conservation across metazoa", status = 'primary', width=NULL, solidHeader=T, collapsible = FALSE, collapsed=FALSE,
              HTML("p-values quantifying human-specific differential coexpression relative to other vertebrates (12 other mammals, chicken and zebrafish:"),
              DT::dataTableOutput("pval_table"), br(),	               						
              fluidRow(
                        column(
                            width = 4,
                            HTML("Distributions of coexpression conservation between human-mammal, mammal-mammal, and vertebrate-mammal"),
                            plotOutput(outputId ="vertBoxplots", height = 350)                            
                        ),                           
                        column(
                            width = 8,                            
                            p(style="text-align: center;",
                              HTML("Gene coexpression conservation across 19 animals")
                              ),                 
                            plotOutput(outputId ="aucHeat", height = 350)                         
                        )
              )                                                        							             
					) #end box.     
				) # end fluidRow
		) # end tabItem primate
	) # end tabItems
) # end DashboardBody



             
             

dashboardPage(
	header, sidebar, body
)