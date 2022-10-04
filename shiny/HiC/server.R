library(Matrix)
library(hdf5r)
library(ComplexHeatmap) # used for pretty and interactive heatmaps
library(plotly)


### TODO fix zoom in on heatmap
# plot_heatmap <- function(region_input,brush = NULL){
#     # TODO test on Rstudio
#     m = matrix(runif(100,0,100) ,10,10)
#     rownames(m) = letters[1:10]
#     colnames(m) = letters[1:10] 
# 
#     g = ComplexHeatmap::Heatmap(m) 
#     
#     # use ranges to filter the heatmap
#     if (!is.null(brush)){
#         lt = ComplexHeatmap:::get_pos_from_brush(brush)
#         pos1 = lt[[1]]
#         pos2 = lt[[2]]
#         pos = selectArea(g, mark = FALSE, pos1 = pos1, pos2 = pos2, 
#             verbose = FALSE, ht_pos = ComplexHeatmap::ht_pos_on_device(g))
#     
# 
#         row_index = unlist(pos[1, "row_index"])
#         column_index = unlist(pos[1, "column_index"])
#         m = g@ht_list[[1]]@matrix
#         ht_select = Heatmap(m[row_index, column_index, drop = FALSE],
#             col = g@ht_list[[1]]@matrix_color_mapping@col_fun,
#             show_heatmap_legend = FALSE,
#             cluster_rows = FALSE, cluster_columns = FALSE)
#         draw(ht_select)
# 
#     }
# 
#     
#     return(g)
# }
# 
# is_valid_region <- function(region) {
#     # do stuff to verify input (should be a named list)
#     valid = TRUE 
#     if (valid){
#         return(region)
#     } else{
#         return(NULL)
#     }
# 
# }

# plot_gene_comp <- function(input, genes) {
#     # if (input$gcomp_submit == 0) {
#     #     placeholder_df <- data.frame(x = 1:5, y = 4:8, gene = rep("gene", 5))
#     #     placeholder_p <- (
#     #         ggplot2::ggplot(data = placeholder_df)
#     #         + ggplot2::geom_point(ggplot2::aes(x = x, y = y))
#     #         + ggplot2::theme_classic()
#     #         # + ggplot2::coord_fixed()
#     #         + ggplot2::theme(axis.text = ggplot2::element_text(size=12),
#     #                          axis.title = ggplot2::element_text(size=22))
#     #         + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::unit(c(10, 0, 0, 0), "mm")),
#     #                          axis.title.y = ggplot2::element_text(margin = ggplot2::unit(c(0, 10, 0, 0), "mm")))
#     #     )
#     #     return(list(plot = placeholder_p, data = placeholder_df))
#     # }
#     selection_metric <- METRICS[input$gcomp_metric1]
#     metric_1 <- METRICS[input$gcomp_metric2]
#     metric_2 <- METRICS[input$gcomp_metric3]
#     hf <- H5File$new(file.path(DATADIR, input$gcomp_species,
#                                paste0(input$gcomp_species, "_gene_gene_metrics.hdf5")),
#                      mode = "r")
# 
#     if (any(c(selection_metric, metric_1, metric_2) %in% c("coexpression"))) {
#         genes <- hf[[file.path(input$gcomp_resolution, "coexpression", "regions", "extra")]][]
#     } else {
#         genes <- hf[[file.path(input$gcomp_resolution, "contacts", "regions", "extra")]][]
#     }
# 
#     target_genes <- hf[[file.path(input$gcomp_resolution, selection_metric, "regions", "extra")]][]
#     m1_genes <- hf[[file.path(input$gcomp_resolution, metric_1, "regions", "extra")]][]
#     m2_genes <- hf[[file.path(input$gcomp_resolution, metric_2, "regions", "extra")]][]
# 
#     target_gene_index <- grep(input$gcomp_gene, target_genes, fixed=TRUE)[1]
#     m1_gene_index <- grep(input$gcomp_gene, m1_genes, fixed=TRUE)[1]
#     m2_gene_index <- grep(input$gcomp_gene, m2_genes, fixed=TRUE)[1]
# 
#     target_vector <- hf[[file.path(input$gcomp_resolution, selection_metric, "matrix")]][target_gene_index, ]
#     names(target_vector) <- target_genes
#     target_vector <- target_vector[genes]
#     m1_vector <- hf[[file.path(input$gcomp_resolution, metric_1, "matrix")]][m1_gene_index, ]
#     names(m1_vector) <- m1_genes
#     m1_vector <- m1_vector[genes]
#     m2_vector <- hf[[file.path(input$gcomp_resolution, metric_2, "matrix")]][m2_gene_index, ]
#     names(m2_vector) <- m2_genes
#     m2_vector <- m2_vector[genes]
# 
# 
#     target_vector_genes <- names(target_vector)[order(target_vector, decreasing=TRUE)[1:(input$gcomp_topngenes+1)]]
#     if (input$gcomp_gene %in% target_vector_genes) {
#         target_vector_genes <- target_vector_genes[target_vector_genes != input$gcomp_gene]
#     } else {
#         target_vector_genes <- target_vector_genes[1:input$gcomp_topngenes]
#     }
#     x <- m1_vector[target_vector_genes]
#     y <- m2_vector[target_vector_genes]
#     hf$close_all()
#     df <- data.frame(x = x, y = y, gene = target_vector_genes)
#     colnames(df) <- c(metric_1, metric_2, "gene")
#     p <- (
#         ggplot2::ggplot(data = df)
#         + ggplot2::geom_point(ggplot2::aes_string(x = colnames(df)[1], y = colnames(df)[2]))
#         + ggplot2::theme_classic()
#         # + ggplot2::coord_fixed()
#         + ggplot2::theme(axis.text = ggplot2::element_text(size=12),
#                          axis.title = ggplot2::element_text(size=22))
#         + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::unit(c(10, 0, 0, 0), "mm")),
#                          axis.title.y = ggplot2::element_text(margin = ggplot2::unit(c(0, 10, 0, 0), "mm")))
#     )
#     return(list(plot = p, data = df))
# }


function(input,output,session){


    {   # home page stuff
        # Switch tabs
        observeEvent(input$to_gdata,{        
            updateTabItems(session,'tabs','gdata')
        })
    
        observeEvent(input$to_gcomp,{        
            updateTabItems(session,'tabs','gcomp')
        })

    }   # end home page stuff


    # {   # tab 1 stuff

    #     species1 <- reactive({ input$gdata_species })

    #     # gene_lists 'react' to the species input
    #     genelist1 <- reactive({
    #         hf <- H5File$new(file.path(DATADIR, input$gdata_species,
    #                                    paste0(input$gdata_species, "_gene_gene_metrics.hdf5")),
    #                          mode = "r")
    #         genes <- hf[[file.path("40kbp", "shared_contacts", "regions", "extra")]][]
    #         hf$close_all()
    #         return(genes)
    #     })

    #     # did the gene_list (ie.species) change?? 
    #     observeEvent(species1(), {
    #         updateSelectizeInput(session,
    #             'gdata_gene', choices = genelist1(), server = TRUE
    #         )
    #     })


    #     # # reacts to 'textinputs'
    #     region_input = reactive({
    #         reg = list(region1 = input$gdata_region1 , region2 = input$gdata_region2  )
    #           # not yet implemented - returns boolean if valid
    #         is_valid_region(reg) # returns the region as a named list if true, else return some fail flag            
    #     })


    #     # update the heat map if the inputs change
    #     observeEvent(region_input(), {
    #         if (! is.null(region_input()) ) {
    #             # use the region to get a matrix and plot the heatmap
    #             
    #             output$tab1_heatmap <- renderPlot({
    #                 plot_heatmap(region_input() , brush = input$heatmap_brush )
    #             })
    #         }

    #     })

    #     # did the user brush and double click on the heatmap?
    #     # NOTE: ranges gives the x and y coordiantes (not indices)
    #     ranges <- reactiveValues(x = NULL, y = NULL)
    #     observeEvent(input$heatmap_dblclick, {
    #         brush <- input$heatmap_brush
    #         if (!is.null(brush)) {  # does a brush exist? if so get values
    #             ranges$x <- c(brush$xmin, brush$xmax)
    #             ranges$y <- c(brush$ymin, brush$ymax)

    #         } else {    # no bush? reset axes
    #             ranges$x <- NULL
    #             ranges$y <- NULL
    #         }
    #     })

    #     # observeEvent(ranges(),
    #     #     output$tab1_heatmap <- renderPlot({
    #     #         plot_heatmap(region_input() , ranges() )   
    #     #     })
    #     # )

    #     # Download execution options 
    #     output$downloader <- downloadHandler(
    #         filename = function() {
    #             'hic_download.tsv'
    #         },
    #         content = function(con) {
    #             if(length(input$downloadOption)==0){}
    #             else{switch(input$downloadOption,
    #                 '1' = {
    #                     tmp = data.frame(col1 = letters[1:5],col2= LETTERS[1:5])
    #                     write.table(tmp, con , row.names = F,sep='\t') 
    #                 },
    #                 '2' = {
    #                     tmp = data.frame(col1 = LETTERS[1:5],col2= letters[1:5])
    #                     write.table(tmp, con , row.names = F,sep='\t') 
    #                 }
    #             )}
    #         } 
    #     )




    # }   # end  tab 1 stuff


    # {   # tab 2 stuff

    #     gcomp_res <- reactive({
    #         hf <- H5File$new(file.path(DATADIR, input$gcomp_species,
    #                                    paste0(input$gcomp_species, "_gene_gene_metrics.hdf5")),
    #                          mode = "r")
    #         resolutions <- names(hf)[endsWith(names(hf), "bp")]
    #         hf$close_all()
    #         return(resolutions)
    #     })

    #     observeEvent(input$gcomp_species, {
    #         updatePickerInput(session,
    #             'gcomp_resolution', choices = gcomp_res()
    #         )
    #     })

    #     gcomp_params <- bindCache(
    #         reactive({
    #             return(c(input$gcomp_species, input$gcomp_metric1 , input$gcomp_metric2, input$gcomp_metric3))
    #         }),
    #         (input$gcomp_metric1 == "Co-Expression"
    #          || input$gcomp_metric2 == "Co-Expression"
    #          || input$gcomp_metric3 == "Co-Expression")
    #     )

    #     # gcomp_params <- reactive({
    #     #     return(c(input$gcomp_species, input$gcomp_metric1,
    #     #              input$gcomp_metric2, input$gcomp_metric3))
    #     # })

    #     genelist2 <- reactive({
    #         if (any(c(input$gcomp_metric1, input$gcomp_metric2, input$gcomp_metric3)
    #                 %in% c("Co-Expression"))) {
    #             gtype <- "coexp"
    #         } else {
    #             gtype <- "general"
    #         }
    #         genes <- read.table(file.path(DATADIR, input$gcomp_species,
    #                                       paste(input$gcomp_species, gtype, "genelist.csv", sep = "_")),
    #                             header = FALSE, sep = ",")[, 1]
    #     })

    #     # did the gene_list (ie.species) change?? 
    #     observeEvent(gcomp_params(), {
    #         updateSelectizeInput(session,
    #             'gcomp_gene', choices = genelist2(), server = TRUE
    #         )
    #     })

    #     # { # placeholder stuff
    #     #     output$tab2_scatter <- renderPlot({
    #     #         plot_gene_comp(input, genelist2())[["plot"]]
    #     #     })

    #     #     # disable("gcomp_img_download")
    #     #     # disable("gcomp_data_download")
    #     # }

    #     observeEvent(input$gcomp_submit,
    #         {        
    #             # output$tab2_scatter <- renderPlot({
    #             #     plot_gene_comp(input, genelist2())[["plot"]]
    #             # })
    #             output$tab2_scatter <- renderCachedPlot(
    #                 {
    #                     plot_gene_comp(input, genelist2())[["plot"]]
    #                 },
    #                 {
    #                     input$gcomp_submit
    #                 }
    #             )
#   #               enable("gcomp_img_download")
#   #               enable("gcomp_data_download")
    #         }
    #     )

    #     output$gcomp_img_download <- downloadHandler(
    #         filename = function() {
    #             "gene_comparison.png"
    #         },
    #         content = function(con) {
    #             ggplot2::ggsave(con, plot = plot_gene_comp(input, genelist2())[["plot"]], dev = "png",
    #                             width = 8, height = 8, units = "in", dpi = 300)
    #         },
    #         contentType = "image/png"
    #     )

    #     output$gcomp_data_download <- downloadHandler(
    #         filename = function() {
    #             paste0("gene_comp_", input$gcomp_species, "_", input$gcomp_gene, 
    #                    "_by_", METRICS[input$gcomp_metric1], ".csv")
    #         },
    #         content = function(con) {
    #             write.table(plot_gene_comp(input, genelist2())[["data"]], file = con, sep = ",",
    #                         quote = FALSE, row.names = FALSE, col.names = TRUE)
    #         },
    #         contentType = "text/csv"
    #     )
    # }   # end tab 2 stuff
    
    {   # download page stuff

        downloadaggres <- reactive({
            resolutions <- c()
            for (entry in list.files(file.path(DATADIR, input$download_aggspecies),
                                     full.names = TRUE, pattern = "[0-9]+kbp_raw")) {
                if (file.info(entry)[1, "isdir"]
                     && file.exists(file.path(entry,
                                              paste(input$download_aggspecies, gsub("_raw", "",
                                                    basename(entry), fixed = TRUE), "aggregate.h5", sep = "_"))) ) {
                    resolutions <- c(resolutions, basename(entry))
                }
            }
            resolutions <- gsub("_raw", "", resolutions)
            resolutions <- names(sort(RESORDER[resolutions]))
            message('resolutions: ', resolutions)
            return(resolutions)
        })

        observeEvent(downloadaggres(), {
            updatePickerInput(session,
                'download_aggresolution', choices = downloadaggres()
            )
        })

        output$download_aggdownload <- downloadHandler(
            filename = function() {
                basename(file.path(DATADIR, input$download_aggspecies, paste0(input$download_aggresolution, "_raw"),
                                   paste(input$download_aggspecies, input$download_aggresolution, "aggregate.h5", sep = "_")))
            },
            content = function(con) {
                download_path <- file.path(DATADIR, input$download_aggspecies,
                                           paste0(input$download_aggresolution, "_raw"),
                                           paste(input$download_aggspecies, input$download_aggresolution,
                                                 "aggregate.h5", sep = "_"))
                file.copy(download_path, con)
            },
            contentType = "application/octet-stream"
        )


        downloadgeneres <- reactive({
            resolutions <- c()
            for (entry in list.files(file.path(DATADIR, input$download_genespecies),
                                     full.names = TRUE, pattern = "[0-9]+kbp_raw")) {
                if (file.info(entry)[1, "isdir"]
                     && file.exists(file.path(entry,
                                              paste(input$download_genespecies, gsub("_raw", "",
                                                    basename(entry), fixed = TRUE), "aggregate.h5", sep = "_")))) {
                    resolutions <- c(resolutions, basename(entry))
                }
            }
            resolutions <- gsub("_raw", "", resolutions)
            resolutions <- names(sort(RESORDER[resolutions]))
            return(resolutions)
        })

        observeEvent(input$download_genespecies, {
            updatePickerInput(session,
                'download_generesolution', choices = downloadgeneres()
            )
        })


        output$download_genedownload <- downloadHandler(
            filename = function() {
                paste0(
                    paste(input$download_gene, input$download_genespecies,
                          input$download_generesolution, sep = "_"),
                    ".bdg"
                )
            },
            content = function(con) {
                df <- read.table(file.path(DATADIR, input$download_genespecies,
                                           paste(input$download_generesolution, "raw", sep = "_"),
                                           paste(input$download_genespecies, input$download_generesolution,
                                                 "intervals.csv", sep = "_")),
                                 header = FALSE, sep = ",", col.names = c("chr", "start", "end"))
                hf_path <- file.path(DATADIR, input$download_genespecies,
                                     paste(input$download_generesolution, "raw", sep = "_"),
                                     paste(input$download_genespecies, input$download_generesolution,
                                           "gene_vectors.hdf5", sep = "_")
                                    )
                hf <- H5File$new(hf_path, mode = "r")
                gene_vector <- hf[[input$download_gene]][]
                df[["score"]] <- gene_vector
                header <- sprintf("track\ttype=bedGraph\tname=%s\n",
                                  paste(input$download_gene, input$download_genespecies,
                                  input$download_generesolution, "hic_contacts", sep = "_"))
                filename <- tempfile(pattern = "hic_download_gene_", tmpdir = file.path(DATADIR, "tmp"), fileext = ".bdg")
                cat(header, file = filename)
                write.table(df, sep = "\t", file = filename, append = TRUE,
                            quote = FALSE, row.names = FALSE, col.names = FALSE)
                file.copy(filename, con)
                file.remove(filename)
            },
            contentType = "text/plain"
        )

        download_genelist <- reactive({
            
            genes = read.csv(file.path(DATADIR, input$download_genespecies,
                               paste(input$download_genespecies, "gene_vector_genes.csv", sep = "_")),
                     col.names = c("gene"))[["gene"]]
            genes    
            
        })

        # did the gene_list (ie.species) change?? 
        observeEvent(input$download_genespecies, {
            tryCatch({
                updateSelectizeInput(session,
                    'download_gene', choices = download_genelist(),selected = download_genelist()[1],server = TRUE
                )
            },error=function(e){})
        })
    }   # end download page stuff
}
