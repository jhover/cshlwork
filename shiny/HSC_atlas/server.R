##################################################################################################*
############################################ Libraries ############################################
##################################################################################################*
# direc = "/srv/shiny-server/sample-apps/CoCoBLAST/" 

library(rhdf5) 
library(Matrix) 
# load("/ftp/CoCoBLAST/gene2go/voc.Rdata")           # descrptions of GO terms
# load(paste0(direc,'www/EGADlite.Rdata'))
 
# library(glue) 
# library(ComplexHeatmap)
# library(rhdf5,quietly=T)
# library(ggplot2)
# library(scattermore)
# library(visNetwork,quietly = T)



##################################################################################################*
############################################ Functions ############################################
##################################################################################################*
savefig <- function(fig){
    ggplot2::ggsave('~/test.png',plot=fig)
}

read_go_enrichment <- function(filepath = paste0(DATADIR,"cluster_GO_enrichment.csv" ) ) {
    go_enrich = read.delim(filepath, sep=",", stringsAsFactors=FALSE)
    colnames(go_enrich)[1] = "Cell.Type"
    return(go_enrich)
}

filter_markers <- function(markers, cell_type = 'all', max_rank = 100 , min_recurence = 2, min_auroc = 0.7) {
    # doesnt work!
    if (cell_type != 'all') {
        markers = markers[ markers$cell_type %in% cell_type]
    } 

    markers = markers %>% 
        filter(rank <= max_rank) %>% 
        filter(recurrence >= min_recurence) %>%
        filter(auroc >= min_auroc)
}

read_metamarkers <- function(filepath = paste0(DATADIR, "droplet_scNym_metamarkers.csv.gz") , clean = TRUE){
    df = as.data.frame(MetaMarkers::read_markers(filepath))
    if (clean){
        df = df[,2:13]
        df[,5:11] = round(df[,5:11] ,3)

        # df =filter_markers(df, cell_type='all', max_rank=100, min_recurence=0, min_auroc=0.0)
    }


    return(df)
}

# currently reads  in a csv file - neewd to convert to hdf5 for fast read
get_go_term_genes <- function(go_term , filepath = paste0(DATADIR,"go2gene.hdf5")){
    genes = rhdf5::h5read(filepath, go_term)
    return( genes)
}

plot_genes_on_umap <- function( gene_input =NULL,loom_file_path = paste0(DATADIR,"erythroid_and_monocyte_lineage_adata_no_gaps.loom") ){
    
    h5ad_path = sub(".loom",".h5ad", loom_file_path)
    umap = rhdf5::h5read(h5ad_path, "obsm/X_umap")
    

    # gene = 'Bin3'
    allgenes = rhdf5::h5read(loom_file_path , name = 'row_attrs/var_names')
    m = match(toupper(gene_input), toupper(allgenes))
    m = m[!is.na(m)]
    if (length(m) ==0 ) {
        gene_input=NULL
    }

    if (!is.null(gene_input)){
        scores = rhdf5::h5read(loom_file_path, 'matrix', index = list(1:ncol(umap) , m))
        umap = data.frame(UMAP_1 = umap[1,], UMAP_2 = umap[2,], score=  scores)

        
        fig = ggplot2::ggplot() + 
            ggrastr::geom_point_rast(data=umap, ggplot2::aes(x = UMAP_1, y = UMAP_2 , col=score), shape=".") + 
            #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
            ggplot2::scale_color_gradient(low = "gray80", 
            high = "darkred") + ggplot2::theme_classic() +
            ggplot2::labs(x = colnames(umap)[1], y =  colnames(umap)[2], col = "log(CPM + 1)") 
    } else{
        
        umap = data.frame(UMAP_1 = umap[1,], UMAP_2 = umap[2,])
        
        fig = ggplot2::ggplot(umap, ggplot2::aes(x = UMAP_1, y = UMAP_2,colour='gray')) +
            ggplot2::geom_point( colour="grey" , size=.005) + #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
            ggplot2::theme_classic()
            
    }

    return(fig)
    # ggplot2::ggsave('~/test.png')
    
}

plot_gene_violin <- function( gene_input , loom_file_path = paste0(DATADIR,"erythroid_and_monocyte_lineage_adata_no_gaps.loom") ,ct_attr = "scNym" ){

    # get expression for the gene of interest. 

 
    # all_genes = rhdf5::h5read(loom_file_path, name = "row_attrs/var_names")
    
    gene_ind = match(toupper(gene_input), toupper(GENELIST))
    gene_ind = gene_ind[!is.na(gene_ind)]

    cell_ids = rhdf5::h5read(loom_file_path, 'col_attrs/obs_names') # these are encoded as bytes
    ct = rhdf5::h5read(loom_file_path, paste0('col_attrs/',ct_attr)) # cell types

    scores = rhdf5::h5read(loom_file_path, 'matrix', index = list(1:length(cell_ids) , gene_ind)) # all cells, select genes   

    df = data.frame( score = scores, Cell_Type = ct )

    fig = ggplot2::ggplot( df, ggplot2::aes(x = Cell_Type, y = score ,fill=Cell_Type ) ) +
        ggplot2::geom_violin() + 
        ggplot2::theme_classic()+
        ggplot2::theme(legend.position = "none") +
        ggplot2::theme(axis.text.x = element_text(angle = 45,hjust = 1))

    return(fig)
    

}

# can be parsed beforehand for speed
plot_go_enrich_violin <- function(go_enrich, go_term, rdspath = ""){
    
    # go_enrich should be the dataframe of CT by GO term
    # currently, the entire violin needs to be computed at each iteration. 
    # Can save the violin instead since it does not vary with user input
    
    # do we have a filepath to the figure? 
    # filepath = paste0(DATADIR,"tab1_go_enrich.rds")
    if (file.exists(rdspath) ){
        fig = readRDS(rdspath)
        
        # get the GO term of interest
        go_oi = go_enrich[, c('Cell.Type',go_term) ,drop=FALSE] 
        colnames(go_oi) = c('Cell_Type',"AUROC")

        fig = fig +
            ggplot2::geom_point(data = go_oi , ggplot2::aes(x= Cell_Type, y = AUROC ) , shape = 8 ,size=5 ) 

        return(fig)
    }
        
    n=ncol(go_enrich)
    colnames(go_enrich)[1] = "Cell.Type"
    go_enrich_melted = tidyr::gather(go_enrich[,2:n])
    colnames(go_enrich_melted) = c("GO_Term","AUROC")

    
    ct = rep(go_enrich$Cell.Type, n-1 )

    go_enrich_melted['Cell_Type'] = ct
    
    # get the GO term of interest
    go_oi = go_enrich[, c('Cell.Type',go_term) ,drop=FALSE] 
    colnames(go_oi) = c('Cell_Type',"AUROC")

    fig = ggplot2::ggplot( go_enrich_melted , ggplot2::aes(x = Cell_Type, y = AUROC, fill = Cell_Type))  + 
        ggplot2::geom_violin() +
        geom_boxplot(width=0.25) +        
        ggplot2::theme_classic()+
        ggplot2::theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))

    fig = fig +
        ggplot2::geom_point(data = go_oi , ggplot2::aes(x= Cell_Type, y = AUROC ) , shape = 8 ,size=5 ) 



    return(fig)
}

plot_go_enrich_dot <- function(go_term ,loom_file_path = paste0(DATADIR,"erythroid_and_monocyte_lineage_adata_no_gaps.loom" ) , ct_attr = "scNym", tab="tab1"){
    # ct_attr = "FACS_labels"
    # loom_file_path = paste0(DATADIR,"mouse_hsc_labeled.loom" )
    # filepath = paste0(DATADIR, "dotplots/",ct_attr,'/',go_term, ".rds")
    filepath = paste0(DATADIR, "dotplots/",tab,"/",go_term, ".rds")
    if (file.exists(filepath)) {
        fig = readRDS(filepath)
        return(fig)
    }

    gene_set = get_go_term_genes(  go_term )

    # currently, loom files contain the dense matrix. can convert to 

    all_genes = rhdf5::h5read(loom_file_path, name = "row_attrs/var_names")
    gene_ind = match(toupper(gene_set), toupper(all_genes))
    gene_ind = gene_ind[!is.na(gene_ind)]

    cell_ids = rhdf5::h5read(loom_file_path, 'col_attrs/obs_names') # these are encoded as bytes
    ct = rhdf5::h5read(loom_file_path, paste0('col_attrs/',ct_attr )) # cell types
    study = rhdf5::h5read(loom_file_path, 'col_attrs/study_id') # cell types

    scores = rhdf5::h5read(loom_file_path, 'matrix', index = list(1:length(cell_ids) , gene_ind)) # all cells, select genes
    # scores = Matrix::Matrix(t(scores),sparse=TRUE)
    
    scores =t(scores)
    
    rownames(scores) =  as.character(all_genes[gene_ind])
    colnames(scores) = NULL

    dat = SummarizedExperiment::SummarizedExperiment( scores )
    fig = MetaNeighbor::plotDotPlot( dat, experiment_labels =as.character(study) ,celltype_labels = as.character(ct), gene_set=as.character(gene_set),normalize_library_size = FALSE,average_expressing_only = FALSE)
    return(fig)
}


plot_pt_genes_umap <- function(gene_input , loom_file_path = paste0(DATADIR,"erythroid_and_monocyte_lineage_adata_no_gaps.loom" ) ){
    umap1 = rhdf5::h5read( loom_file_path, 'col_attrs/monocle_UMAP1')
    umap2 = rhdf5::h5read( loom_file_path, 'col_attrs/monocle_UMAP2')

    study = as.character(rhdf5::h5read( loom_file_path, 'col_attrs/study_id'))

    df = data.frame(UMAP_1  = umap1, UMAP_2 = umap2, tmp =study)

    df = df %>% mutate(study = str_extract(tmp,"[0-9]+$"))
    df$study[is.na(df$study)] =""
    df$study = paste0(toupper(substr(df$tmp,0,1)) , df$study)


    m = match(toupper(gene_input), toupper(GENELIST))
    m = m[!is.na(m)]
    if (length(m) ==0 ) {
        gene_input=NULL
    }

    df$score = rhdf5::h5read(loom_file_path, 'matrix', index = list(1:nrow(df) , m))


    fig = ggplot2::ggplot() + 
        ggrastr::geom_point_rast(data=df, ggplot2::aes(x = UMAP_1, y = UMAP_2 , col=score), shape=".") + 
        ggplot2::facet_wrap(~study ,ncol=3 ,scales='free' ) + 
        ggplot2::scale_color_gradient(low = "gray80", high = "darkred") + 
        ggplot2::theme_void() +
        ggplot2::labs(x = colnames(df)[1], y =  colnames(df)[2], col = "log(CPM + 1)") 

    return(fig)
}

orthologMapper <- function(speciesA,genes, speciesB){
    ##  Maps genes in speciesA to orthologous genes in speciesB
    #   Input:  genes           - genes of SpeciesA
    #           speciesA        - name of reference species
    #           speciesB        - name of target species
    #   Output: temp            - orthologous genes in speciesB

    ## get NM orthologs  
    consDirec= "/var/ftp/data/CoCoBLAST/coexp_cons/"
    file1  = paste0( consDirec ,speciesA,"_",speciesB,"_CoCoBLAST_scores.hdf5")
    file2  = paste0( consDirec ,speciesB,"_",speciesA,"_CoCoBLAST_scores.hdf5")

    

    if (file.exists(file1) ) {
        orthoDat = rhdf5::h5read(file1, name = 'ortho_map' ) # read.table(file1,sep=",",stringsAsFactors=F,header=T)
    } else if(file.exists(file2)) {
        orthoDat = rhdf5::h5read(file2, name = 'ortho_map' ) #read.table(file2,sep=",",stringsAsFactors=F,header=T)
    } else {
        # no ortholog map found ....  for this pair...
        orthoDat = data.frame(matrix(NA, 1,2 ))
        colnames(orthoDat) =c( speciesA, speciesB)
    }
    
    tarorthogenes =  orthoDat[toupper( orthoDat[,speciesA]) %in% toupper(genes), ,drop=FALSE]  

    m = match(tarorthogenes[,speciesA] , unique(tarorthogenes[,speciesA]))
    n = match(tarorthogenes[,speciesB] , unique(tarorthogenes[,speciesB]))

    # as sparse matrix
    if (length(m) >0  & length(n) >0){
        orthoMat=Matrix::sparseMatrix(i = m, j=n )
        rownames(orthoMat) =unique(tarorthogenes[,speciesA]) #  getDisplayNames(speciesA, )
        colnames(orthoMat) =unique(tarorthogenes[,speciesB])  # getDisplayNames(speciesB, unique(tarorthogenes[,speciesB]))
    } else {
        orthoMat= NULL
    }



    return(list(tarorthogenes,orthoMat))
}


# keep --
geneConversionTable <- function(species){
    dataDirec= "/var/ftp/data/supplementData/"

    return( read.delim(paste0(dataDirec, "geneInfo/",species,"_info.csv"), stringsAsFactors=F, sep=","))
}


getNetworkID <- function(species, gene_list,flag='meta',verbose = FALSE){
    # given list of genes in any format , get the corresponding ensembl ID if it exists. otherwise, throw it away! Gene not included in STAR output
    # input genes is allowed to be anything now. EnsemblID , orthoID, proteinID, entrezID, genesymbol, etc  or any synonyms
    # outputs gene symbol if it exists, otherwise outputs ensemblID

    gene_list = gene_list[!is.na(gene_list)]
    gene_list = as.character(gene_list)

    # get gene info
    # geneInfo = read.delim(paste0("/var/ftp/CoCoBLAST/geneInfo/",species,"_geneInfo.tab"), stringsAsFactors=F, sep="\t")
    geneInfo = geneConversionTable(species)
    # second file are updated versions
  
    
    df = data.frame(input = gene_list,stringsAsFactors=FALSE)
    sets =c('NetworkIDs','EntrezID','GeneSymbol','EnsemblID','LocusTag')
    for (i in sets){
        
        k = match(toupper( gene_list) , toupper(geneInfo[,i]))
        df[i] =geneInfo[k, 'NetworkIDs']    # match everything to network IDs
    
    }

    df['found'] = df$NetworkIDs
    for (i in sets){
        ind = is.na(df$found)   # not found
        df$found[ind] = df[ind,i]
    }

    df['notFound']= is.na(df$found)

    # search in synonyms
    temp            = paste0("\\|", df$input[df$notFound] ,"\\|")
    not_found       = lapply(temp , grep , geneInfo$Synonyms,ignore.case=TRUE )
    not_found       = lapply(not_found, function(x) geneInfo[x,'NetworkIDs'])
    names(not_found)= df$input[df$notFound] 

    not_found = unlist(lapply(not_found, function(x) x[1] ))
    df$found[is.na(df$found)] = not_found
    df['notFound']= is.na(df$found)

    # search in xdb
    temp            = paste0(":", df$input[df$notFound] ,"\\|")
    not_found       = lapply(temp , grep , geneInfo$dbXrefs,ignore.case=TRUE )
    not_found       = lapply(not_found, function(x) geneInfo[x,'NetworkIDs'])
    names(not_found)= df$input[df$notFound] 

    not_found = unlist(lapply(not_found, function(x) x[1] ))
    df$found[is.na(df$found)] = not_found
    df['notFound']= is.na(df$found)



    # filter genes to only those in the HC/AA nets
    if (flag =='prio'){
        allGenes = h5read(paste0('/var/ftp/data/CoCoCoNet/networks/',species,'_',flag,'AggNet.hdf5'), name= 'row')
        m=match( allGenes, df$found)
        m=m[!is.na(m)]

        df = df[sort(m),]
    }


    N = length(unique(df$found))
    count = sum(df$notFound)
    if(verbose){
        if (count > 1 ) {
            
            if (count > 5)  {
                temp = paste(df$input[df$notFound][1:5],collapse=", ")
                temp = paste(temp ,"\n and",  count -5 , "others." )

                # showModal(paste("The following input genes were not found! \n ",temp), type='warning',closeOnClickOutside=T)
            } else {
                temp = paste(df$input[df$notFound],collapse=", ")
                # shinyalert(paste("The following input genes were not found! \n ",temp), type='warning',closeOnClickOutside=T)
            }
        }
    }


    # ensemblGenes should be used to match to the data, 
    # genes_not_found should be reported to the user
    # genes_found should be used as display names
    return(df)  
}

getDisplayNames <- function(species,ensGenes){
    # given the network IDs, convert them to symbols for display
    geneInfo = geneConversionTable(species)
    ind = match(toupper(ensGenes), toupper(geneInfo$NetworkIDs) )

    # get gene symbol
    displayIDs = geneInfo$GeneSymbol[ind ]
    

    # if gene symbol is not available, return what was given.
    ind2 = is.na(displayIDs)
    na_inds = ind[ind2]
    displayIDs[ind2] = ensGenes[ind2]


    return(displayIDs)
}


##################################################################################################*
############################################# Driver ##############################################
##################################################################################################*

# shiny_env=new.env()

function(input,output,session){
{   # home page

    # Switch tabs
    observeEvent(input$to_tab_1,{        
        updateTabItems(session,'tabs','tab_1')
    })
 
    observeEvent(input$to_tab_2,{        
        updateTabItems(session,'tabs','tab_2')
    })

    observeEvent(input$to_tab_3,{        
        updateTabItems(session,'tabs','tab_3')
    })
    
    observeEvent(input$to_tab_4,{        
        updateTabItems(session,'tabs','tab_4')
    })

}

# can wrap this in a function if desired.
{   # tab 1 - clusters
    updateSelectizeInput( session , 'tab_1_gene', choices = GENELIST,selected= 'Klf1' )
    golist =read.delim(paste0(DATADIR, 'go2description.tsv'), sep="\t" )[,1]

    updateSelectizeInput( session , 'tab_1_go', choices = golist, selected = golist[1] )

    # render plots
    output$tab1_umap =  renderPlot({
        tryCatch({
            fig = plot_genes_on_umap(gene_input =input$tab_1_gene, TAB1_LOOM_FILE_PATH  )
            fig
        }, error=function(e){})
    })

    output$tab1_gene_violin = renderPlot({
        tryCatch({
            fig = plot_gene_violin( gene_input =input$tab_1_gene , TAB1_LOOM_FILE_PATH,"scNym" )
            fig
        }, error=function(e){})
    })

    output$tab1_marker_table  = DT::renderDataTable({
        tryCatch({
            df = read_metamarkers(TAB1_MARKER_PATH , clean=TRUE)
            df[order(df$rank),]

        },error =function(e){}) 
        
        
        },rownames=FALSE,options = list(pageLength = 5, scrollX = TRUE)
    )

    output$tab1_go_enrich_violin  =  renderPlot({
        tryCatch({
            df = read_go_enrichment(TAB1_GO_ENRICH )
            
            goterm = strsplit(input$tab_1_go , ' - ')[[1]][1]
            goterm = sub('\\:','\\.' ,goterm)
            fig = plot_go_enrich_violin(df, goterm , paste0(DATADIR,"tab1_go_enrich.rds"))
            fig
        }, error= function(e){})
    })

    # takes a while to load Summarized Experiment package
    output$tab1_go_dot  =  plotly::renderPlotly({
        # if (!is.null(input$tab_1_go)) {

            tryCatch({
                
                goterm = strsplit(input$tab_1_go , ' - ')[[1]][1]
                goterm = sub('\\:','\\.' ,goterm)
                fig = plot_go_enrich_dot( goterm , TAB1_LOOM_FILE_PATH , 'scNym',tab="tab1")
                fig
            }, error= function(e){})
        # }
    })

} # end tab 1 section



{ # tab 2 - 
    updateSelectizeInput( session , 'tab_2_gene', choices = GENELIST,selected= 'Pf4' )
    golist =read.delim(paste0(DATADIR, 'go2description.tsv'), sep="\t" )[,1]
    updateSelectizeInput( session , 'tab_2_go', choices = golist, selected = golist[1] )

    # render plots
    output$tab2_umap = renderPlot({
        tryCatch({ #gene_input =input$tab_2_gene
            fig = plot_genes_on_umap(input$tab_2_gene, TAB2_LOOM_FILE_PATH  )
            fig
        }, error=function(e){})
    })

    output$tab2_gene_violin = renderPlot({
        tryCatch({
            fig = plot_gene_violin( gene_input =input$tab_2_gene , TAB2_LOOM_FILE_PATH ,"FACS_labels")
            fig
        }, error=function(e){})
    })

    output$tab2_marker_table  = DT::renderDataTable({
        tryCatch({
            df = read_metamarkers(TAB2_MARKER_PATH, clean=TRUE)
            df[order(df$rank),]

        },error =function(e){}) 
        
        
        },rownames=FALSE,options = list(pageLength = 5, scrollX = TRUE)
    )


    output$tab2_go_enrich_violin  = renderPlot({
        tryCatch({
            df = read_go_enrichment(TAB2_GO_ENRICH )
            goterm = strsplit(input$tab_2_go , ' - ')[[1]][1]
            goterm = sub('\\:','\\.' ,goterm)
            fig = plot_go_enrich_violin(df, goterm , paste0(DATADIR,"tab2_go_enrich.rds") )
            fig
        }, error= function(e){})
    })

    # takes a while to load Summarized Experiment package
    output$tab2_go_dot  = plotly::renderPlotly({
        tryCatch({
            goterm = strsplit(input$tab_2_go , ' - ')[[1]][1]
            goterm = sub('\\:','\\.' ,goterm)

            fig = plot_go_enrich_dot( goterm , TAB2_LOOM_FILE_PATH,"FACS_labels",tab='tab2')
            fig
        }, error= function(e){})
       
    })
} # end tab 2


{ # tab 3

    updateSelectizeInput( session , 'tab_3_gene', choices = GENELIST,selected= 'Klf1' )
    # updateSelectizeInput( session , 'tab_3_go', choices = GO_LIST, selected = "GO.0000002" )

    # render plots
    output$tab3_umaps = renderPlot({
        tryCatch({ #gene_input =input$tab_2_gene
            fig = plot_pt_genes_umap(input$tab_3_gene, TAB1_LOOM_FILE_PATH  )
            fig
        }, error=function(e){})
    })


    output$tab3_pt_stats  = DT::renderDataTable({
        tryCatch({
            df= read.delim(paste0(DATADIR ,'pseudotime_statistics.csv') ,sep=",",row.names="X")
            df$fc = format( as.numeric(df$fc),scientific =FALSE ,digits=3)
            df$p_adj = format( as.numeric(df$p_adj),scientific =TRUE ,digits=3)
            df$p_adj_eryth = format(as.numeric(df$p_adj_eryth),scientific =TRUE ,digits=3)
            df$p_adj_mono = format(as.numeric(df$p_adj_mono),scientific =TRUE ,digits=3) 
            
            ind = order(df$p_adj)
            df[ind,]
        },error =function(e){}) 
        
        
        },rownames=FALSE,options = list(pageLength = 10, scrollX = TRUE)
    )
    
    output$tab3_go_enrich  = DT::renderDataTable({
        tryCatch({            
            df = read.delim(paste0(DATADIR ,'pseudotime_go_enrichment.csv') ,sep=",",row.names="X")
            ind = order(df$Padj)
            df$Padj = format(as.numeric(df$Padj),scientific = TRUE, digits=3)
            df[ind,]
        },error =function(e){} ) 
        
        
        },rownames=FALSE,options = list(pageLength = 10, scrollX = TRUE)
    )

} # end tab 3


{
    updateSelectizeInput( session , 'tab_4_mouse_gene', choices = GENELIST,selected= 'Cd74' )
    
    observeEvent(input$tab_4_mouse_gene,{        
        # change the input gene - 
        tryCatch({
            ens=getNetworkID('mouse',input$tab_4_mouse_gene)
            hsgenes = orthologMapper('mouse',ens$found[!ens$notFound] ,'human')[[1]]['human']
            hsnetids = h5read('/home/ftp/data/HSC_atlas/human_scRNAseq_data.loom','row_attrs/NetworkIDs')

            ind = match(hsgenes[,1] , hsnetids )
            ind = ind[!is.na(ind)]
            varnameshs = h5read('/home/ftp/data/HSC_atlas/human_scRNAseq_data.loom','row_attrs/var_names')
            varnameshs = varnameshs[ind]
            if (length(varnameshs) > 0){
                updateSelectizeInput( session , 'tab_4_human_gene', choices = varnameshs,selected =varnameshs[1]  )
            } else {
                updateSelectizeInput( session , 'tab_4_human_gene', choices = 'No Orthologs Available')
            }

            drgenes = orthologMapper('mouse',ens$found[!ens$notFound] ,'zebrafish')[[1]]['zebrafish']
            drnetids = h5read('/home/ftp/data/HSC_atlas/zebrafish_scRNAseq_data.loom','row_attrs/NetworkIDs')

            ind = match(drgenes[,1] , drnetids )
            ind = ind[!is.na(ind)]
            varnamesdr = h5read('/home/ftp/data/HSC_atlas/zebrafish_scRNAseq_data.loom','row_attrs/var_names')
            varnamesdr = varnamesdr[ind]
            if( length(varnamesdr) > 0){
                updateSelectizeInput( session , 'tab_4_zebrafish_gene', choices = varnamesdr,selected =varnamesdr[1]  )
            } else {
                updateSelectizeInput( session , 'tab_4_zebrafish_gene', choices = 'No Orthologs Available')
            }

        }, error= function(e){} )
        
    })

    # if (length(drgenes) > 0 ){
    #     updateSelectizeInput( session , 'tab_4_zebrafish_gene', choices = drgenes,selected = drgenes[1] )
    # }


    # updateSelectizeInput( session , 'tab_4_human_gene', choices = hsgenes )


    output$tab4_mouse_umap = renderPlot({
        tryCatch({ #gene_input =input$tab_2_gene
            fig = plot_genes_on_umap(input$tab_4_mouse_gene, TAB1_LOOM_FILE_PATH)
            fig
        }, error=function(e){})
    })
        
    output$tab4_human_umap = renderPlot({
        tryCatch({ #gene_input =input$tab_2_gene
            fig = plot_genes_on_umap(input$tab_4_human_gene, paste0(DATADIR, 'human_scRNAseq_data.loom'))
            fig
        }, error=function(e){})
    })

    output$tab4_zebrafish_umap = renderPlot({
        tryCatch({ #gene_input =input$tab_2_gene
            fig = plot_genes_on_umap(input$tab_4_zebrafish_gene, paste0(DATADIR, 'zebrafish_scRNAseq_data.loom'))
            fig
        }, error=function(e){})
    })

} # end tab 4


} # end function







