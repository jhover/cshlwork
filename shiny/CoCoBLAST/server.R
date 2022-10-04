##################################################################################################*
############################################ Libraries ############################################
##################################################################################################*
# direc = "/srv/shiny-server/sample-apps/CoCoBLAST/" 

# library(rhdf5) 
# library(Matrix)
# load("/ftp/CoCoBLAST/gene2go/voc.Rdata")           # descrptions of GO terms
# load(paste0(direc,'www/EGADlite.Rdata'))
 
# library(SingleCellExperiment)
library(tidyverse)
# library(glue)
library(ComplexHeatmap)
library(rhdf5,quietly=T)
library(ggplot2)
# library(scattermore)
# library(visNetwork,quietly = T)



noMeta =c('human','atlanticsalmon','cotton','mustard','rainbowtrout','rapeseed','tea','tobacco')
##################################################################################################*
############################################ Functions ############################################
##################################################################################################*
#keep --
getColor <- function(species){
    #' Obtain colors matching the paper
    

    specList = list(
        fungi       = c('yeast','neurospora','fissionyeast','candida'),
        triclads    = c('smed'),
        fish        = c('atlanticsalmon','rainbowtrout','zebrafish','stickleback'),
        frogs       = c('afrog','wfrog'),
        ungulate     = c('dog','horse','cow','sheep','goat','boar'),
        rodents     = c('rabbit','rat','mouse'),
        primates    = c('rhesusm','crabm','chimp','human'),
        inverts     = c('bee','fruitfly','mosquito','bombyx'),
        roundworm   = c('roundworm'),
        algae       = c('chlam'),
        dicots      = c('rice','sorghum','brome','indica','maize'),
        eudicots    = c('grape','cotton','arabidopsis','mustard','rapeseed','cucumber','soybean','medicago','peanut', 'apple', 'tobacco', 'tomato', 'potato','tea'),
        unicell     = c('tryp','plasmo')    
    )

    group =names(specList[unlist(lapply(specList, function(x)  species%in% x ))])

    cols = list(
        fungi       ='#dcb966',
        triclads    ='#7e6072',
        fish        ='#af8c90',
        frogs       ='#01bbbe',
        birds       ='#a2d7c4',
        ungulates   ='#9697cb',
        rodents     ='#4156a6',
        primates    ='#282e58',
        inverts     ='#dd716e',
        roundworm   ='#d46827',
        algae       ='#8cc63f',
        dicots      ='#d7df23',
        eudicots    ='#009444',
        unicell     ='#58585b'
    )

    return(cols[[group]])
}


# keep --
geneConversionTable <- function(species){
    return( read.delim(paste0(dataDirec, "geneInfo/",species,"_info.csv"), stringsAsFactors=F, sep=","))
}

# keep --
specDisplay<-function(spec){
    speciesList = c( 
        "Human"         ="human",
        "Mouse"         ="mouse",
        "Arabidopsis"   ="arabidopsis",                   
        "Cow"           ="cow", 		   
        "C. elegans"    ="roundworm",
        "Zebrafish"     ="zebrafish",
        "Fruitfly"      ="fruitfly",
        "Chicken"       ="chicken",
        "Soybean"       ="soybean",
        "Rat"           ="rat",
        "Budding yeast" ="yeast",
        "Boar"          ="boar",
        "Rice"          ="rice",
        "Maize"         ="maize",
        "Atlantic Salmon"   ='atlanticsalmon',
        'Honey Bee'     ='bee',
        'Silkworm'      ='bombyx',
        'Brome'         ='brome',
        'Chimpanzee'    ='chimp',
        'C. reinhardtii'    ='chlam',
        'Crab-eating macaque'   ='crabm',
        'Dog'   ='dog',
        'Fission yeast' ='fissionyeast',
        'Goat'          ='goat',
        'Grape'         ='grape',
        'Horse'         ='horse',
        'Medicago'      ='medicago',
        'Mosquito'      ='mosquito',
        'Potato'        ='potato',
        'Rabbit'        ='rabbit',
        'Rhesus macaque'    ='rhesusm',
        'Sheep'         ='sheep',
        'Sorghum'       ='sorghum',
        'Tomato'        ='tomato',
        'Western clawed frog'   ='wfrog',
        'Neurospora'    = 'neurospora',
        'Mustard'       = 'mustard',
        'Apple'         = 'apple',
        'Tobacco'       = 'tobacco',
        'Rainbow Trout' = 'rainbowtrout'
    )
    return(names(speciesList[speciesList==spec] ) )
}

# discard
order_rows_according_to_cols = function(M, alpha = 1) {
    M <- M**alpha
    row_score <- colSums(t(M)*seq_len(ncol(M)), na.rm=TRUE)/rowSums(M, na.rm=TRUE)
    return(order(row_score))
}


# toss
cluster_scores_by_ortholog <- function(aurocs, orthoMat, alpha_row=10,alpha_col=1){
    tmpAurocs = aurocs * orthoMat
    # auroc_cols <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
    # breaks <- seq(-5, 5, length=99)
    # breaks = 1/(1+ exp(breaks) )
    # breaks = c(1, breaks, 0)
    # breaks = breaks[length(breaks):1]

    auroc_no_na <- tmpAurocs
    auroc_no_na[is.na(tmpAurocs)] <- 0
    col_order <- stats::as.dendrogram(
        stats::hclust(stats::dist(t(auroc_no_na)**alpha_col), method = "average")
    )

    row_order = order_rows_according_to_cols(auroc_no_na[,labels(col_order) ] ,alpha_row)
    aurocs = aurocs[row_order, labels(col_order)]
    return(aurocs)
}




# keep --
plotHeatMap <- function(aurocs,orthoMap, specA, specB, alpha_col = 1, alpha_row = 10, cex = 1, margins = c(8, 8))  {

    auroc_cols <- rev((grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, 
        "RdYlBu")))(100))

    # for gplots::
    # breaks <- seq(0, 1, length = 101)


    # for ComplexHeatmap::
    breaks <- seq(0, 1, length = 100)
    ht_cols = circlize::colorRamp2(breaks = breaks, colors = auroc_cols)


    # breaks <- seq(-4, 4, length=99)
    # breaks = 1/(1+ exp(breaks) )
    # breaks = c(1, breaks, 0)
    # breaks = breaks[length(breaks):1]
    
    

    # auroc_no_na <- aurocs
    # col_order = 
    # auroc_no_na[is.na(aurocs)] <- 
    # col_order <- stats::as.dendrogram(stats::hclust(stats::dist(t(auroc_no_na)^alpha_col), 
        # method = "average"))
    # row_order <- order_rows_according_to_cols(auroc_no_na[, labels(col_order)], 
        # alpha_row)
    # aurocs <- aurocs[row_order, ]

    inp     = (orthoMap$input_displayGenes)
    para    = (orthoMap[[paste0(specA,'_displayGenes')]])
    ortho   = (orthoMap[[paste0(specB,'_displayGenes')]])


    rr = match(rownames(aurocs),para)
    cc = match(colnames(aurocs),ortho)


    # m = match(pmaps[m,2],unique(pmaps[,2] ))

    
    
    # n = match(getDisplayNames(specA,orthoMap[n,specA])  ,para)
    # n = match(inp[n],unique(inp) )




    
    # sideCols <- RColorBrewer::brewer.pal(8, "Dark2")
    # sideCols = sideCols[1:length(unique(inp))]
    # if (sum(is.na(n))) {
    #     sideCols = c(sideCols, '#FFFFFF')
    #     n[is.na(n)] = length(sideCols)
    # }

    # if (sum(is.na(m))) {
    #     sideCols = c(sideCols, '#FFFFFF')
    #     m[is.na(m)] = length(sideCols)
    # }



    ComplexHeatmap::Heatmap(aurocs, name = "AUROC",row_split = inp[rr],column_split = inp[cc] ,col=ht_cols ,border=TRUE ,height = 200)
    # ComplexHeatmap::draw(ht)
    

    # getDisplayNames(specA,pmaps$input[n])


    # n = match( colnames(aurocs) ,getDisplayNames(specB,orthoMap[,specB] ) )
    # n = match( orthoMap[n,specA],pmaps[,1] )
    # n = match(pmaps[n,2] ,unique(pmaps[,2]))

    # gplots::heatmap.2(x = aurocs, margins = margins, ylab=specDisplay(specA),xlab = specDisplay(specB), key = TRUE, 
    #     keysize = 1, key.xlab = "AUROC", key.title = "NULL", 
    #     offsetRow = 0.1, offsetCol = 0.1, trace = "none", density.info = "none", 
    #     dendrogram = "none", col = auroc_cols, breaks = breaks, 
    #     na.color = '#FFFFFF', Rowv = FALSE, Colv = FALSE, 
    #     cexRow = cex, cexCol = cex,
    #     RowSideColors = sideCols[m], ColSideColors = sideCols[n]
    # ) 

    # legend(x =.2, y =1.125,xpd=TRUE,
    #     title='Input Gene',
    #     legend = unique(inp),
    #     col = sideCols,
    #     lty= 1,
    #     lwd = 5,
    #     cex=1, ncol = ceiling(length(sideCols)/3 )
    # )    
}

# keep --
getNetworkID <- function(species, gene_list,flag='meta',verbose = TRUE){
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

# keep -- 
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


# tossing
orthoGroup <- function(speciesA,genes, speciesB,style="NM"){
    ##  Maps genes in speciesA to orthologous genes in speciesB
    #   Input:  genes           - genes of SpeciesA
    #           speciesA        - name of reference species
    #           speciesB        - name of target species
    #   Output: temp            - orthologous genes in speciesB

    ## get NM orthologs  
    file1  = paste0( dataDirec,'orthologMaps/' ,speciesA,"_",speciesB,"_orthoNM.csv")
    file2  = paste0( dataDirec,'orthologMaps/' ,speciesB,"_",speciesA,"_orthoNM.csv")



    if (file.exists(file1) ) {
        orthoDat = read.table(file1,sep=",",stringsAsFactors=F,header=T)
    } else if(file.exists(file2)) {
        orthoDat = read.table(file2,sep=",",stringsAsFactors=F,header=T)
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
        orthoMat=sparseMatrix(i = m, j=n )
        rownames(orthoMat) =  getDisplayNames(speciesA, unique(tarorthogenes[,speciesA]))
        colnames(orthoMat) =  getDisplayNames(speciesB, unique(tarorthogenes[,speciesB]))
    } else {
        orthoMat= NULL
    }

    return(list(tarorthogenes[,1:2],orthoMat))
}



# keep --
orthologMapper <- function(speciesA,genes, speciesB,style="NM"){
    ##  Maps genes in speciesA to orthologous genes in speciesB
    #   Input:  genes           - genes of SpeciesA
    #           speciesA        - name of reference species
    #           speciesB        - name of target species
    #   Output: temp            - orthologous genes in speciesB

    ## get NM orthologs  
    file1  = paste0( consDirec ,speciesA,"_",speciesB,"_CoCoBLAST_scores.hdf5")
    file2  = paste0( consDirec ,speciesB,"_",speciesA,"_CoCoBLAST_scores.hdf5")

    

    if (file.exists(file1) ) {
        orthoDat = h5read(file1, name = 'ortho_map' ) # read.table(file1,sep=",",stringsAsFactors=F,header=T)
    } else if(file.exists(file2)) {
        orthoDat = h5read(file2, name = 'ortho_map' ) #read.table(file2,sep=",",stringsAsFactors=F,header=T)
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
        orthoMat=sparseMatrix(i = m, j=n )
        rownames(orthoMat) =unique(tarorthogenes[,speciesA]) #  getDisplayNames(speciesA, )
        colnames(orthoMat) =unique(tarorthogenes[,speciesB])  # getDisplayNames(speciesB, unique(tarorthogenes[,speciesB]))
    } else {
        orthoMat= NULL
    }



    return(list(tarorthogenes,orthoMat))
}

# keep --
getParalogs <- function(species,genes){
    #edited to include oringal gene
    fname = paste0(dataDirec, 'paralogMaps/',species,'_paralogMaps.csv')
    pm = read.delim(fname, sep=',',stringsAsFactors=FALSE)


    nm = unique(c(pm[,1],pm[,2]))
    i = match(pm[,1] , nm)
    j = match(pm[,2],nm)

    n = length(nm)

    mat = sparseMatrix(i=i,j=j,dims = c(n,n) , dimnames= list(nm,nm),symmetric=TRUE)

    m = match(genes,  nm )
    m = m[!is.na(m)]


    mat = mat[,m , drop=FALSE]
    mat = mat[rowSums(mat)>0,,drop=FALSE]
        rn = rownames(mat)
    
    
    df = summary(mat)
    df = data.frame(paralog = rownames(mat)[df$i] ,input= colnames(mat)[df$j] ,stringsAsFactors=FALSE)

    df = rbind(df , data.frame(paralog=genes,input=genes,stringsAsFactors=FALSE))
    return(list(c(rn,genes),df))
}

# keep --
getCoexpConsScores <- function(specA, specB ,geneMaps){
    
        file1  = paste0( consDirec ,specA,"_",specB,"_CoCoBLAST_scores.hdf5")
        file2  = paste0( consDirec ,specB,"_",specA,"_CoCoBLAST_scores.hdf5")

        if (file.exists(file1)){

            # get specA genes
            specAgenes = h5read(file1, name = 'sp1_netid')
            specBgenes = h5read(file1, name = 'sp2_netid')

            specAind = match(toupper(geneMaps[,specA]),toupper(specAgenes ))
            specAind= unique(specAind[!is.na(specAind)])
            specBind = match(toupper(geneMaps[,specB]),toupper(specBgenes) )
            specBind= unique(specBind[!is.na(specBind)])

            cons = h5read(file1, name = 'fc_scores',index = list(specAind,specBind) ,drop=FALSE)
            # rownames(cons) = getDisplayNames (specA, specAgenes[specAind] )
            # colnames(cons) = getDisplayNames (specB, specBgenes[specBind] )
            rownames(cons) = specAgenes[specAind]
            colnames(cons) = specBgenes[specBind]
            
            spec = h5read(file1, name = 'sc_scores',index = list(specAind,specBind),drop=FALSE ) 
            rownames(spec) = rownames(cons)
            colnames(spec) = colnames(cons)


        } else if (file.exists(file2) ) {

            # get specA genes
            specAgenes = h5read(file2, name = 'sp2_netid')
            specBgenes = h5read(file2, name = 'sp1_netid')

            specAind = match(toupper(geneMaps[,specA]),toupper(specAgenes ))
            specAind= unique(specAind[!is.na(specAind)])
            specBind = match(toupper(geneMaps[,specB]),toupper(specBgenes ))
            specBind= unique(specBind[!is.na(specBind)])

            cons = h5read(file2, name = 'fc_scores',index = list(specBind,specAind) ,drop=FALSE)
            # colnames(cons) = getDisplayNames (specA, specAgenes[specAind] )
            # rownames(cons) = getDisplayNames (specB, specBgenes[specBind] )
            colnames(cons) = specAgenes[specAind]
            rownames(cons) = specBgenes[specBind]

            cons = t(cons)
            spec = h5read(file2, name = 'sc_scores',index = list(specBind,specAind),drop=FALSE) # all zeros?!
            colnames(spec) = rownames(cons)
            rownames(spec) = colnames(cons)
            spec=t(spec)

        } 

        
    return(list(cons,spec))
}

clusterByGroup <-function(scores,paraMap, orthoMap ){
    # given a scores matrix, reutns a rearranged version 
    return()
}

plot_genes_on_umap <- function (scores, umap_coordinates, normalize_scores = FALSE ) {
    umap = umap_coordinates
    colnames(umap_coordinates) = c("umap_1", "umap_2")
    to_plot =cbind(score= scores ,umap_coordinates)

    result = ggplot2::ggplot() + geom_point(data=to_plot, ggplot2::aes(x = umap_1, 
        y = umap_2, col = score), size=.5) + #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
        ggplot2::scale_color_gradient(low = "gray80", 
        high = "darkred") + ggplot2::theme_classic() +
        ggplot2::labs(x = colnames(umap)[1], y =  colnames(umap)[2], col = "log(CPM + 1)") 
        


        # ggplot2::annotate("text", x = 4, y = 25, label = "Some text")
    
        

    
    return(result)
}


load_tSNE_data <- function(species,genes){ # allow user to input TM or arabidopsis
    # geneList should be input genes converted to either mouse or arabidopsis

    # load the data and 
    if (species == 'mouse'){ # compare to tm
        filepath = paste0(consDirec ,'../TM_data.hdf5')
    } else if (species =='arabidopsis' ) { # compare to arabidopsis
        filepath = paste0(consDirec ,'../Ara_data.hdf5')
    } 



    geneDat = h5read(filepath,name = 'rowdata')
    m = match(genes , geneDat$symbol)
    m= m[!is.na(m)]
    
     
    cellDat =  h5read(filepath,name = 'coldata')
    
    counts = h5read(filepath, name = 'normalized_counts',index =list(m,1:nrow(cellDat)))
    rownames(counts) = geneDat$symbol[m]
    colnames(counts) = cellDat$cell

    embedding = h5read(filepath, 'embedding')

    if (species =='arabidopsis'){ 
        # need to remove cells with obsurd umap coord
        ind=which(abs(embedding[,2]) > 16 ) 
        cellDat = cellDat[-ind,]
        counts = counts[,-ind]
        embedding = embedding[-ind,]

        # color and re name columns for now - eventually, clusters would be annotated
        tissues = c('Atrichoblast','Trichoblast','Xylem','Phloem','Cortex','Meristematic','Endodermis','Phloem','Endodermis')
        # tissues = 1:9
        # green ->[6,5,1]
        # reds/brown -> [4,2,7]
        # purple -> [3]
        # 8
        colors = RColorBrewer::brewer.pal(8,'Dark2')[c(1,5, 2,7,3, 6 , 8,7,8)]
        cellDat$tissue = tissues[cellDat$metaclust]
        cellDat$color = colors[cellDat$metaclust]
        

    }

    
    return(list(counts=counts,embedding=embedding,cellDat=cellDat))
}


getTopRankedGenes <- function(speciesA, genes, speciesB,topN = 50 ) {
    # input should be network IDs of the input genes
    # genes = as.list(genes)


    file1  = paste0( consDirec ,speciesA,"_",speciesB,"_CoCoBLAST_scores.hdf5")
    file2  = paste0( consDirec ,speciesB,"_",speciesA,"_CoCoBLAST_scores.hdf5")

    if (file.exists(file1)){
        specAgenes = h5read(file1, name = 'sp1_netid')
        m = match(genes, specAgenes)
        m = m[!is.na(m)]

        specBgenes = h5read(file1, name = 'sp2_netid')
        cons = h5read(file1, name = 'fc_scores', index =list(m,1:length(specBgenes) ) )
        colnames(cons) = specBgenes
        rownames(cons) = specAgenes[m]
        
        df = apply(cons,1, function(x){
            ind = order(x,decreasing=TRUE)
            a = names(x[ind]) [ 1:topN ]
            b = x[ind][1:topN]
            return(cbind(a,b))
        })

        # df = data.frame(a=df[1:topN,], b=df[(topN+1):(2*topN),] )
        df = data.frame(gene=df[1:topN,], score=df[(topN+1):(2*topN),] ,stringsAsFactors=FALSE )
        colnames(df)[1:(length(genes)) ] = genes
        colnames(df )[(length(genes)+1): (2*length(genes)) ] = paste0(getDisplayNames(speciesA, genes) ,'_coexp_cons')
        df = rbind(colnames(df), df )


        # recip ranks
        df2 = apply(df[,1:length(genes),drop=FALSE],2 ,function(x){
            gene = x[1]
            m = match(x[2:length(x)], specBgenes)
            m = m[!is.na(m)]

            
            cons2 = h5read(file1, name = 'fc_scores', index =list( 1:length(specAgenes),m ) )
            colnames(cons2) = specBgenes[m]
            rownames(cons2) = specAgenes 
            temp= apply(cons2, 2, rank)[gene,]
            nrow(cons2) - temp +1

        })
        colnames(df2) = paste0(getDisplayNames(speciesA, genes),'_recip_ranks' )

        df = cbind(df[2:nrow(df),] , df2)
        t=data.frame(apply( df[,1:length(genes),drop=FALSE ] , 2,getDisplayNames, species = speciesB ))
        colnames(t) = getDisplayNames(speciesA, colnames(t) )

        df = cbind(df ,t)



        df = cbind(Rank = 1:nrow(df),df)
        return(df)
    } else if (file.exists(file2) ){
        specAgenes = h5read(file2, name = 'sp2_netid')
        m = match(genes, specAgenes)
        m = m[!is.na(m)]

        specBgenes = h5read(file2, name = 'sp1_netid')
        cons = h5read(file2, name = 'fc_scores', index =list(1:length(specBgenes), m )  ,drop=FALSE)
        rownames(cons) = specBgenes
        colnames(cons) = specAgenes[m]
        
        df = apply(cons,2, function(x){
            ind = order(x,decreasing=TRUE)
            a = names(x[ind]) [ 1:topN ]
            b = x[ind][1:topN]
            return(cbind(gene = a,score =b))
        })

        df = data.frame(gene=df[1:topN,], score=df[(topN+1):(2*topN),] ,stringsAsFactors=FALSE )
        # df$rank= 1:topN

        colnames(df)[1:(length(genes)) ] = genes
        colnames(df )[(length(genes)+1): (2*length(genes)) ] = paste0(getDisplayNames(speciesA, genes) ,'_coexp_cons')
        df = rbind(colnames(df), df )

        df2 = apply(df[,1:length(genes),drop=FALSE],2 ,function(x){
            gene = x[1]
            m = match(x[2:length(x)], specBgenes)
            m = m[!is.na(m)]

            
            cons2 = h5read(file2, name = 'fc_scores', index =list(m, 1:length(specAgenes) ) )
            rownames(cons2) = specBgenes[m]
            colnames(cons2) = specAgenes 
            temp= apply(cons2, 1, rank)[gene,]
            ncol(cons2) - temp +1

        })
        colnames(df2) = paste0(getDisplayNames(speciesA, genes),'_recip_ranks' )

        df = cbind(df[2:nrow(df),] , df2)
        t=data.frame(apply( df[,1:length(genes),drop=FALSE ] , 2,getDisplayNames, species = speciesB ))
        colnames(t) = getDisplayNames(speciesA, colnames(t) )

        df = cbind(df ,t)



        df = cbind(Rank = 1:nrow(df),df)
        return(df)
        
    }   

}


##################################################################################################*
############################################# Driver ##############################################
##################################################################################################*

# shiny_env=new.env()

function(input,output,session){
{ # home page
    # observeEvent(input$toComparison,{
        
    #     updateTabItems(session,'tabs','coexnet')
    # })
    observeEvent(input$toConservation,{
        
        updateTabItems(session,'tabs','cons')


    })
 

}


# conservaton score page
{
    
    # consSpecA = reactiveVal()
    # consSpecB = reactiveVal()
    # consGenes = reactiveVal() 


    getResultsCons<- function(speciesA, speciesB, geneList){
        
        # verify that speciesA and species B are different
        if (speciesA == speciesB){ 
            showModal(modalDialog(title='Uh Oh... ','Please choose two different species!',fade=TRUE ))
                Sys.sleep(5)
            removeModal()
            return()
        }

        # clean up gene list
        geneList = trimws(strsplit(geneList,',')[[1]],which='both')
        

        # given the list of genes, need to convert them to the ENS IDs, extract paralogs, get orthologs, get scores between orthologs, display plots
        geneDF = getNetworkID( speciesA, geneList ,flag= 'meta', verbose= FALSE  )
        geneList = unique(geneDF$found[!geneDF$notFound  ])

        # need to verify that there are genes found - otherwise quit.
        if (length(geneList) == 0 ) { # no genes found  - do nothing
            showModal(modalDialog(title='Uh Oh... ',"We didn't find any of your input genes! Try other genes or contact John Lee at johlee@cshl.edu for help.",fade=TRUE ))
                Sys.sleep(5)
            removeModal()
            return()
        }


        # look for para/orthologs

        # given the list of genes, need to convert them to the ENS IDs, extract paralogs, get orthologs, get scores between orthologs, display plots
        geneDF = getNetworkID( speciesA, geneList ,flag= 'meta', verbose= FALSE  )
        geneList = unique(geneDF$found[!geneDF$notFound  ])

        # get orthologs
        tempOM = orthologMapper(speciesA, geneList,speciesB)
        
        # do we have orthologs???
        if (nrow(tempOM[[1]]) ==0 ){
            showModal(modalDialog(title='Uh Oh... ',paste0('None of your input genes have orthologs in ', speciesB,'. Try a different species, other genes, or to report an issue, contact John Lee at johlee@cshl.edu ' ),fade=TRUE ))
                Sys.sleep(5)
            removeModal()
            return()      
        }

        orthoMaps = orthologMapper(speciesB, tempOM[[1]][,speciesB],speciesA) 
        tempOM[[1]] = tempOM[[1]][,c(speciesB,speciesA)]
        colnames(tempOM[[1]])  = c(speciesB, 'input_gene')
        orthoMat = t(orthoMaps[[2]])

        orthoMap = orthoMaps[[1]][,c(speciesA,speciesB)]

        # verify that we get orthologs 
        if (nrow(orthoMap) ==0 ){
            showModal(modalDialog(title='Uh Oh... ',paste0('None of your input genes have orthologs in ', speciesB,'. Try a different species, other genes, or to report an issue, contact John Lee at johlee@cshl.edu ' ),fade=TRUE ))
                Sys.sleep(5)
            removeModal()
            return()      
        }


        orthoMap = merge(tempOM[[1]],orthoMap, by.x = speciesB)
        tmp = data.frame(a = getDisplayNames(speciesA,orthoMap[,'input_gene'] ), b=getDisplayNames(speciesA,orthoMap[,speciesA] ) ,c= getDisplayNames(speciesB,orthoMap[,speciesB] ) )
        colnames(tmp) =  paste0(c('input',speciesA,speciesB), '_displayGenes')
        orthoMap = cbind(orthoMap, tmp)


        inputGenes =unique(orthoMap[,'input_gene'] ) 
        paraGenes = unique(orthoMap[,speciesA] )
        orthoGenes = unique(orthoMap[,speciesB] )
        inputSyms =unique(orthoMap[,'input_displayGenes'])
        paraSyms = unique(orthoMap[,paste0(speciesA,'_displayGenes')])
        orthoSyms =unique(orthoMap[,paste0(speciesB,'_displayGenes')])



        # map genes to mouse and arabidopsis as well.
        tsnegenes=list()
        if (speciesA == 'mouse'){ 
            tsnegenes[['Tabula Muris']] =  paraSyms
            tsnegenes[['Arabidopsis Root']] = getDisplayNames('arabidopsis', orthologMapper(speciesA,inputGenes,'arabidopsis')[[1]][,1])
        } else if (speciesA =='arabidopsis'){
            tsnegenes[['Tabula Muris']] = getDisplayNames('mouse', orthologMapper(speciesA,inputGenes,'mouse')[[1]][,1])
            tsnegenes[['Arabidopsis Root']] = paraSyms
        } else {
            tsnegenes[['Tabula Muris']] = getDisplayNames('mouse',orthologMapper(speciesA,inputGenes,'mouse')[[1]][,1])
            tsnegenes[['Arabidopsis Root']] = getDisplayNames('arabidopsis',orthologMapper(speciesA,inputGenes,'arabidopsis')[[1]][,1])
        }
        tsnegenes = lapply(tsnegenes, unique)


        # load tsne data
        tsneDat=list()
        tsneDat[['Tabula Muris']]= load_tSNE_data('mouse',tsnegenes[['Tabula Muris']])      # displa names are given
        tsneDat[['Arabidopsis Root']]= load_tSNE_data('arabidopsis',tsnegenes[['Arabidopsis Root']])
  
        tsnegenes[['Tabula Muris']] = rownames(tsneDat[['Tabula Muris']]$counts )
        tsnegenes[['Arabidopsis Root']] = rownames(tsneDat[['Arabidopsis Root']]$counts )


        # get average expression for each gene in each tissue
        # tsneMeanExp =list()
        for (dataSet in c('Arabidopsis Root','Tabula Muris')) {

            tmp = tsneDat[[dataSet]]
            tissues = unique(tmp$cellDat$tissue)

            df = lapply(as.list(tissues) , 
                function(tissue ){

                    genes = tsnegenes[[dataSet]]
                    cellInd = tmp$cellDat$tissue == tissue
                    color = tmp$cellDat$color[cellInd]

                    df=  data.frame(tissue = rowMeans(tmp$counts[,cellInd,drop=FALSE]))
                    colnames(df) = tissue
                
                    df = as.data.frame(t(df) )
                    df$Color = color[1]
                    df

                }
            )

            tsneDat[[dataSet]]$meanExpression = do.call(rbind,df )

        }



        aucs = getCoexpConsScores(speciesA, speciesB, orthoMap)
        consMat = aucs[[1]]
        specMat = aucs[[2]]


        rr = match(rownames(orthoMat),rownames(consMat)  )
        cc = match(colnames(orthoMat),colnames(consMat)  )
        consMat = consMat[rr,cc,drop=FALSE]
        specMat = specMat[rr,cc,drop=FALSE]



        rs = match(toupper(rownames(consMat)) , toupper(orthoMap[[speciesA]]) )
        cs = match(toupper(colnames(consMat)) , toupper(orthoMap[[speciesB]]) )
        rownames(consMat) = orthoMap[[paste0(speciesA,'_displayGenes')]][rs]
        colnames(consMat) = orthoMap[[paste0(speciesB,'_displayGenes')]][cs]

        rownames(orthoMat) = rownames(consMat)
        colnames(orthoMat) = colnames(consMat)

        rownames(specMat) = rownames(consMat)
        colnames(specMat) = colnames(consMat)


        
        consFlat = consMat * orthoMat
        specFlat = specMat * orthoMat
        ind = which(consFlat > 0,arr.ind=TRUE)
        # ind2 = which(specFlat > 0,arr.ind=TRUE)

        # should not lose precision here. Allow user download
        tableDF = data.frame(geneA = rownames(consFlat)[ind[,1]], geneB = colnames(consFlat)[ind[,2]],  cons = consFlat[ind ], spec = specFlat[ind] )
        colnames(tableDF) = c(specDisplay(speciesA), specDisplay(speciesB), 'Coexp conservation','Specificity')

        m = match(tableDF[,1] , paraSyms)
        tableDF['inputGene'] = orthoMap[m,'input_displayGenes']

        # score table
        tmptable = tableDF
        m = match(tmptable[,1] , orthoMap[[paste0(speciesA,'_displayGenes')]])
        tmptable[,3:4] = round(tmptable[,3:4],4)
        
        topDF = getTopRankedGenes(speciesA, inputGenes, speciesB )

        # tryCatch({ # cluster all data throwing away off diagonal blocks
        #     consMat=cluster_scores_by_ortholog(consMat, orthoMat)
        #     specMat=cluster_scores_by_ortholog(specMat, orthoMat)

            
        # },error=function(e){})


        

        output$viewText11 <- renderUI({h4('View')})
        output$viewText22 <- renderUI({h4('for')})

        output$geneSelector1 <- renderUI({
            # drop down
            # if (length(inputGenes )>1 ){
            tryCatch({

                pickerInput('inputgene',label=NULL, 
                    choices =  inputSyms,
                    selected= inputSyms[1] ,
                    multiple=FALSE,
                    options = list(
                        'actions-box' = FALSE,
                        'none-selected-text'='Choose at least one'
                    )

                )
            },error= function(e){})
            
        })



        output$tableSelector <- renderUI({
            # drop down
            tryCatch({
                pickerInput('table',label=NULL, 
                    choices =c('Top hits','Ortholog Groups') ,
                    selected='Top hits'
                )
            },error=function(e){})
        })



        # render table  

        output$cons_table <- DT::renderDataTable({
            tryCatch({
                if (input$table =='Top hits'){
                    
                    tmpdf= topDF[,c(1,grep(input$inputgene, colnames(topDF) )) , drop=FALSE]


                    
                    tmpdf = tmpdf[,c(4,1,3)]
                    tmpdf[,2] = as.numeric(tmpdf[,2])
                    tmpdf[,3] = as.numeric(tmpdf[,3])
                    colnames(tmpdf) = c(paste(specDisplay( speciesB),'genes') , 'Rank',  'Reciprocal Rank') 
                    tmpdf['Mutual Rank'] = round(sqrt(as.numeric(tmpdf$Rank) * as.numeric(tmpdf$'Reciprocal Rank')),2)
                    # colnames(tmpdf) = c('Rank',paste(specDisplay( speciesB),'genes'),  'Reciprocal Rank' ,'Coexp cons')
                    tmpdf
                } else if (input$table =='Ortholog Groups'){
                    tmptable[tmptable$inputGene == input$inputgene,1:4]
                }
            }, error=function(e){}
            )
        },
            rownames=FALSE,options = list(pageLength = 6)
        )



        # render user input view
        output$viewText1 <- renderUI({h4('View')})



        output$geneSelector <- renderUI({
            # drop down
            tryCatch({
                if(input$fig =='Embeddings'){
                    

                    pickerInput('gene',label=NULL, 
                        choices =  c('Tissues',sort(tsnegenes[[input$tsneSpec]])) ,
                        selected= 'Tissues' ,
                        multiple=FALSE,
                        options = list(
                            'actions-box' = TRUE,
                            'none-selected-text'='Choose at least one'
                        )

                    )
                } else {
                    pickerInput('gene',label=NULL, 
                        choices = inputSyms, # getDisplayNames(speciesA, unique(pmaps$input)) ,
                        selected=  inputSyms, #getDisplayNames(speciesA, unique(pmaps$input)),
                        multiple=TRUE,
                        options = list(
                            'actions-box' = TRUE,
                            'none-selected-text'='Choose at least one'
                        )
                    )
                }
            },error=function(e){}
            )

        })



        output$figSelector <- renderUI({
            # drop down
            tryCatch({
                pickerInput('fig',label=NULL, 
                    choices =c('Co-exp Conservation Scores','Specificity Scores','Embeddings') ,
                    selected='Co-exp Conservation Scores'
                )
            },error=function(e){})
        })


        output$viewText2 <- renderUI({
            tryCatch({

                h4('for')

            },error=function(e){})
        })



        output$tnseSpecSelector <- renderUI({
            # drop down
            tryCatch({
                if(input$fig =='Embeddings'){
                    pickerInput('tsneSpec',label=NULL, 
                        choices =  c('Tabula Muris','Arabidopsis Root'),
                        selected= 'Tabula Muris',
                        multiple=FALSE,
                        options = list(
                            'actions-box' = FALSE,
                            'none-selected-text'='Choose at least one'
                        )

                    )
                }
            },error=function(e){}
            )

        })

        output$viewText3 <- renderUI({
            tryCatch({
                if (input$fig !='Embeddings' ){
                    h4('paralogs and orthologs.')
                } else{
                    h4('in:')
                }
            },error=function(e){})

        })
        


        output$mainFig = renderPlot({
            tryCatch({ 
                if (input$fig =='Embeddings'){
                    
                    tsneDF = data.frame(tsneDat[[input$tsneSpec]]$embedding[,1:2], color=tsneDat[[input$tsneSpec]]$cellDat$color, tissue=tsneDat[[input$tsneSpec]]$cellDat$tissue )
                    colnames(tsneDF) = c('tSNE_1','tSNE_2','color','tissue')
                    if (input$gene =='Tissues'){
                        # tsneDF = data.frame(tsneDat[[input$tsneSpec]]$embedding[,1:2], color=tsneDat[[input$tsneSpec]]$cellDat$color, tissue=tsneDat[[input$tsneSpec]]$cellDat$tissue )
                        # colnames(tsneDF) = c('tSNE_1','tSNE_2','color','tissue')
                    


                        ggplot(tsneDF,aes(x=tSNE_1,y=tSNE_2,col=color))+
                            geom_point(size=0.5)+
                            scale_color_identity(breaks = tsneDF$color,labels = tsneDF$tissue,guide = "legend")+
                            guides(colour = guide_legend(title = 'Tissue',override.aes = list(size=8)))+theme_classic() 

                    

                    } else{
                        score=tsneDat[[input$tsneSpec]][['counts']] [input$gene,]
                        umap =tsneDat[[input$tsneSpec]][['embedding']][,1:2]


                        # # get centroids
                        # tisList=as.list(unique(tsneDF$tissue))
                        # A = lapply(tisList, 
                        #     function(x) {
                        #         tmpdf = tsneDF[tsneDF$tissue == x ,]
                        #         c(median(tmpdf[,1]) , median(tmpdf[,2]))
                        #     }
                        # )
                        # A = as.data.frame(do.call(rbind,A))
                        # rownames(A) = unlist(tisList)
                        # colnames(A) = colnames(tmpdf)[1:2]
                        # A$Tissue = rownames(A)
                        

                        plot_genes_on_umap(score,umap, FALSE  )#, annots = A   )
                        
                    }             
                    

                }else { # not embeddings
                    if (input$fig =='Co-exp Conservation Scores'){
                        tmpCons = consMat
                    } else if (input$fig == 'Specificity Scores'){  
                        tmpCons = specMat 
                    } 
                    tmpDF = orthoMap
                    tmpDF = tmpDF[tmpDF$input_displayGenes %in% input$gene , ]
                    # rows
                    m = sort(match(unique(tmpDF[,paste0(speciesA,'_displayGenes')]) ,  rownames(tmpCons)))
                    # cols
                    n =sort(match(unique(tmpDF[,paste0(speciesB,'_displayGenes')]) ,  colnames(tmpCons)))
                    
                    tmpCons = tmpCons[m,n,drop=FALSE]
            
                    ComplexHeatmap::draw(plotHeatMap(tmpCons,orthoMap ,speciesA,speciesB))
                    # ht_pos <<- ht_pos_on_device(shiny_env$ht)    
                    
                    
                }

            }, error = function(e){})
        })

      
        output$bipart_slider <- renderUI({     
            tryCatch({
                if (input$fig !='Embeddings'){
                    if (input$fig =='Co-exp Conservation Scores'){
                        name = 'Coexp conservation'
                    } else {
                        name ='Specificity'
                    }


                    if (nrow(tmptable) <= 2){ 
                        v = round(min(tmptable[,name]),2 ) -0.05
                    } else if (nrow(tmptable) <= 5){
                        v = round(quantile(tmptable[,name],.25 ),2)
                    } else if(nrow(tmptable) <= 10){
                        v = round(quantile(tmptable[,name],.4 ),2 )
                    } else if(nrow(tmptable) <= 20) {
                        v =round(quantile(tmptable[,name],.6) ,2)
                    } else if(nrow(tmptable) <= 50) {
                        v =round(quantile(tmptable[,name],.8) ,2)
                    }else {
                        v = round(quantile(tmptable[,name],.95) ,2)
                    }

                    sliderTextInput("thresInput","Score Threshold:",choices=0:99/100, selected = v,hide_min_max=T )
                } else {
                    if (input$gene !='Tissues'){
                        tmp= round(tsneDat[[input$tsneSpec]]$meanExpression[,input$gene],2)
                        
                        ch = round(seq(from = min(tmp), to = max(tmp) , length.out=51),2)
                        v = round(median(tmp),2)
                        

                        sliderTextInput("thresInput","Expression Threshold:",choices=ch, selected = ch[which.min(abs(ch - v))],hide_min_max=T )
                    }
                }
                


            },error=function(e){})
        })

        # output$spec_heatmap = renderPlot({
        #     tryCatch({ 

        #         # MetaNeighbor::plotHeatmapPretrained(as.matrix(aucMat))
        #         plotHeatMap(as.matrix(specMat),pmaps,orthoMap,speciesA,speciesB )
        #     }, error = function(e){})
        # })
       


        # plot a bipartite graph
        output$bipartite <- renderVisNetwork({

            tryCatch( {
                if (input$fig !='Embeddings'){
                    df = tmptable[tmptable$inputGene %in% input$gene,,drop=FALSE]
                    df2 = data.frame(  colnames(df)[1],colnames(df)[2], 1 ,1,NA  )
                    colnames(df2) =colnames(df)
                    df = rbind( df2, df )
                    n = length(unique(df[,specDisplay(speciesA)])) 
                    m = length(unique(df[,specDisplay(speciesB)])) 
                    # shinyalert(c(n,m))
                    inc = matrix(0, n,m)
                    rownames(inc ) = paste0(unique(df[,specDisplay(speciesA)]) ,speciesA )
                    colnames(inc) =  paste0(unique(df[,specDisplay(speciesB)]) ,speciesB)


                    FCind = as.matrix(df[,1:2 ])
                    FCind[,1] = paste0(FCind[,1],speciesA)
                    FCind[,2] = paste0(FCind[,2],speciesB)
                    # shinyalert(FCind)
                    if (input$fig =='Co-exp Conservation Scores'){
                        colmn='Coexp conservation'
                        # inc[FCind ] = df[,'Coexp conservation']
                    } else if(input$fig=='Specificity Scores'){
                        colmn='Specificity'
                        # inc[FCind ] = df[,'Specificity']
                    }

                    # remove edges with score < thres
                    edge_vals =df[df[,colmn ]>=input$thresInput,colmn ,drop=TRUE]

                    inc[FCind] = as.numeric(df[[colmn]])
                    # filter by threshold
                    inc[inc < input$thresInput] = 0
                    inc = inc[rowSums(inc) >0, colSums(inc) >0 ,drop=FALSE]


                    g = igraph::graph_from_incidence_matrix(inc)
                    # lays = igraph::layout_as_bipartite(g)[,c(2,1)]
                    # lays[,1] = abs(lays[,1]-1)  # make it vertical

                    E(g)$width=10*edge_vals **5
                    E(g)$width[1] =0
                    E(g)$color = c('white', rep('black',length(edge_vals)-1 ))
                    E(g)$label = edge_vals
                    E(g)$label[1] =NA
                    V(g)$color=c(rep(getColor(speciesA),nrow(inc)) , rep(getColor(speciesB),ncol(inc)))
                    # V(g)$id = 1:(nrow(inc)*ncol(inc))
                    # V(g)$label=c(rinc ,cinc)
                    # V(g)$size = size[temp] 

                    visnet = visNetwork::visIgraph(g,layout = 'layout_with_sugiyama' )

                    coord_y <- visnet$x$nodes$y  
                    visnet$x$nodes$y <- visnet$x$nodes$x
                    visnet$x$nodes$x <- coord_y 
                    visnet$x$edges$font.align = 'bottom'
                    visnet$x$edges$font.size = 30
                    visnet$x$nodes$font.size = 30
                    visnet$x$nodes$shape = 'dot'
                    visnet$x$nodes$label = sub(speciesA, '',visnet$x$nodes$label)
                    visnet$x$nodes$label = sub(speciesB, '',visnet$x$nodes$label)
                    

                    
                    m = match(colnames(df)[1:2] ,visnet$x$nodes$label ) 
                    visnet$x$nodes$shape[m] = 'text'

                    # visnet$x$nodes$shape = 'ellipse'

                    # create a white face with black border
                    # visnet$x$nodes$font.color= 'white'
                    # visnet$x$nodes$font.strokeWidth = 3
                    # visnet$x$nodes$font.strokeColor = 'black'


                    visnet %>% 
                        visOptions( highlightNearest = list(enabled=T,hover=T,degree=1)  ) 
                        # visInteraction(zoomView=FALSE,navigationButtons=TRUE) 
                } else {    # is embeddings - show average expression
                    ds = input$tsneSpec
                    tgene = input$gene
                
                    if (tgene != 'tissue'){
                        tmp =tsneDat[[ds]]$meanExpression[, c(tgene,'Color') ,]
                        tsnedf = data.frame(Tissue = rownames(tmp) , Mean_Expression = tmp[,tgene],Color = tmp$Color,stringsAsFactors=FALSE)
                        tsnedf = rbind(colnames(tsnedf), tsnedf)
                        minTmp = min(tmp[,1])
                        maxTmp = max(tmp[,1])
                        

                        n = nrow(tsnedf)
                        
                        # shinyalert(c(n,m))
                        inc = matrix(0, n,n)
                        rownames(inc ) = tsnedf[['Tissue']]
                        colnames(inc) =  tsnedf[['Mean_Expression']]

                        FCind = as.matrix(tsnedf[,1:3 ])
                        ord =order(FCind[,2],decreasing=TRUE)
                      

                        # remove edges with score < thres
                        inc[FCind[,1:2]] = as.numeric(c(100,tsnedf[['Mean_Expression']][2:length(tsnedf[['Mean_Expression']])]))
                        FCind = FCind[ord ,]
                        inc = inc[ord,ord] 
                        # filter by threshold
                        inc[inc < input$thresInput] = 0
                        inc = inc[rowSums(inc) >0, colSums(inc) >0 ,drop=FALSE]
                        FCind = FCind[which(FCind[,2] >= input$thresInput) ,]
                        
                        
                        g = igraph::graph_from_incidence_matrix(inc)
                        # lays = igraph::layout_as_bipartite(g)[,c(2,1)]
                        # lays[,1] = abs(lays[,1]-1)  # make it vertical

                        E(g)$width= 1
                        E(g)$width[1] =0
                        E(g)$color = c('white', rep('black',nrow(FCind)-1 ))
                        # E(g)$label = edge_vals
                        # E(g)$label[1] =NA

                        V(g)$color= c(FCind[,3],FCind[,3])
                        V(g)$color[1]='black'

                        # V(g)$id = 1:(nrow(inc)*ncol(inc))
                        # V(g)$label=c(rinc ,cinc
                        # V(g)$size = size[temp] 

                        visnet = visNetwork::visIgraph(g,layout = 'layout_with_sugiyama' )

                        coord_y <- visnet$x$nodes$y /2
                        visnet$x$nodes$y <- visnet$x$nodes$x 
                        visnet$x$nodes$x <- coord_y 
                        visnet$x$edges$font.align = 'bottom'
                        visnet$x$edges$font.size = 30
                        visnet$x$edges$dashes = TRUE
                        visnet$x$edges$dashes = TRUE
                        visnet$x$nodes$font.size = 30
                        visnet$x$nodes$shape = 'dot'
                        
                        m = match(colnames(tsnedf)[1:2] ,visnet$x$nodes$id ) 
                        visnet$x$nodes$shape[m] = 'text'
                        visnet$x$nodes$shape[grepl('[0-9]',visnet$x$nodes$id )] = 'text'
                        visnet$x$nodes$label=gsub('_',' ', visnet$x$nodes$label)
                        # visnet$x$nodes$title=gsub('_',' ', visnet$x$nodes$label)

                        m = match(FCind[,1],visnet$x$nodes$id )
                        visnet$x$nodes$size = 10
                        visnet$x$nodes$size[m[!is.na(m)]] = 10 + 10*( (as.numeric(FCind[,2]) - minTmp ) /  (0.0001 + maxTmp-minTmp) )
                        visnet$x$nodes$size[is.na(visnet$x$nodes$size )] =max(visnet$x$nodes$size,na.rm=T) +5
                        visnet$x$nodes$label[grepl('[0-9]',visnet$x$nodes$label )]  = round(as.numeric( visnet$x$nodes$label[grepl('[0-9]',visnet$x$nodes$label )]) ,4)
                        # showModal( modalDialog( paste(diag(inc ) ,collapse=', ') ) )
                        # visnet$x$nodes$color.border = 'black'
                        # visnet$x$nodes$font.color = visnet$x$color
                        

                        # need to sort nodes by value and adjust size to be proportional to value.
                        # occasional missing edge

                        visnet %>% 
                            visOptions( highlightNearest = list(enabled=T,hover=T,degree=1)  ) 
                            







                    } else {
                        # plot nothing
                    }
                }
                    


                



            }, error=function(e){}
            )
        })


        # render downloadables
        
        # download selection radio buttons
        output$option <- renderUI({
            # if(length(input$downloadCheck1)==0){}
            
            prettyRadioButtons(
                inputId = "downloadOption",
                label = "Download:",
                thick = TRUE,
                choiceNames = c('Top Hits','Ortholog Maps','Co-expression Conservation Scores','Specificity Scores'),
                choiceValues = c(0,1,2,3),
                animation = "pulse",
                status = "info",outline = T
            )
        
        })

        # button to download the top connections of a gene
        output$downloadButton <- renderUI({
            tryCatch({
                shinyWidgets::downloadBttn('downloader',
                    label='Download:',color='primary',
                    style='stretch',
                    size='sm',block=T)  
            }, error = function(e){})
        })

        # download executor
        output$downloader <- downloadHandler(
            filename = function() {
                paste('CoCoBLAST.csv', sep='')
            },
            content = function(con) {
                if(length(input$downloadOption)==0){}
                else{switch(input$downloadOption,
                    '0' = {
                        tmp = topDF
                        
                        pt1 = inputSyms
                        pt2 = paste0(inputSyms,'_recip_ranks')
                        tmp = tmp[,c('Rank',sort(c(pt1,pt2) )) ]
                        
                        write.table(tmp, con , row.names = F,sep=',') 
                    },

                    '1' = {
                        tmp = orthoMap
                        tmp = tmp[,c(speciesA, speciesB,paste0(c(speciesA, speciesB),'_displayGenes' ))]
                        spA = specDisplay(speciesA)
                        spB = specDisplay(speciesB)
                        colnames(tmp) = c(spA, spB,paste0(c(spA, spB),'_Symbol' ))
                        


                        write.table(tmp, con , row.names = F,sep=',') 
                    },
                    '2' = {
                        tmp = consMat
                        tmp = cbind(rownames(tmp) , tmp)
                        write.table(tmp, con , row.names = F,sep=',') 
                    },
                    '3' = {
                        tmp = specMat
                        tmp = cbind(rownames(tmp) , tmp)
                        write.table(tmp, con , row.names = F,sep=',') 
                    }

                )}
            } 
        )



    }

    clearResults <- function(){
        # output$cons_heatmap <- renderPlot({})
        # output$cons_table <- DT::renderDataTable({})
        # # download selection radio buttons
        # output$option <- renderUI({})

        # # button to download the top connections of a gene
        # output$downloadButton <- renderUI({})

        # # download executor
        # output$downloader <- downloadHandler()

        # output$bipartite <- renderVisNetwork({})

        # output$cons_heatmap <- renderPlot({})


        output$viewText11   <- renderUI({})
        output$viewText22   <- renderUI({})
        output$viewText2   <- renderUI({})

        output$geneSelector1   <- renderUI({})
        output$tableSelector   <- renderUI({})
        output$cons_table   <- DT::renderDataTable({})
        output$viewText1   <- renderUI({})
        output$figSelector   <- renderUI({})
        output$geneSelector    <- renderUI({})
        output$viewText3   <- renderUI({})
        output$tnseSpecSelector   <- renderUI({})
        output$mainFig   <- renderPlot({})
        output$bipart_slider   <- renderUI({})
        # spec_heatmap   <- renderPlot({})
        output$bipartite   <- renderVisNetwork({})
        output$option   <- renderUI({})
        output$downloadButton   <- renderUI({})
        output$downloader <- downloadHandler()



    }

    

    observeEvent(input$cons_submit,
        tryCatch({

            clearResults()
            showModal(modalDialog(title = "Querying results" ,"Please wait.", footer=NULL))
                tryCatch({
                getResultsCons(input$cons_specA,input$cons_specB, input$cons_input) 
                },error=function(e) {})
            removeModal()
        },error=function(e){})

    )


  

}

}







