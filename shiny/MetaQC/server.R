##################################################################################################*
############################################ Libraries ############################################
##################################################################################################*
# direc = "/srv/shiny-server/sample-apps/CoCoBLAST/" 

library(rhdf5) 
library(Matrix)  
library(visNetwork)
library(dplyr)
 
#  class_label class_color

### TODO 
# 10x / ss cell counts
# Filter cells based on histograms
# PR curves in the adata (statistics.py)




##################################################################################################*
############################################ Functions ############################################
##################################################################################################*
savefig <- function(fig){
    ggplot2::ggsave('~/test.png',plot=fig)
}

filter_markerset <- function(meta_markers,cell_type_name, min_recurrence = 1, min_auroc = 0.5, 
    min_fc = 0, fc_threshold = 3, auroc_threshold = 0.8){
    
    filtered_markers = meta_markers %>% dplyr::select(.data$cell_type, 
        .data$gene, .data$recurrence, .data$auroc, .data$fold_change_detection) %>% 
        dplyr::filter(.data$cell_type == cell_type_name & .data$recurrence >= 
            min_recurrence) %>% dplyr::filter(.data$auroc > min_auroc & 
        .data$fold_change_detection > min_fc) %>% dplyr::mutate(recurrence = factor(.data$recurrence)) %>% 
        dplyr::group_by(.data$cell_type) %>% dplyr::mutate(is_pareto = is_pareto_front(.data$auroc, 
        .data$fold_change_detection))
    return(filtered_markers)
}


is_pareto_front = function(x, y, tolerance = 0) {
    ids = seq_along(x)
    result = data.frame(x, y, ids) %>%
        dplyr::arrange(dplyr::desc(.data$x), dplyr::desc(.data$y)) %>%
        dplyr::mutate(max_y = cummax(.data$y), is_pareto = .data$y >= .data$max_y - tolerance) %>%
        dplyr::select(.data$ids, .data$is_pareto) %>%
        dplyr::arrange(.data$ids) %>%
        dplyr::pull(.data$is_pareto)
    return(result)
}

read_h5ad <- function(h5path, groups = c( 'obs','var','obsm','uns') ,updateProgress =NULL ){

    

    if (FALSE){"
    Currently only reads 1D vectors from obs and var
    obs/var     - 1D vectors for the cell and gene respectively (N / M)

    TODO: include option to read the obsm/varm/obsp/varp/uns, store as lists.
    "}

    # 'var' , 'obs' etc
    output = vector(mode='list' ,length= length(groups) )
    names(output) = groups
    for (vn in groups){
        if (is.function(updateProgress)) {
            text = paste0('Please wait.')
            updateProgress(detail = text)
        }
        flag = 0
        tryCatch({
            outputlist = rhdf5::h5read(h5path, name = vn )
            flag = 1
        },error = function(e){}) 

        if (flag == 1 ){
            if (vn == 'obsm') {
                # ensure that the dimensions match
                # should be  - x nCells
                for (n in names(outputlist)) {
                    if (any(class(outputlist[[n]]) == 'list')) {
                        output[[vn]][[n]] = do.call(cbind.data.frame, outputlist[[n]] )
                    } else if (any(class(outputlist[[n]]) == 'matrix' )) { 
                        output[[vn]][[n]] = data.frame(t(outputlist[[n]]))
                    }
                }

            } else if (endsWith(vn, 'p')){
                #TODO
            } else if(vn =='uns') {
                for (n in names(outputlist)){   
                    # which ones should we convert to dataframes?
                    # if (endsWith(n ,'_df') ){   # should only pick up samp, ext, run dfs
                    #     output[[vn]][[n]] = do.call(cbind.data.frame, outputlist[[n]] )
                    # }else {
                    output[[vn]][[n]] = outputlist[[n]]
                    # }
                }

            } else { # obs or var
                # transform categorical data
                cats = outputlist$'__categories'
                for (cat in names(cats)){
                    ind = outputlist[[cat]]
                    outputlist[[cat]] = outputlist$'__categories'[[cat]][ind+1]
                }

                l = lengths(outputlist)
                output[[vn]] = data.frame(outputlist[l == l['_index'] ])
                # rownames(output[[vn]]) =output[[vn]][['X_index']]
            }
        }
    }

    for (g in groups) {
        tryCatch({
            output[g] = taRifx::remove.factors(output[g])
        },error = function(e){})
    }
        


    rhdf5::h5closeAll()
    return(output)
}

plot_umap_cts <- function(adata, value = NULL ,ranges=NULL){ #'log1p_n_genes_by_counts'
    if (FALSE) {"
    obsname should be a valid column header under obs/ 
    
    TODO specify by <normalization>_<variable>
    "}

    umap = adata$obsm$X_umap
    umap['v'] = value
    colnames(umap) = c('x','y','v')
    
    # umap = data.frame( x = coords[1,] , y = coords[2,] , v = adata$obs[[obsname]])
    
    # TODO how many cells are in the window?? change point size accordingly
    if (is.null(ranges$x ) | is.null(ranges$y )   ){
        N = nrow(umap)
    } else{
        N = sum(umap$x > ranges$x[1] & umap$x < ranges$x[2] & umap$y > ranges$y[1] & umap$y < ranges$y[2] ) 
    }


    # if (N > 100000) {
    #     ptsize = .01
    #     ptshape =16
    # } else{
        ptsize = (100/(N+1))^.5 +.2    # maxes at 10 when N = 1
                            # for N = 100k, s = 0.02 ??
        ptshape =16
    # }

    
    # ptsize = (100/N)^.5     # maxes at 10 when N = 1
                            # for N = 100k, s = 0.25

    
    fig = ggplot2::ggplot(umap, ggplot2::aes(x = x, y = y, col = v )) +
        ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y)+
        # ggplot2::geom_point(shape = 16,  size=ptsize) + #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
        ggplot2::theme_classic() + 
        ggplot2::scale_color_gradient(low = "grey90", high = "dodgerblue4") + #,limits = c(0,NA)
        ggplot2::xlab('UMAP 1') +
        ggplot2::ylab('UMAP 2') + 
        ggplot2::labs(colour = NULL) +
        ggplot2::theme(legend.position='right')
    if (N > 10000)  {
        fig= fig +  ggrastr::geom_point_rast(shape = 16,  size=ptsize,alpha = .75 )
    } else {
        fig = fig +  ggplot2::geom_point(shape = 16,  size=ptsize,alpha = .75)
    }

    return(fig)
}

plot_umap_discrete <- function(adata, value = NULL,  colormap = NULL, ranges=NULL){ #'log1p_n_genes_by_counts'
    if (FALSE) {"
    obsname should be a valid column header under obs/ 
    
    TODO specify by <normalization>_<variable>
    "}

    umap = adata$obsm$X_umap
    umap['v'] = value
    colnames(umap) = c('x','y','v')
    umap$v[umap$v =='NA'] = 'unassigned'
    
    # sort by number of colors - visually aesthetic
    a =names(sort(table(umap$v),decreasing=TRUE))
    umap=umap[order(match( umap$v,a)),]    
    # umap=umap[order(-as.numeric(factor(umap$v))),]


    # umap = data.frame( x = coords[1,] , y = coords[2,] , v = adata$obs[[obsname]])
    if (is.null(colormap)){
        n = length(unique(umap$v))
        colormap = data.frame(label = a, color = DISTINCT_COLORS[1:n] ,stringsAsFactors=FALSE)
    }
    
    umap = merge(umap, colormap ,by.x = 'v',by.y = 'label' ,sort=FALSE)
    # umap$v <- factor(umap$v , levels =a)
    if (is.null(ranges$x ) | is.null(ranges$y )   ){
        N = nrow(umap)
    } else{
        N = sum(umap$x > ranges$x[1] & umap$x < ranges$x[2] & umap$y > ranges$y[1] & umap$y < ranges$y[2] ) 
    }
    

    ptsize = (100/(N+1))^.5 +.2
    ptshape =16
    

    fig = ggplot2::ggplot(umap, ggplot2::aes(x = x, y = y,col = v )) +
        # ggplot2::geom_point(alpha = .75, shape = ptshape, size = ptsize) + #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
        ggplot2::theme_classic() + 
        ggplot2::scale_colour_manual(breaks=umap$v, values=umap$color)+
        # ggplot2::scale_color_gradient(low = "azure1", high = "dodgerblue4") + 
        ggplot2::xlab('UMAP 1') +
        ggplot2::ylab('UMAP 2') + 
        ggplot2::labs(colour = NULL) +
        ggplot2::theme(legend.position='right')+
        ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=1)))

    if (N > 10000)  {
        fig= fig +  ggrastr::geom_point_rast(shape = 16,  size=ptsize,alpha = .75 )
    } else {
        fig = fig +  ggplot2::geom_point(shape = 16,  size=ptsize,alpha = .75)
    }
    return(fig)
}

plot_pca <- function(pca_var_ratio){
    xlab = paste0('PC' , 1:(length(pca_var_ratio)))
    df =  data.frame('PC' = xlab , 'ev'  =pca_var_ratio * 100,stringsAsFactors=FALSE)
    df$PC <- factor(df$PC, levels=xlab)

    fig = ggplot2::ggplot(df, ggplot2::aes(x = PC, y = ev) ) +
        ggplot2::geom_point() + 
        ggplot2::theme_classic() +         
        ggplot2::xlab('PC component') +
        ggplot2::ylab('PCA Explained Variance (%)') 
    return(fig)
}

plot_pr_curve <- function(pr_dfs ,ranges =NULL){
    # df = data.frame(Recall = adata$uns[[id]]['Recall'] , Precision = adata$uns[[id]]['Precision'], Enrichment = adata$uns[[id]]['value'] )
    # pr_df is a list of dataframes.
    for (i in names(pr_dfs)){
        pr_dfs[[i]]['Label_Type'] = i
    }

    
    pr_df = do.call('rbind' , pr_dfs )
    pr_df = filter_PR(pr_df )
    if (is.null(names(pr_dfs) )){
        fig = ggplot2::ggplot(pr_df, ggplot2::aes(x = Recall, y= Precision)) + 
                ggplot2::geom_line()+
                ggplot2::theme_classic() +
                ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    } else{
        fig = ggplot2::ggplot(pr_df, ggplot2::aes(x = Recall, y= Precision,color=Label_Type)) + 
                ggplot2::geom_line()+
                ggplot2::theme_classic() +
                ggplot2::scale_colour_manual(values = DISTINCT_COLORS[1:2])+
                ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y)+
                ggplot2::guides(colour=ggplot2::guide_legend(title='Labelset'))

    }
    return(fig)
}

filter_discrete_labels <- function(adata,discrete_labs){
    # which of the discrete labels are plot-able?
    ncells = nrow(adata$obs)
    counts = lapply(discrete_labs , function(x) {length(unique(adata$obs[,x]))} )
    return(discrete_labs[(counts  <= N_COLORS ) & (counts >1)] )
}

# TODO dataframes aren;t saved to the h5ad files??? Fix - load and replace this function
build_sample_attr_table <- function(obs_df, samp_df ){

    # break down the dictionary format to a sample    
    tmp = paste( gsub("'", '"', samp_df$attribute ) ,collapse = ',')
    sdf = jsonlite::fromJSON(paste('[',tmp,']'))
    rownames(sdf) = samp_df$samp_id
    
    outdf_tmp_cell = vector(mode = 'list', length=length(colnames(sdf)))
    names(outdf_tmp_cell) =colnames(sdf)

    # count how often the pair appear in all samples
    cell_count = merge( obs_df[c('cell_id' ,'samp_id' )] , sdf , by.x = 'samp_id'  ,by.y = 0, all.x = TRUE,all.y=FALSE)

    m = match( obs_df$cell_id , cell_count$cell_id)
    cell_count = cell_count[m, ]
    cell_count[is.na(cell_count) ] = 'NA'
    # cell_count = cell_count[3:ncol(cell_count)]
    for (n in colnames(sdf)){
        tmp_cell = table(cell_count[n])
        outdf_tmp_cell[[n]] = data.frame(Attribute = n, tmp_cell  )
    }

    
    outdfcell = do.call('rbind', outdf_tmp_cell)
    colnames(outdfcell) = c('Attribute','Value','Cell Count')
    outdfcell = outdfcell[outdfcell$'Cell Count' >1,]
    outdfcell=outdfcell[order(outdfcell$'Cell Count',decreasing=TRUE),]

    # which ones can we plot?


    lens = table(outdfcell$Attribute)
    add_to_obs = names(lens[ lens >1 & lens <= N_COLORS ])
    return(list(outdfcell, cell_count[,add_to_obs,drop=FALSE] ))
}

# TODO make prettier 
plot_histogram_bar <- function(adata, val ,geom='histogram',ranges = NULL, xlabel = NULL ,obsname=NULL){
    
    umap = adata$obsm$X_umap
    colnames(umap) = c('x','y')
    umap['v'] = val 
    umap['Group'] = 'Dataset'
    df_filter= NULL
    
    

    if (!is.null(ranges$x) & !is.null(ranges$y) ){
        df_filter = umap[  (umap$x > ranges$x[1]) & (umap$x < ranges$x[2])    , ,drop=FALSE] 
        df_filter = df_filter[  (df_filter$y > ranges$y[1]) & (df_filter$y < ranges$y[2])    , , drop=FALSE] 
        
        
        if (nrow(df_filter) == 0 ) {
            df_filter=NULL
        } else{
            df_filter['Group'] = 'Zoomed'
        }

        
    }

    if (!is.null(df_filter) ){
        umap = rbind(umap, df_filter)
    }
    
    if ( ! is.null(obsname) ){
        if (obsname == 'total_counts' || obsname =='n_genes_by_counts'){
            umap$v = log(umap$v)
            # df_filter$v = log(df_filter$v)

            # xlabel = paste0('log(',xlabel,')')
            
        }
    }
    

    fig = ggplot2::ggplot(umap )+ #  ..count.. )) +
        # ggplot2::geom_histogram( alpha = .5,ggplot2::aes(y = ..density.. ),position = 'identity' ) + #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
        ggplot2::theme_classic() + 
        ggplot2::xlab(xlabel) 

    if (geom =='histogram'){
        
        stat_order = rhdf5::h5read(paste0(DATADIR, '/metadata/proj_bin_heights.hdf5') , name = 'stat_order')
        global_bin_heights = data.frame(t(rhdf5::h5read(paste0(DATADIR, '/metadata/proj_bin_heights.hdf5') , name = 'global_bin_heights' )))
        brks = rhdf5::h5read(paste0(DATADIR, '/metadata/proj_bin_heights.hdf5') , name = 'breaks')
        colnames(global_bin_heights) = stat_order
        rownames(brks) = stat_order
        global_bin_heights = global_bin_heights[,obsname]
        brks = brks[obsname,]


        df = data.frame(global_height =global_bin_heights  , 
            brks = brks[1:(length(brks)-1) ] 
        )



        binsize = df$brks[2] - df$brks[1]
        # df$brks = df$brks + binsize/2
        # global distr
        fig = fig+ 
                ggplot2::geom_col(data = df ,ggplot2::aes(x = brks , y = global_height),  color = NA ,fill = 'lightgrey',alpha =0.7 ,width = binsize, position = ggplot2::position_nudge(binsize/2)) + 
                ggplot2::theme_classic()

        fig = fig + 
            ggplot2::geom_histogram( data=umap,ggplot2::aes(x=v , y = ..density..,fill = Group ), alpha = .6,bins = 50, breaks =brks, color=NA,position = 'identity' )+
            ggplot2::scale_fill_manual(values=DISTINCT_COLORS[1:2]) +
            ggplot2::ylab('Density')
        
 
        # cols = c("LINE2"="#0000cd","BAR"="#8b0000","LINE1"="lightgrey",)

        # fig = fig + 
        #     ggplot2::scale_fill_manual(name="Bar",values=cols) 

        auroc = ALL_AUROCS[adata$uns$exp_df$proj_id[1],obsname ]
        y_pos = max(max(hist(umap$v,breaks = brks,plot=FALSE)$density,na.rm=TRUE) ,max(df$global_height,na.rm=TRUE ),na.rm=TRUE )

        # if (!is.null(ranges$x) & !is.null(ranges$y) ){
        #     maxheight_zoom = max(hist(df_filter$v,breaks = breaks,plot=FALSE)$density)
        #     y_pos = max(y_pos)
        # } 

        
        if (!is.null(auroc)){
            fig = fig + ggplot2::annotate('text',x = quantile(brks,.05) ,y = y_pos  ,label =paste0('AUROC: ', round(auroc,4)))
        }

    } else if(geom =='bar'){
        
        fig = fig + ggplot2::geom_bar(data= umap, ggplot2::aes(x = v, y=(..count..), fill = Group ) , alpha = .5,position = 'identity' ) + #/sum(..count..)
            ggplot2::ylab('Number of Cells')+
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,  vjust = 0.5)) +
            ggplot2::scale_fill_manual(values=c('#0000cd','#8b0000')) 

            # ggplot2::geom_text(aes(label=(..count..)), position=position_dodge(width=0.9), vjust=-0.25)



    }


    return(fig)
}

filter_PR <- function(pr_df){
    # used to make ggplot faster - give it fewer points
    # pr_df = as.data.frame(pr_df)

    # pr_df = data.frame('Precision' = prec, Recall = recall)
    n = nrow(pr_df)
    prec = pr_df$Precision
    recall = pr_df$Recall

    slopes = (prec[3:n] - prec[1:(n-2)]) / (recall[3:n] - recall[1:(n-2)])
    to_keep = slopes != 0 & slopes != -Inf & slopes != Inf

    
    pr_df = pr_df[ c(TRUE, to_keep, TRUE) ,]
    return(pr_df)
}

plot_confusion <- function(val1, val2,val1name,val2name ){
    # for comparing two discrete values
    
    cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, 
        "YlOrRd"))(100)


    df = data.frame(x=val1,y=val2)

    um = unique(df$x)
    m = match(df$x, um)

    un = unique(df$y)
    n = match(df$y, un)
    mn = paste(m,n,sep="|")
    t = table(mn)

    mat = matrix(0, length(um),length(un))
    for (i_j in names(t)){
        ij = as.numeric(strsplit(i_j ,"\\|")[[1]])
        mat[ij[1] , ij[2]] = t[i_j]
    }
    
    rownames(mat) = um
    colnames(mat) = un

    mat = log(mat+1) 
    breaks <- seq(min(mat), max(mat), length = 100)
    ht_cols = circlize::colorRamp2(breaks = brks, colors = cols)
    # png('~/test.png')
    fig = ComplexHeatmap::Heatmap(mat, name = "Log(counts+1)",col=ht_cols ,border=TRUE,
        column_names_rot = 45,row_names_rot = 45,
        row_title=val1name, column_title=val2name,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.1f", mat[i, j] ), x, y, gp = grid::gpar(fontsize = 10))
        }
    )
    # dev.off()

    return(fig)
}

plot_roc <- function(cval, dval, cvalname, dvalname ){
    # cval is continuous
    # dval is discrete
    #TODO filter where slopes are 0 or inf

    df = data.frame(cval__ = cval )
    for (i in unique(dval)){
        df[i] = dval == i
    }
    # df = df[
    
    df = df[order(df$cval__,decreasing=TRUE),]
    
    auc = data.frame(fpr= character(), tpr = character(), group=character())
    for (g in  unique(dval)){
        
        h1 = df[,g]
        h2 = !df[,g]

        auc = rbind( auc, 
            data.frame( FPR= cumsum(h2) / sum(h2),TPR = cumsum(h1) / sum(h1), group = g))
    }
    
    auc$group <- factor(auc$group, levels = colnames(df)[1+order(colSums(df[2:ncol(df)] ),decreasing=TRUE)])
    
    fig = ggplot2::ggplot(auc, ggplot2::aes(x = FPR, y= TPR, group=group,color=group)) +
        ggplot2::geom_line() + ggplot2::theme_classic() +
        ggplot2::guides(colour=ggplot2::guide_legend(title=dvalname))+
        ggplot2::scale_colour_manual( values = DISTINCT_COLORS  )

    return(fig)
    # fig$labels$colour <- "New Legend Title"

}

plot_scatter <- function(val1, val2, val1name,val2name){
    if (startsWith(val1name,'log')) val1 = log(val1)
    if (startsWith(val2name,'log')) val2 = log(val2)
    df = data.frame( x = val1, y = val2 )
    

    N = nrow(df)
    ptsize = (100/(N+1))^.5 +.2


    fig = ggplot2::ggplot(df , ggplot2::aes(x = x , y=y) ) + 
        ggplot2::theme_classic()+
        ggplot2::xlab(val1name) + 
        ggplot2::ylab(val2name)+
        ggplot2::geom_smooth(method = "lm",formula = y ~ x, se = TRUE)+
        ggpubr::stat_cor(method="pearson",na.rm =TRUE, label.x.npc = 'center',label.y.npc = 'bottom')

        if (N > 10000)  {
            fig= fig +  ggrastr::geom_point_rast(shape = 16,  size=ptsize)
        } else {
            fig = fig +  ggplot2::geom_point(shape = 16,  size=ptsize)
        }

    return(fig)
}

graph_it <- function(aurocs,clusts, low_thrs = 0, high_thrs = 1){
    # aurocs
    # thrs = 0.9

    color = c(  # obtained from sampling 25 colors from github:iwanthue using the default preset
        "#aad537","#6a57df","#d2bf00","#3783ff","#0d8900",
        "#fc52c8","#26e190","#fe2982","#01c1ae","#c43c00",
        "#02aaff","#b65d00","#3a47b1","#956000","#c68bff",
        "#5d5700","#931c7d","#00754d","#ff774e","#58d5fb",
        "#962a44","#005f9e","#ffa287","#ebb2fb","#d48da6"
    )

    
    aurocs[is.na(aurocs)] =0
    diag(aurocs) = 0    # removes loops
    studyCell = rownames(aurocs)
    temp = strsplit(studyCell,'\\|')
    study_id= unlist(lapply(temp , function(x) x[1]   ))
    cell_type= unlist(lapply(temp , function(x) x[2]   ))
    
    df = which((aurocs > low_thrs)& (aurocs <high_thrs),arr.ind=TRUE)
    widths = aurocs[df] *2.5
    mass = rowSums(aurocs + t(aurocs)) 
    mass = (mass - min(mass)) / (max(mass) - min(mass)+0.0000001 )# in 10,100
    size = mass*20 + 10 
    
    
    ind = match(study_id,unique(study_id))
    # color = sample(color, length(unique(study_id)))
    nodeColor = color[ind]
    
    edgeColor = rep('gray',length(widths))
    edgeColor[aurocs[df] <=.5] ='#ff8033'

    # clusts = MetaNeighbor::extractMetaClusters(aurocs, threshold=0.7)
    # clustScores = MetaNeighbor::scoreMetaClusters(clusts,aurocs)


    CT2clust  =  rep('',nrow(aurocs))

    names(CT2clust) = studyCell
    for (cls in names(clusts)){
        CT2clust[clusts[[cls]]] = cls
    }



    # arrows =list(to = list(enabled = TRUE, scaleFactor = 2)),
    edges = data.frame(
        to = studyCell[df[,1]], 
        from = studyCell[df[,2]], 
        color = 'grey',
        width = widths,
        smooth=FALSE,
        # arrows = 'from',
        stringsAsFactors=FALSE
    )
    # showModal(modalDialog(paste(edges$color)))
    
    nodes = data.frame(
        id = studyCell, 
        label = cell_type,
        shape='dot',
        group = study_id,
        size = size+10, 
        mass= mass*2 + 1,
        metaCls = CT2clust,
        borderWidth=2, 
        font.size = 18,
        color.border='black',
        color.background = nodeColor,
        stringsAsFactors=FALSE
    )


    legendNodes = data.frame(
        id = unique(ind) , 
        label = unique(study_id),
        color = color[unique(ind)],
        shape='dot'
    )

    
    return(list(nodes,edges,legendNodes))


}

# not used
mn_heatmap = function(aurocs) {
    subclass_to_col = c(
        Lamp5='#DA808C', Sncg='#D633FF', Vip='#B864CC', Sst='#FF9900',
        Pvalb='#D93137', "L2/3 IT"='#C4EC04', "L5 IT"='#50B2AD',
        "L6 IT"='#A19922', "L5 ET"='#0D5B78', "L6 IT Car3"='#5100FF',
        "L5/6 NP"='#3E9E64', "L6 CT"='#2D8CB8', L6b='#53377D', CR='#00FF66',
        Astro='#665C47', OPC='#2E3E39', Oligo='#2E3E39', Endo='#8D6C62',
        SMC='#807059', Peri='#665547', VLMC='#697255', Macrophage='#94AF97'
    )
    dataset_to_col = c(scSS='#EC2028', scCv3='#F15A25',
                       scCv2='#F68921',
                       snSS='#7E9C3D', snCv3M='#B5D554',
                       snCv3Z='#D3E3A6', snCv2='#F2F5C3')
    cell_type = MetaNeighbor::getCellType(colnames(aurocs))
    subclass_label = cell_type
    for (s in names(subclass_to_col)) {subclass_label[startsWith(cell_type, s)] = s}
    row_order <- order(factor(subclass_label, levels = names(subclass_to_col)),
                      na.last = TRUE)
    aurocs = aurocs[row_order, row_order]
    subclass_label = subclass_label[row_order]

    gplots::heatmap.2(
        x = aurocs, margins = c(11,11), key = TRUE, keysize = 1, key.xlab="AUROC",
        key.title="NULL", offsetRow=0.1, offsetCol=0.1, trace = "none",
        Rowv = FALSE, Colv = FALSE, labRow = NA, labCol = NA, na.color = gray(0.95), 
        dendrogram = "none",  density.info = "none",
        col = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20)),
        breaks = seq(0, 1, length = 21),
        RowSideColors = subclass_to_col[subclass_label],
        ColSideColors = dataset_to_col[MetaNeighbor::getStudyId(colnames(aurocs))]
    )

    par(lend = 1)
    legend("topleft", inset = c(0.82, 0.28), legend = names(subclass_to_col),
           col = subclass_to_col, cex = 0.7, lwd = 10, bty="n")
    legend("topleft", inset = c(0.82, 0.03), legend = names(dataset_to_col),
           col = dataset_to_col, pt.cex = 1, cex = 0.7, lwd = 10, bty="n")
}

create_HTML_link <- function(main_url, id, displayname ) {
    
    return( paste0("<a href='",main_url,"%22", id,"%22",  "' target='_blank'>",displayname ,"</a>"))
    # paste0("<a href='https://gillisweb.cshl.edu/FAIRit/?_inputs_&projid=%22" ,projid,"%22> target='_blank' ",projid,'</a>' ) 
}

create_HTML_img <- function(  id ){
    
    sprintf('<img src="%s_umap.png"  width="80" height="80">', id)
    # paste0('<img src="' ,id,  '_umap.png" height="52">')
}

################################################################################
#################################### Driver ####################################
################################################################################


function(input,output,session){

### TAB ITEM
output$bookmarkurl = renderUI({
    bookmarkButton()    
})



{ # search datasets 

    search_df_out = reactive({
        tryCatch({
        outdf = SEARCH_DF
        rownames(outdf) = create_HTML_link('https://gillisweb.cshl.edu/FAIRit/?_inputs_&tabs=%22project%22&projid=', outdf$proj_id,create_HTML_img(outdf$proj_id))
        outdf['pid'] = outdf$proj_id
        outdf$proj_id = create_HTML_link('https://gillisweb.cshl.edu/FAIRit/?_inputs_&tabs=%22project%22&projid=', outdf$proj_id,outdf$proj_id) 
        


        # filter by tech
        if (input$is10x ){
            is10x = outdf$chromium > 0 
        } else{
            is10x = outdf$chromium < 0 
        }
        if (input$isSS ){
            isSS = outdf$smartseq > 0 
        } else{
            isSS = outdf$smartseq < 0 
        }
        # filter by datasource
        if (input$isNemo ){
            isNemo = outdf$datasource =='NeMO'
            
        } else{
            isNemo = outdf$datasource == '' 
        }
        if (input$isSRA ){
            isSRA = outdf$datasource =='SRA'
        } else{
            isSRA = outdf$datasource ==''
        }

        # filter by class/subclass
        class_flag = outdf$primary_class %in% input$primary_class 
        subclass_flag = outdf$primary_subclass %in% input$primary_subclass 
        
        # filter by date 
        message(zoo::as.Date(zoo::as.yearmon(input$date_slider[1])  ))
        message(zoo::as.Date(zoo::as.yearmon(input$date_slider[1]) ,frac=1  ))

        message(class(input$date_slider[1]))
        message(as.Date(input$date_slider[1] ,format = '%b %Y' ))
        
        date_flag = SEARCH_DF$pdate >= zoo::as.Date(zoo::as.yearmon(input$date_slider[1]),frac=-1  )& 
            SEARCH_DF$pdate <=  zoo::as.Date(zoo::as.yearmon(input$date_slider[2]) ,frac=1  )

        # input$date_slider[2]

        # filter by cts values 
        # ctsvals = c('gini_ds','ncells','percent_dropout','mt','class_pr','subclass_pr')

        # gini
        gini_flag = (outdf$gini_ds >= input$gini_ds_slider[1]) &(outdf$gini_ds <= input$gini_ds_slider[2])
        ncell_flag = (outdf$ncells >= as.numeric(gsub(',','',input$ncell_slider[1]))) &
            (outdf$ncells <= as.numeric(gsub(',','',input$ncell_slider[2])))
        percent_dropout_flag = (outdf$percent_dropout >= input$percent_dropout_slider[1]) &(outdf$percent_dropout <= input$percent_dropout_slider[2])
        mt_flag = (outdf$mt >= input$mt_slider[1]) &(outdf$mt <= input$mt_slider[2])
        class_pr_flag = (outdf$class_pr >= input$class_pr_slider[1]) &(outdf$class_pr <= input$class_pr_slider[2])
        subclass_pr_flag = (outdf$subclass_pr >= input$subclass_pr_slider[1]) &(outdf$subclass_pr <= input$subclass_pr_slider[2])
    
        tech_flags = is10x | isSS
        datasource_flags = isNemo | isSRA
        label_flag = subclass_flag & class_flag

        all_flags = datasource_flags & tech_flags & label_flag & gini_flag & date_flag &
                ncell_flag&percent_dropout_flag&mt_flag &class_pr_flag&subclass_pr_flag

        outdf = outdf[all_flags,,drop=FALSE]
        outdf$ncells = format((outdf$ncells), big.mark=",")
        outdf$chromium = format((outdf$chromium), big.mark=",")
        outdf$smartseq = format((outdf$smartseq), big.mark=",")
        colnames(outdf) = c('Project','Title','Source','Date',
            'Class','Subclass', ' #Cells','#10x','#SS',
            'Class PR','Subclass PR', 'Gini', '%Dropout', 'Mito','pid'
        )
        outdf[order(outdf$Date,decreasing=TRUE),]
        # outdf
        },error = function(e){})
    })

    output$dataset_search <- DT::renderDataTable({
        # tryCatch({
            search_df_out()[,1:(ncol(search_df_out() )-1) ]
        # },error = function(e){})
    },
    
    rownames = TRUE,
    escape = FALSE,
    server=FALSE,
    
    options = list( 
        rowCallback = DT::JS(
            "function(row, data) {",
            "var full_text = data[2]",
            "$('td', row).attr('title', full_text);",
            "}"),
        columnDefs = list(list( 
            targets = 15,   # drop column 15 -> abstract bu still make it searchable.
            searchable = TRUE,
            visible = FALSE
        )
        )
    ),


    ) #


    # download dataset
    
    output$downloadoptions <- renderUI({
        prettyRadioButtons(
            inputId = "download_option",
            label = paste0("Download filtered metadata (",length(unique(search_df_out()$pid)) ,' projects):'),
            thick = TRUE,
            choiceNames = c('Projects', 'Samples','Experiments','Runs'),
            choiceValues = c('projects', 'samples','experiments','runs'),
            animation = "pulse",
            status = "info",outline = T
        )
    })


    # button to download the top connections of a gene
    output$downloadButton <- renderUI({
        tryCatch({
            downloadButton('downloader',
            label='Download')#color='primary',
            # style='stretch') #size='m',block=FALSE)  
        }, error = function(e){})
    })

    # download executor
    output$downloader <- downloadHandler(
        filename = function() {
            paste0(input$download_option ,'_metadata.tsv')
        },
        content = function(con) {
            if(length(input$download_option)==0){}
            else{switch(input$download_option,
                'runs' = {
                    rdf = read.delim(paste0(DATADIR, '/metadata/runs_.tsv'),sep="\t")
                    rdf = rdf[rdf$proj_id %in% search_df_out()$pid,,drop=FALSE]
                    write.table(rdf, con , row.names = FALSE,sep='\t',quote =FALSE) 
                },
                'samples' = {
                    sdf = read.delim(paste0(DATADIR, '/metadata/samples_.tsv'),sep="\t")
                    sdf = sdf[sdf$proj_id %in% search_df_out()$pid,,drop=FALSE]
                    write.table(sdf, con , row.names = FALSE,sep='\t',quote =FALSE) 
                },
                'experiments' = {
                    edf = read.delim(paste0(DATADIR, '/metadata/experiments_.tsv'),sep="\t")
                    edf = edf[edf$proj_id %in% search_df_out()$pid,,drop=FALSE]
                    write.table(edf, con , row.names = FALSE,sep='\t',quote =FALSE) 
                },
                'projects' = {
                    pdf = read.delim(paste0(DATADIR, '/metadata/projects_.tsv'),sep="\t")
                    pdf = pdf[pdf$proj_id %in% search_df_out()$pid,,drop=FALSE]
                    write.table(pdf, con , row.names = FALSE,sep='\t',quote =FALSE) 
                }
            )}
        } 
    )


}


{   # START MAIN PROJECT TAB
    # render bookmark button for testing purposes. Bookmark used in product chart 
    ## TODO : REMOVE BOOKMARK BUTTON IN FINAL VERSION


    #### DEFINE REACTIVES ####
    # observeEvent(input$projid ,{
    #     discrete_stat_options 
    # })

    h5path <- reactive( # path to the h5ad file.
        paste0(DATADIR, '/',input$projid, '.h5ad')
    )


    # takes at least 25 seconds for Arlotta dataset
    adata <- reactive({# read the h5ad and store as adata
        tryCatch({
            if (input$tabs=='project' ){

                progress <- shiny::Progress$new(style = 'notification')
                progress$set(message = "Reading in the dataset", value = 0)
                # Close the progress when this reactive exits (even if there's an error)
                on.exit(progress$close())

                updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                    value <- progress$getValue()
                    value <- value + (progress$getMax() - value) / 5
                }
                progress$set(value = value, detail = detail)
                }

                read_h5ad(
                    h5path(),
                    groups = c('obs','obsm','uns'),
                    updateProgress
                )
            }
        },error=function(e){})
    })

    ext_ids <- reactive({ # external IDs (GEO/)
        tryCatch({
        tmp = gsub("\\{","", adata()$uns$header_box$ext_ids )
        tmp = gsub("\\}","", tmp)
        tmp = gsub("\\'","", tmp)
        tmp = unlist(strsplit(tmp,', '))
        strsplit(tmp," ")
        # vals = paste(lapply(ext_ids,function(x) x[2]),collapse="\n")
        },error=function(e){})

    })


    #### TEXT/TABLE OUTPUTS ####
    # get title/abstract text
    output$title <- renderText({ 
        tryCatch({
            adata()$uns$header_box$title 
        },error=function(e){})

    })

    output$abstract <- renderText({
        tryCatch({
            adata()$uns$header_box$abstract
        },error=function(e){})

    })

    output$labshare_url <- renderUI({
        tryCatch({

            labshareurl = paste0('https://labshare.cshl.edu/shares/gillislab/resource/FAIRit/',input$projid)
            
            HTML('Download parsed data for',input$projid , paste0("<a href='",labshareurl, "' target='_blank'>here</a>.") )
        },error=function(e){})

    })

    output$ext_ids_sub <- renderUI({
        tryCatch({
            subs = unlist(lapply(ext_ids(), '[[',1))
            HTML(paste(subs,collapse="<br/><br/>") )
        },error=function(e){})

    })

    output$ext_ids <- renderUI({
        tryCatch({
            ids = unlist(lapply(ext_ids(), '[[',2))

            prj_ids = grepl('PRJ',ids)
            geo_ids =grepl('GSE',ids)

            # include urls if prjna or geo. More urls eventually.
            ids[prj_ids]  = paste0("<a href='",BIOPROJ_URL,ids[prj_ids],  "' target='_blank'>",ids[prj_ids] ,"</a>")
            ids[geo_ids]  = paste0("<a href='",GEO_URL,    ids[geo_ids],  "' target='_blank'>",ids[geo_ids] ,"</a>")

            HTML(paste(ids,collapse="<br/><br/>") )
        },error=function(e){})

    })

    # get basic stats
    output$basic_info <- renderTable({
        tryCatch({

            pdat = as.Date(adata()$uns$header_box$publish_date)
            pdat = format(pdat, "%d %b %Y" )

            tech_counts = adata()$uns$header_box$tech_count
            tx_counts = sum(as.numeric(tech_counts[startsWith( names(tech_counts) ,'10x' )]))
            ss_counts = sum(as.numeric(tech_counts[startsWith( names(tech_counts) ,'smartseq' )]))


            data.frame(col1 = c('Publish Date:','Sample count:','Experiment count:','Run count:') ,
                col2 = c(as.character(pdat),        
                    format(as.numeric(adata()$uns$header_box$n_samp), big.mark=","),
                    format(as.numeric(adata()$uns$header_box$n_exps), big.mark=","),
                    format(as.numeric(adata()$uns$header_box$n_runs), big.mark=",")
                ),
                col3 = c('','Cell count:','SmartSeq cells:','10x cells:'),
                col4 = c('',
                    format(as.numeric(adata()$uns$header_box$n_cell) , big.mark=","),
                    format(as.numeric(ss_counts), big.mark=","),
                    format(as.numeric(tx_counts), big.mark=",")
                )
            )

        },error=function(e){})

        
    },colnames=FALSE)

    # parse sample attributes
    # TODO get number of cells per attribute/value 
    samp_attr <- reactive({
        tryCatch({
        build_sample_attr_table(adata()$obs , adata()$uns$samp_df )
        # add_to_obs = samp_attr[[2]]
        # samp_df = samp_attr[[1]]
        },error=function(e){})

    })


    # output$downloadButton_ds <- renderUI({
    #     tryCatch({
    #         downloadButton('downloader_ds',
    #         label='Download')#color='primary',
    #         # style='stretch') #size='m',block=FALSE)  
    #     }, error = function(e){})
    # })


    # output$downloader_ds <- downloadHandler(
    #     filename = function() {
    #         paste0(input$projid ,'_data.h5ad')
    #     },
    #     content = function(con) {
    #         file.copy('ERP122596.h5ad')
    #      } #,contentType ='hdf5'
    # )

    # uiOutput('downloadButton_ds')  

    output$sample_attr <- DT::renderDataTable({
        tryCatch({
            sdf = samp_attr()[[1]]
            sdf['Cell Count'] = format(as.numeric(sdf[,'Cell Count']), big.mark=",")
            sdf[,c('Attribute','Value','Cell Count')]
        },error=function(e){})

        },
        rownames=FALSE, server=FALSE,
        options = list(autoWidth = FALSE, pageLength = 5)

    )

    # Get available options to plot.
    discrete_stat_options <- reactive({
        tryCatch({

        dl = filter_discrete_labels(adata() , DISCRETE_LABS)
        sa = colnames(samp_attr()[[2]])

        names(sa) = sa

        dso = c( dl, sa[sa!='samp_id'] )
        dso[! duplicated(dso)]
        },error=function(e){})

    })

    observeEvent( discrete_stat_options(),{
        updatePickerInput( session , 'obsname', choices = list(' '=CTS_LABS[1:4] , 'Gene set AUROCs' =   CTS_LABS[5:length(CTS_LABS)], 'Author annotations' = discrete_stat_options() ),selected= CTS_LABS[1] )
    })

    observeEvent( discrete_stat_options(),{
        updatePickerInput( session , 'obsname2', choices = list(' '=CTS_LABS[1:4] , 'Gene set AUROCs' =   CTS_LABS[5:length(CTS_LABS)], 'Author annotations' = discrete_stat_options() ),selected= CTS_LABS[2] )
    })
    
    
    # adata()$obs = merge(adata()$obs,  samp_attr()[[2]], by = 'cell_id' )

    #### FIGURE OUTPUTS ####
    
    # brush window - zoomed in umap
    # ranges1 <- reactiveValues(x = NULL, y = NULL)
    ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$umap_dblclick, {
        brush <- input$umap_brush
        if (!is.null(brush)) {
            ranges$x <- c(brush$xmin, brush$xmax)
            ranges$y <- c(brush$ymin, brush$ymax)

        } else {
            ranges$x <- NULL
            ranges$y <- NULL
        }

    })

    # vals <- reactiveValues(
    #     keeprows = rep(TRUE, nrow(adata()$obs ))
    # )   

    # observeEvent(input$exclude_toggle, {
    #     res <- brushedPoints(adata()$obs, input$plot1_brush, allRows = TRUE)

    #     vals$keeprows <- xor(vals$keeprows, res$selected_)
    # })

    output$umap <- renderPlot({
        

        tryCatch({
            progress <- shiny::Progress$new(style = 'notification')
            progress$set(message = "Drawing the UMAP", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                    value <- progress$getValue()
                    value <- value + (progress$getMax() - value) / 5
                }
                progress$set(value = value, detail = detail)
            }

            if (input$obsname %in% discrete_stat_options()) {
                obsname = input$obsname
                if (obsname %in% colnames(samp_attr()[[2]])){
                    val = as.character(samp_attr()[[2]][,obsname])
                } else{
                    val = as.character(adata()$obs[,obsname])
                }
                fig = plot_umap_discrete(adata(), value = val , colormap= NULL, ranges )
            } else if ( input$obsname %in% colnames(adata()$obs) ) {
                obsname = input$obsname

                vals = adata()$obs[obsname][,1]
                if (obsname == 'total_counts' | obsname == 'n_genes_by_counts'){
                    vals = log(vals)
                }
                    
                fig = plot_umap_cts(adata(), value =vals ,ranges )
                # message('cts var')
            } else {
                # message(input$obsname)
                fig = NULL
            }

            
            updateProgress(detail = 'UMAP drawn')              
            fig
        } , error = function(e){e})
    })

    output$hist <- renderPlot({
        tryCatch({
            
            progress <- shiny::Progress$new(style = 'notification')
            progress$set(message = "Drawing the histogram", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                    value <- progress$getValue()
                    value <- value + (progress$getMax() - value) / 5
                }
                progress$set(value = value, detail = detail)
            }

            if (input$obsname %in% discrete_stat_options()) {
                obsname = input$obsname
                disname = names(discrete_stat_options()[discrete_stat_options() == obsname])
                if (obsname %in% colnames(samp_attr()[[2]])){
                    val = as.character(samp_attr()[[2]][,obsname])
                } else{
                    val = as.character(adata()$obs[,obsname])
                }
                fig = plot_histogram_bar(adata(), val ,geom='bar',ranges = ranges, xlabel = disname , obsname=NULL)

                # fig = plot_umap_discrete(adata(), value = val , colormap= NULL, ranges )
            } else if ( input$obsname %in% colnames(adata()$obs) ) {
                obsname = input$obsname

                disname = names(CTS_LABS[CTS_LABS == obsname])
                # fig = plot_umap_cts(adata(), value = ,ranges )
                val=adata()$obs[,obsname]
                fig = plot_histogram_bar(adata(),val ,geom='histogram',ranges = ranges, xlabel = disname,obsname=obsname)
                
            } else {
                obsname = NULL
            }

            # TODO dropdown for varnames plus normalization
            # adata, val ,geom='histogram',ranges = NULL, xlabel = NULL )

            
            # fig = plot_histogram(adata(), obsname = obsname,ranges = ranges)
            updateProgress(detail = 'Histogram drawn')              

            fig
        } , error = function(e){e})
    })


    # TEXT BOX

    output$stat1_info <- renderUI({
        tryCatch({
            if (input$obsname %in% CTS_LABS){
                auc = ALL_AUROCS[input$projid, input$obsname ]
                cts_stat = names(which(CTS_LABS==input$obsname))
                text = list('Comparing the distribution of <b>',cts_stat,'</b> (in blue) 
                    in this project to the global distribution (in lightgray), 
                    we see an AUROC of <b>', auc,'</b>, ' )
                if (input$obsname == 'total_counts'){
                    if (auc > .6 ){
                        tech = '<b>SmartSeq</b>'
                    } else{
                        tech = '<b>Chromium 10x</b>'
                    }

                    text = append(text ,paste0(' which is within expectations for ', tech, ' experiments.') )
                } else if (input$obsname %in% c('mt','ribo','essential','housekeeping','highly_variable','top_100') ) {
                    if (auc > .8 ){
                        adj = 'significantly higher than'                        
                    } else if (auc > .6) {
                        adj = 'slightly higher than'                        
                    } else if (auc > .4) {
                        adj = 'comparable to'                        
                    } else if (auc > .2) {
                        adj = 'slightly lower than'                                               
                    } else if (auc > 0) {
                        adj = 'significantly lower than'                        
                    } 
                    text = append(text ,paste0(' which is ',adj,' the global distribution for ', cts_stat,'.') )
                # } else if () {
                    
                }
                
                HTML(paste0(text,collapse = '') )
            }
        },error = function(e){})
    })

    output$pca_expvar <- renderPlot({
        tryCatch({
            plot_pca( adata()$uns$pca$variance_ratio )
        },error = function(e){})
    })

    output$hist2 <- renderPlot({
        
        tryCatch({

            progress <- shiny::Progress$new(style = 'notification')
            progress$set(message = "Drawing the histogram", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                    value <- progress$getValue()
                    value <- value + (progress$getMax() - value) / 5
                }
                progress$set(value = value, detail = detail)
            }

            if (input$obsname2 %in% discrete_stat_options()) {
                obsname = input$obsname2
                disname = names(discrete_stat_options()[discrete_stat_options() == obsname])
                
                if (obsname %in% colnames(samp_attr()[[2]])){
                    val = as.character(samp_attr()[[2]][,obsname])
                } else{
                    val = as.character(adata()$obs[,obsname])
                }
                fig = plot_histogram_bar(adata(), val ,geom='bar',ranges = ranges, xlabel = disname,obsname=NULL)

                # fig = plot_umap_discrete(adata(), value = val , colormap= NULL, ranges )
            } else if ( input$obsname2 %in% colnames(adata()$obs) ) {
                obsname = input$obsname2
                disname = names(CTS_LABS[CTS_LABS == obsname])
                # fig = plot_umap_cts(adata(), value = ,ranges )
                val=adata()$obs[,obsname]
                fig = plot_histogram_bar(adata(),val ,geom='histogram',ranges = ranges, xlabel = disname,obsname=obsname)

            } else {
                obsname = NULL
            }

            # TODO dropdown for varnames plus normalization
            # adata, val ,geom='histogram',ranges = NULL, xlabel = NULL )

            
            # fig = plot_histogram(adata(), obsname = obsname,ranges = ranges)
            updateProgress(detail = 'Histogram drawn')              
            fig
        } , error = function(e){e})
    })

    output$stat_comp <- renderPlot({
        tryCatch({
        
            progress <- shiny::Progress$new(style = 'notification')
            progress$set(message = "Drawing stat comparison", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            updateProgress <- function(value = NULL, detail = NULL) {
                if (is.null(value)) {
                    value <- progress$getValue()
                    value <- value + (progress$getMax() - value) / 5
                }
                progress$set(value = value, detail = detail)

            }

            
            if (input$obsname %in% discrete_stat_options()) {
                obsname = input$obsname
                disname1 = names(discrete_stat_options()[discrete_stat_options() == obsname])

                type1 = 'discrete'
                
                if (obsname %in% colnames(samp_attr()[[2]])){
                    val1 = as.character(samp_attr()[[2]][,obsname])
                } else{
                    val1 = as.character(adata()$obs[,obsname])
                }

            # is a cts variable
            } else if ( input$obsname %in% colnames(adata()$obs) ) {
                type1 = 'cts'
                obsname = input$obsname
 
                disname1 = names(CTS_LABS[CTS_LABS == obsname])
                # fig = plot_umap_cts(adata(), value = ,ranges )
                val1=adata()$obs[,obsname]
            }

            
            if (input$obsname2 %in% discrete_stat_options()) {
                obsname2 = input$obsname2
                disname2 = names(discrete_stat_options()[discrete_stat_options() == obsname2])
                                
                type2 = 'discrete'

                if (obsname2 %in% colnames(samp_attr()[[2]])){
                    val2 = as.character(samp_attr()[[2]][,obsname2])
                } else{
                    val2 = as.character(adata()$obs[,obsname2])
                }

            # is a cts variable
            } else if ( input$obsname2 %in% colnames(adata()$obs) ) {
                type2= 'cts'
                obsname2 = input$obsname2
                disname2 = names(CTS_LABS[CTS_LABS == obsname2])
                # fig = plot_umap_cts(adata(), value = ,ranges )
                val2=adata()$obs[,obsname2]
            }

            
            if (type1 == 'discrete' & type2 =='discrete'){
                fig = plot_confusion(val1,val2, disname1, disname2)
            } else if (type1 == 'cts' & type2 =='cts'){
                fig = plot_scatter(val1,val2, disname1, disname2)
            } else if (type1 =='cts' & type2 =='discrete' ){
                fig = plot_roc(val1, val2, disname1, disname2 )
            } else if (type1 =='discrete' & type2 =='cts' ){
                fig = plot_roc(val2, val1, disname2, disname1 )
            }

            updateProgress(detail = 'Stat comparision drawn')
            fig

        }, error = function(e){})
    })
} # end tab 1 section


{   # START METAMARKERS section
    # update picker for labels
    # updatePickerInput( session , 'obsname', choices = c('log1p_n_genes_by_counts','b','c') ,selected= 'log1p_n_genes_by_counts' )

    class_choices <- reactive({
        # if (input$labelset == 'Class'){        
        choices = c('All',gsub('Other', 'Non-neuronal',unique(adata()$uns$MetaMarkers$class_assign$predicted))) 
        # } else if (input$labelset == 'Subclass'){
        #     choices = c('Predicted',unique(adata()$obsm$subclass_pred$predicted) )
        # } else{
        #     choices = NULL
        # }
        choices = gsub("/", "-",choices)
        choices[!grepl('NA|unassigned',choices) ]
    })

    subclass_choices <- reactive({

        choices = c('All',unique(adata()$uns$MetaMarkers$subclass_assign_hier$predicted) )

        choices = gsub("/", "-",choices)
        choices[!grepl('NA|unassigned',choices) ]
    })


    observeEvent(class_choices(),{
        output$class_label_box <- renderUI({
            pickerInput(
                inputId="class_label",
                label="Class:", multiple=FALSE,
                choices=class_choices(),
                selected = class_choices()[1],
                width='fit',
                options = pickerOptions(
                    liveSearch=FALSE
                ),
                inline=TRUE
            )
        })                    
    })


    observeEvent(subclass_choices(),{
        output$subclass_label_box <- renderUI({
            pickerInput(
                inputId="subclass_label",
                label="Subclass:", multiple=FALSE,
                choices=subclass_choices(),
                selected = subclass_choices()[1],
                width='fit',
                options = pickerOptions(
                    liveSearch=FALSE
                ),
                inline=TRUE
            )
        })                    
    })



    observeEvent(input$class_label,{
        if (input$class_label != 'All' ){
            output$class_mm_metric_box <- renderUI({
                pickerInput(
                    inputId="class_mm_metric",
                    label="Metric:", multiple=FALSE,
                    choices=c('Predicted','Enrichment','Scores'),
                    selected = 'Predicted' ,
                    width='fit',
                    options = pickerOptions(
                        liveSearch=FALSE
                    ),
                    inline=TRUE
                )
            })       
        } else{
            output$class_mm_metric_box <- renderUI({
            })   
        }
    })


    class_ranges_mm <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$class_umap_mm_dblclick, {
        brush <- input$class_umap_mm_brush
        if (!is.null(brush)) {
            class_ranges_mm$x <- c(brush$xmin, brush$xmax)
            class_ranges_mm$y <- c(brush$ymin, brush$ymax)

        } else {
            class_ranges_mm$x <- NULL
            class_ranges_mm$y <- NULL
        }
    })

    output$class_umap_mm <- renderPlot({
        tryCatch({
            # if (input$labelset =='Class' ){

            progress <- shiny::Progress$new(style = 'notification')
            progress$set(message = "Drawing the UMAP", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            updateProgress <- function(value = NULL, detail = NULL) {
            if (is.null(value)) {
                value <- progress$getValue()
                value <- value + (progress$getMax() - value) / 5
            }
            progress$set(value = value, detail = detail)
            }

            labs = 'class'
            colormap = CLASS_COLORS


            if (input$class_label == 'All' ){
                met= '_pred'
                vals = gsub('Other', 'Non-neuronal', adata()$uns$MetaMarkers$class_assign$predicted)
                fig = plot_umap_discrete(
                    adata(), 
                    value = vals, 
                    colormap=colormap,
                    range  = class_ranges_mm
                ) 
            } else if ( input$class_mm_metric == 'Predicted'){
                tmp = adata()$uns$MetaMarkers$class_assign
                tmp$predicted= gsub("/","-",tmp$predicted)
                tmp$predicted = gsub('Other', 'Non-neuronal',tmp$predicted)
                values = tmp$predicted == input$class_label
                values[values == FALSE] = 'unassigned'
                values[values == TRUE] = input$class_label

                fig = plot_umap_discrete(
                    adata(), 
                    value = values,
                    colormap=colormap,
                    range  = class_ranges_mm) 
                    
                                    
            } else if (input$class_mm_metric == 'Enrichment') {
                tmp = adata()$uns$MetaMarkers$class_enr
                col = grep(paste0('\\|',input$class_label,'$') ,names(tmp) )
                fig = plot_umap_cts(
                    adata(), 
                    value = as.numeric(tmp[[col]]), # find the cell type at the end
                    range = class_ranges_mm) 
            
                
                # vals = adata()$uns$MetaMarkers$class_enr

            } else if (input$class_mm_metric == 'Scores') {
                tmp = adata()$uns$MetaMarkers$class_score
                col = grep(paste0('\\|',input$class_label,'$') ,names(tmp) )
                fig = plot_umap_cts(
                    adata(), 
                    value = as.numeric(tmp[[col]]), # find the cell type at the end
                    range = class_ranges_mm) 
            
                
                # vals = adata()$uns$MetaMarkers$class_score
            } 

            updateProgress(detail='UMAP drawn.')
            fig
        },error=function(e){})
    })

    # render PR
    class_ranges_mm_pr <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$class_umap_mm_pr_dblclick, {
        brush <- input$class_umap_mm_pr_brush
        if (!is.null(brush)) {
            class_ranges_mm_pr$x <- c(brush$xmin, brush$xmax)
            class_ranges_mm_pr$y <- c(brush$ymin, brush$ymax)

        } else {
            class_ranges_mm_pr$x <- NULL
            class_ranges_mm_pr$y <- NULL
        }
    })

    output$class_mm_pr <- renderPlot({
        tryCatch({

            pr_dfs = list(as.data.frame(adata()$uns$MetaMarker$class_PR),as.data.frame(adata()$uns$MetaMarker$subclass_PR_hier) )
            names(pr_dfs) = c('Class' ,'Subclass')
            
            fig = plot_pr_curve(
                pr_dfs, 
                ranges = class_ranges_mm_pr
                ) 


            fig
        },error=function(e){})
    })



    observeEvent(input$subclass_label,{
        if (input$subclass_label != 'All' ){
            output$subclass_mm_metric_box <- renderUI({
                pickerInput(
                    inputId="subclass_mm_metric",
                    label="Metric:", multiple=FALSE,
                    choices=c('Predicted','Enrichment','Scores'),
                    selected = 'Predicted' ,
                    width='fit',
                    options = pickerOptions(
                        liveSearch=FALSE
                    ),
                    inline=TRUE
                )
            })       
        } else{
            output$subclass_mm_metric_box <- renderUI({
            })   
        }
    })


    subclass_ranges_mm <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$subclass_umap_mm_dblclick, {
        brush <- input$subclass_umap_mm_brush
        if (!is.null(brush)) {
            subclass_ranges_mm$x <- c(brush$xmin, brush$xmax)
            subclass_ranges_mm$y <- c(brush$ymin, brush$ymax)

        } else {
            subclass_ranges_mm$x <- NULL
            subclass_ranges_mm$y <- NULL
        }
    })

    output$subclass_umap_mm <- renderPlot({
        tryCatch({
        
            progress <- shiny::Progress$new(style = 'notification')
            progress$set(message = "Drawing the UMAP", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            updateProgress <- function(value = NULL, detail = NULL) {
            if (is.null(value)) {
                value <- progress$getValue()
                value <- value + (progress$getMax() - value) / 5
            }
            progress$set(value = value, detail = detail)
            }
            # if (input$labelset =='Class' ){            labs = 'class'
            colormap = SUBCLASS_COLORS


            if (input$subclass_label == 'All' ){
                
                vals = adata()$uns$MetaMarkers$subclass_assign_hier$predicted
                fig = plot_umap_discrete(
                    adata(), 
                    value = vals, 
                    colormap=colormap,
                    range  = subclass_ranges_mm
                ) 
            } else if ( input$subclass_mm_metric == 'Predicted'){
                tmp = adata()$uns$MetaMarkers$subclass_assign_hier
                tmp$predicted=gsub("/","-",tmp$predicted)

                values = tmp$predicted == input$subclass_label
                values[values == FALSE] = 'unassigned'
                values[values == TRUE] = input$subclass_label

                fig = plot_umap_discrete(
                    adata(), 
                    value = values,
                    colormap=colormap,
                    range  = subclass_ranges_mm) 
                    
                                    
            } else if (input$subclass_mm_metric == 'Enrichment') {
                tmp = adata()$uns$MetaMarkers$subclass_enr_hier
                col = grep(paste0('\\|',input$subclass_label,'$') ,names(tmp) )
                fig = plot_umap_cts(
                    adata(), 
                    value = as.numeric(tmp[[col]]), # find the cell type at the end
                    range = subclass_ranges_mm) 
            
            } else if (input$subclass_mm_metric == 'Scores') {
                tmp = adata()$uns$MetaMarkers$subclass_scores_hier
                col = grep(paste0('\\|',input$subclass_label,'$') ,names(tmp) )
                fig = plot_umap_cts(
                    adata(), 
                    value = as.numeric(tmp[[col]]), # find the cell type at the end
                    range = subclass_ranges_mm) 
            } 
            updateProgress(detail ='UMAP drawn')
            fig
        },error=function(e){})
    })

    # render PR
    subclass_ranges_mm_pr <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$subclass_umap_mm_pr_dblclick, {
        brush <- input$subclass_umap_mm_pr_brush
        if (!is.null(brush)) {
            subclass_ranges_mm_pr$x <- c(brush$xmin, brush$xmax)
            subclass_ranges_mm_pr$y <- c(brush$ymin, brush$ymax)

        } else {
            subclass_ranges_mm_pr$x <- NULL
            subclass_ranges_mm_pr$y <- NULL
        }
    })

    # output$class_mm_pr <- renderPlot({
    #     # tryCatch({

    #         id1 = 'MetaMarker_subclass_PR'
    #         id2 = 'MetaMarker_subclass_PR_heir'

    #         pr_dfs = list(as.data.frame(adata()$uns$MetaMarker$subclass_PR),
    #                       as.data.frame(adata()$uns$MetaMarker$subclass_PR_hier))
    #         names(pr_dfs) = c('classical','hierarchical')
    #         fig = plot_pr_curve(
    #             pr_dfs, 
    #             ranges = subclass_ranges_mm_pr
    #             ) 


    #         fig
    #     # },error=function(e){})
    # })



    # output$enrich_slider <-renderUI({
    #     tryCatch({

    #     })
    # })
}


{   # START ABOUT SECTION

    # pickerInput
    
    # read aurocs 
    aurocs1v1_mat <-reactive({
        a = as.matrix(read.delim(paste0(DATADIR, '/biccn_metaneighbor_scores/biccn_',input$label_set,'_aurocs.csv'),sep=",",row.names=1))
        colnames(a) =rownames(a)
        a
    })

    hits <-reactive({
        MetaNeighbor::topHitsByStudy(aurocs1v1_mat(),threshold=0.8 , n_digits = 4)[,1:3]  
        
    })    
    
    clusts <-reactive({
        MetaNeighbor::extractMetaClusters(aurocs1v1_mat(), threshold=0.7)
    })

    clustScores <-reactive({
        cscores = MetaNeighbor::scoreMetaClusters(clusts(),aurocs1v1_mat())
        cscores['score'] = round(cscores['score'],4)
        tmp = strsplit(cscores$clusters,'; ')
        names(tmp) = cscores$meta_cluster
        cscores['joint_type'] = unlist(lapply(tmp, function(x){
            a = strsplit(x,'\\|')
            if (length(a) == 0) {a = list(c('outliers','outliers')) }
            a = table(unlist(lapply(a,function(x) x[2])))
            B = names(which.max(a))
            C = names(which.max(a[length(a):1]))
            if (B ==C){
                return(B)
            } else {
                return('Multiple')
            }

        
        }))

        cscores$joint_type[cscores$meta_cluster == 'outliers'] = 'outliers'
        colnames(cscores) = c('Meta Cluster','Clusters','N Studies','Score','Joint Type')

        cscores
    })


    nodes_edges <- reactive({ 
        graph_it(aurocs1v1_mat(), clusts() , low_thrs= 0.3)
    })
    # nodes_edges[[1]]$metaCls = unlist(lapply(strsplit(rownames(nodes_edges[[1]]),'\\|'),function(x) x[2]))


    output$biccn_MN_heatmap <- renderPlot({
        MetaNeighbor::plotHeatmap(aurocs1v1_mat())
    },height = 550
    )

    output$biccn_MN_network <- visNetwork::renderVisNetwork({
        
        visNetwork::visNetwork( nodes = nodes_edges()[[1]], edges = nodes_edges()[[2]] ) %>%
            visNetwork::visEdges(arrows = list( to = list(enabled =TRUE) ) )%>%
            visNetwork::visOptions(
                highlightNearest = list(enabled=TRUE, hover=T,degree=1) ,
                selectedBy = list(variable ="metaCls", main= "Choose a label", multiple = FALSE)
            )   %>% 
            visNetwork::visInteraction(
                zoomView =FALSE, 
                navigationButtons=TRUE
            ) %>% 
            visNetwork::visLegend(main = 'Meta Cluster', useGroups=FALSE, addNodes = nodes_edges()[[3]] 
            )   %>%
            visNetwork::visPhysics(
                enabled = TRUE,
                maxVelocity=10, 
                minVelocity=2,
                barnesHut = list(damping=.05),
                stabilization = list(iterations = 500) 
            )
    })


    output$biccn_MN_top_hits <- DT::renderDataTable({
        hits()
    },rownames = FALSE ,server=FALSE)


    output$biccn_MN_metaclust <- DT::renderDataTable({
        clustScores()
    },rownames = FALSE ,server=FALSE)


    ### START metamarker about section
    cmarkerset = MetaMarkers::read_meta_markers(paste0(DATADIR,'/biccn_metamarkers/biccn_MoP_class_marker_set.csv.gz'))
    cmarkerset = as.data.frame(cmarkerset)
    scmarkerset = MetaMarkers::read_meta_markers(paste0(DATADIR,'/biccn_metamarkers/biccn_MoP_subclass_marker_set.csv.gz'))
    scmarkerset = as.data.frame(scmarkerset)
    scmarkerset$cell_type = gsub('/','-',scmarkerset$cell_type)
    

    

    
    filtereddf_class <- reactive({
        df= as.data.frame(filter_markerset(cmarkerset,input$mm_class_label))
        df['log2foldchange'] = log2(df$fold_change_detection)
        df
    })

    output$class_par_hover_near <- renderUI({

        res = nearPoints(filtereddf_class(), 
            input$class_par_hover,
            xvar ='log2foldchange',
            yvar = 'auroc',
            threshold = 15,
            maxpoints = 1,
            addDist = FALSE
        )
        if (nrow(res) > 0 ){
            HTML(paste0('<b>',res$gene[1],':</b>',
                '</br>&emsp;AUROC: ',round(res$auroc[1],4), 
                '</br>&emsp;Log2(Fold Change): ',round(res$log2foldchange[1],4) ,
                '</br>&emsp;Recurrence: ',res$recurrence[1] ,
                '</br>&emsp;Is pareto: ',res$is_pareto[1] 
            ))        
        }else { 
            HTML("Hover over a gene to view the gene's data")        
        }
    })

    class_par_range <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$class_par_dblclick, {
        brush <- input$class_par_dblclick_brush
        if (!is.null(brush)) {
            class_par_range$x <- c(brush$xmin, brush$xmax)
            class_par_range$y <- c(brush$ymin, brush$ymax)

        } else {
            class_par_range$x <- NULL
            class_par_range$y <- NULL
        }

    })



    
    filtereddf_subclass <- reactive({
        df= as.data.frame(filter_markerset(scmarkerset,input$mm_subclass_label))
        df['log2foldchange'] = log2(df$fold_change_detection)
        df
    })

    output$subclass_par_hover_near <- renderUI({

        res = nearPoints(filtereddf_subclass(), 
            input$subclass_par_hover,
            xvar ='log2foldchange',
            yvar = 'auroc',
            threshold = 15,
            maxpoints = 1,
            addDist = FALSE
        )
        if (nrow(res) > 0 ){
            HTML(paste0('<b>',res$gene[1],':</b>',
                '</br>&emsp;AUROC: ',round(res$auroc[1],4), 
                '</br>&emsp;Log2(Fold Change): ',round(res$log2foldchange[1] ,4),
                '</br>&emsp;Recurrence: ',res$recurrence[1] ,
                '</br>&emsp;Is pareto: ',res$is_pareto[1] 
            ))        
        }else { 
            HTML("Hover over a gene to view the gene's data")        
        }
    })
 


    


    subclass_par_range <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$subclass_par_dblclick, {
        brush <- input$subclass_par_dblclick_brush
        if (!is.null(brush)) {
            subclass_par_range$x <- c(brush$xmin, brush$xmax)
            subclass_par_range$y <- c(brush$ymin, brush$ymax)

        } else {
            subclass_par_range$x <- NULL
            subclass_par_range$y <- NULL
        }

    })



    output$pareto_front_class <- renderPlot({

            MetaMarkers::plot_pareto_markers(cmarkerset[nrow(cmarkerset):1,], input$mm_class_label) + 
                ggplot2::theme_classic() +
                ggplot2::ylim(c(0.5,1)) + 
                ggplot2::coord_cartesian(xlim = class_par_range$x, ylim = class_par_range$y)

    })

    output$pareto_front_subclass <- renderPlot({
        MetaMarkers::plot_pareto_markers(scmarkerset[nrow(scmarkerset):1,], input$mm_subclass_label) + 
            ggplot2::theme_classic() +
            ggplot2::ylim(c(0.5,1))
    })

}  # end tab




} # end of function(session,input,output)
