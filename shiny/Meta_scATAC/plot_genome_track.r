# require(GenomicRanges)
# require(GenomicFeatures)
# require(Gviz)
# require(stringr)
# require(BSgenome.Mmusculus.UCSC.mm10)
# options(ucscChromosomeNames=FALSE)
source("global.R")
require(EnrichedHeatmap)

extractIdtrack <- function(grange) {
    # cat(file=stderr(), paste('ichr', as.character(grange@seqnames), '\n'))
    itrack <- Gviz::IdeogramTrack(genome = 'mm10', chromosome = paste0('chr', gsub('chr', '', as.character(grange@seqnames))))
    # cat(file=stderr(), 'extract ')
    # cat(file=stderr(), paste(grange@seqnames, as.character(grange@seqnames)))
    return(itrack)
}

setDbObjects <- function() {
    cat(file=stderr(), paste0("Initial set: ", Sys.time(), "\n"))    
    gtf.gene.file = './www/genome_ann/Mus_musculus.GRCm38.96.gtf'
    gtxdb <<- GenomicFeatures::makeTxDbFromGFF(gtf.gene.file, format="gtf")
    GenomeInfoDb::seqlevels(gtxdb) <- paste0('chr', GenomeInfoDb::seqlevels(gtxdb))
    gtf.tss.file = './www/genome_ann/Mus_musculus_trans_tss.GRCm38.96.gtf'
    txdb <<- GenomicFeatures::makeTxDbFromGFF(gtf.tss.file, format="gtf")
    GenomeInfoDb::seqlevels(txdb) <- paste0('chr', GenomeInfoDb::seqlevels(txdb))
    gene.dict.file = './www/genome_ann/gene_dict_sorted.tsv'
    gdf <<- read.table(gene.dict.file, header=T, sep="\t", stringsAsFactors=FALSE)
    gdf[,'chr'] <- paste0('chr', gdf[,'chr'])
    gene.dict <<- GenomicRanges::makeGRangesFromDataFrame(gdf)
    names(gene.dict) <- gdf[,'symbol']
}

writeDbObjects <- function() {
    saveDb(gtxdb, file='./www/gtxdb.txdb')
    saveDb(txdb, file='./www/txdb.txdb')
    save(gdf, gene.dict, tss.dict, file="./www/genomic.data.rds")
}

readDbObjects <- function() {
    cat(file=stderr(), paste0("Initial set: ", Sys.time(), "\n"))    
    # cat(file=stderr(), paste0(getwd(), '\n'))
    gtxdb <<- AnnotationDbi::loadDb('./www/gtxdb.txdb')
    # cat(file=stderr(), paste0('genomedb', length(transcripts(gtxdb)), '\n'))
    GenomeInfoDb::seqlevels(gtxdb) <<- paste0('chr', as.character(GenomeInfoDb::seqlevels(gtxdb)))
    txdb <<- AnnotationDbi::loadDb('./www/txdb.txdb')
    GenomeInfoDb::seqlevels(txdb) <<- paste0('chr', as.character(GenomeInfoDb::seqlevels(txdb)))
    load("./www/genomic.data.rds")
    assign("gdf", gdf, envir = .GlobalEnv)
    assign("gene.dict", gene.dict, envir = .GlobalEnv)
    assign("tss.dict", tss.dict, envir = .GlobalEnv)
    setTrackObject()
}

extractTssDict <- function() {
    tss.dict <<- Gviz::GeneRegionTrack(txdb)
    tss.dict@chromosome <<- paste0('chr', tss.dict@chromosome)
    tss.dict@genome <<- 'mm10'
}

setTrackObject <- function() {
    grtrack <<- Gviz::GeneRegionTrack(gtxdb, name = "Gene Model", transcriptAnnotation = "symbol")
    gtrack <<- Gviz::GenomeAxisTrack()
}


setRegionFromInput <- function(seq, offset) {
    vec = unlist(stringr::str_split(seq, ":"))
    loc.vec = unlist(stringr::str_split(vec[2], "-"))
    range = GenomicRanges::GRanges(vec[1], IRanges::IRanges(max(1, as.integer(loc.vec[1])-offset), as.integer(loc.vec[2])+offset))
    # cat(file=stderr(), paste(vec[1], vec[2], loc.vec[1], loc.vec[2], '\n'))
    # cat(file=stderr(), paste(GenomeInfoDb::seqlevels(gene.dict@ranges), '\n'))
    GenomeInfoDb::seqlevels(range, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(gene.dict)
    return(range)
}

setRegionFromGeneName <- function(gname, offset) {
    range = gene.dict[which(names(gene.dict) == gname),]
    cat(file=stderr(), paste(max(start(range@ranges)-offset, 1), end(range@ranges)+offset, '\n'))
    cat(file=stderr(), '\n')
    start(range@ranges) = max(start(range@ranges)-offset, 1)
    end(range@ranges) = end(range@ranges)+offset
    return(range)
}

readOverlappedGene <- function(range) {
    ogrange = IRanges::subsetByOverlaps(gene.dict, range)
    return(ogrange)
}

readOverlappedSignalCluster <- function(range, dataset, clusters, dtrack.var='horizon') {
    dtrack <- list()
    key = dataset
    signals = IRanges::subsetByOverlaps(cgr[[key]], range)
    chr <- c()
    start <- c()
    end <- c()
    if (length(signals) == 0) return(list())
    index = unlist(sapply(clusters, function(x) {
        return(which(x == cdata.label[[key]]))
    }))
    for (i in index) {
        score = c()
        for (j in 1:length(signals)) {
            base.acc = unlist(stringr::str_split(signals[j]$score, ','))
            score = c(score, base.acc[i])
        }
        score <- as.numeric(score)
        dtrack[[paste0(key, '_', i)]] <- Gviz::DataTrack(data = score, start = start(signals@ranges), end = end(signals@ranges), fontcolor.group='black', 
                                                            chromosome = GenomeInfoDb::seqnames(signals)[1], genome = 'mm10', name = paste0(cdata.label[[key]][i]), type=dtrack.var)
    }
    return(dtrack)
}

readOverlappedSignal <- function(range, index, dtrack.var='horizon') {
    dtrack <- list()
    for (key in names(gr)) {
        signals = IRanges::subsetByOverlaps(gr[[key]], range)
        chr <- c()
        start <- c()
        end <- c()
        if (length(signals) == 0) return(list())
        for (i in index) {
            score = c()
            for (j in 1:length(signals)) {
                base.acc = unlist(stringr::str_split(signals[j]$score, ','))
                score = c(score, base.acc[i])
            }
            score <- as.numeric(score)
            dtrack[[paste0(key, '_', i)]] <- Gviz::DataTrack(data = score, start = start(signals@ranges), end = end(signals@ranges), 
                                                             chromosome = GenomeInfoDb::seqnames(signals)[1], genome = 'mm10', name = paste0(key, "_", data.label[[key]][i]), type=dtrack.var)
        }
    }
    return(dtrack)
}

obtainDtrackTypeLabel <- function(dtrack.type) {
    if (dtrack.type == 'line') return('b')
    return(dtrack.type)
}

plotTargetRegionCluster <- function(dataset, clusters, sequence='', gene.name='', offset=10000, dtrack.type='line', header='output') {
    # Make sure it closes when we exit this reactive, even if there's an error
    rrtrack <- NULL
    if (nchar(gene.name) > 0) {        
        cat(file=stderr(), 'gene ')
        raw.range = setRegionFromGeneName(gene.name, 0)
        range = setRegionFromGeneName(gene.name, offset)
    } else {
        cat(file=stderr(), 'seq ')
        raw.range = setRegionFromInput(sequence, 0)
        range = setRegionFromInput(sequence, offset)
    }
    cat(file=stderr(), paste('seqname', GenomeInfoDb::seqnames(range), start(range), end(range), "?", "\n"))
    dtrack <- readOverlappedSignalCluster(range, dataset, clusters, obtainDtrackTypeLabel(dtrack.type))
    itrack <- extractIdtrack(range)
    rtrack <- Gviz::AnnotationTrack(raw.range, name='target')
    track <- list()
    track <- c(track, list(itrack))
    track <- c(track, list(gtrack))
    if (length(rtrack) > 0) track <- c(track, list(rtrack))
    track <- c(track, list(grtrack))
    track <- c(track, dtrack)
    pdf(paste0('/data/rkawaguc/data/211104_enhancer/gtrack_', header, '.pdf'), width=6, height=31)
    # png(paste0('/data/rkawaguc/data/211104_enhancer/gtrack_', header, '.png'), width=700, height=2200)
    Gviz::plotTracks(track, from=start(range), to=end(range), chromosome=as.character(GenomeInfoDb::seqnames(range)))
    # plot(g)
    dev.off()
}

plotTargetRegion <- function(sequence='', gene.name='', clusters=c(1), offset=10000, dtrack.type='line') {
    # Make sure it closes when we exit this reactive, even if there's an error
    if (nchar(gene.name) > 0) {        
        cat(file=stderr(), 'gene ')
        range = setRegionFromGeneName(gene.name, offset)
    } else {
        cat(file=stderr(), 'seq ')
        range = setRegionFromInput(sequence, offset)
    }
    cat(file=stderr(), paste('seqname', GenomeInfoDb::seqnames(range), start(range), end(range), "\n"))
    dtrack <- readOverlappedSignal(range, clusters, obtainDtrackTypeLabel(dtrack.type))
    itrack <- extractIdtrack(range)
    cat(file=stderr(), paste('seqname', GenomeInfoDb::seqnames(range), start(range), end(range), "\n"))
    rtrack <- Gviz::AnnotationTrack(range, name='target')
    track <- list()
    track <- c(track, list(itrack))
    track <- c(track, list(gtrack))
    if (length(rtrack) > 0) 
        track <- c(track, list(rtrack))
    track <- c(track, list(grtrack))
    track <- c(track, dtrack)
    return(Gviz::plotTracks(track, from=start(range), to=end(range), chromosome=as.character(GenomeInfoDb::seqnames(range))))
}

getMarkerEnrichmentHeatmap <- function(input) {
    marker = list("Major celltype"='SF', "Subtype"='SM')[[input$hm.marker]]
    dataset = dataid.list[which(dataset.list == input$hm.data)]
    fname = paste0('www/marker/gene_scluster_', marker, '_', dataset, '_jaccard.csv')
    if (file.exists(fname)) {
        table <- as.matrix(read.table(fname, header=T, sep=",", row.names=1, stringsAsFactors=F))
        colnames(table) <- sapply(colnames(table), function(x) {
            return(unlist(stringr::str_split(x, 'cluster_'))[2])
        })
        scaled.table <- apply(table, 2, function(x)(x-min(x))/(max(x)-min(x)))
        rownames(scaled.table) <- rownames(table)
        colnames(scaled.table) <- colnames(table)
        return(ComplexHeatmap::Heatmap(scaled.table, col=viridis::viridis(100)))
    }
}

extractOverlappedGenes <- function(sequence='', gene.name='', offset=10000) {
    if (nchar(gene.name) > 0) {
        range = setRegionFromGeneName(gene.name, offset)
    } else {
        range = setRegionFromInput(sequence, offset)
    }
    ogenes <- unique(names(readOverlappedGene(range)))
    table <- gdf[sapply(ogenes, function(x) { return(which(gdf[,4] == x)[1]) }),]
    return(as.data.frame(table))
}

setATACSignal <- function(fname.list) {
    tgr <- list()
    tdata.label <- list()
    for (i in 1:length(fname.list)) {
        df <- read.table(fname.list[[names(fname.list)[i]]], header=F, sep="\t")
        colnames(df) <- c('chr', 'start', 'end', 'score')
        df[,1] <- paste0('chr', gsub('chr', '', df[,1]))
        tgr[[names(fname.list)[i]]] <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
        tdata.label[[names(fname.list)[i]]] <- unlist(stringr::str_split(unlist(stringr::str_split(readLines(fname.list[[names(fname.list)[i]]], n=1), ":"))[2], ","))
    }
    return(list(tgr, tdata.label))
}

setATACSignalCluster <- function(fname.list) {
    tcgr <- list()
    tcdata.label <- list()
    for (i in 1:length(fname.list)) {
        df <- read.table(fname.list[[names(fname.list)[i]]], header=F, sep="\t")
        colnames(df) <- c('chr', 'start', 'end', 'score')
        df[,1] <- paste0('chr', gsub('chr', '', df[,1]))
        tcgr[[names(fname.list)[i]]] <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
        tcdata.label[[names(fname.list)[i]]] <- unlist(stringr::str_split(unlist(stringr::str_split(readLines(fname.list[[names(fname.list)[i]]], n=1), ":"))[2], ","))
    }
    return(list(tcgr, tcdata.label))
}

setDataFiles <- function(marker, set.global=FALSE) {
    data.dir = "/data/public_labshare/resource/Meta_scATAC/500/"
    sample.file.list = list(BICCN=file.path(data.dir, "output_BICCN2_gene_global_index__all_scanpy_obj_clust_ave_average_SM_", marker, "_500.bed"), 
                            GSE100033=file.path(data.dir, "output_GSE100033_gene_global_index_5000__all_scanpy_obj_clust_ave_average_SM_", marker, "_500.bed"),
                            GSE123576=file.path(data.dir, "output_GSE123576_gene_global_index_1000__all_scanpy_obj_clust_ave_average_SM_", marker, "_500.bed"),
                            GSE126074=file.path(data.dir, "output_GSE126074_gene_global_index__all_scanpy_obj_clust_ave_average_SM_", marker, "_500.bed"),
                            GSE111586=file.path(data.dir, "output_GSE111586_gene_global_index_5000__all_scanpy_obj_clust_ave_average_SM_", marker, "_500.bed"),
                            GSE130399=file.path(data.dir, "output_GSE1303990_gene_global_index_5000__all_scanpy_obj_clust_ave_average_SM_", marker, "_500.bed"))
    results <- setATACSignal(sample.file.list)
    if (set.global) {
        gr <<- results[[1]]
        data.label <<- results[[2]]
    } else {
        return(results)
    }
}

setDataFilesCluster <- function(dataset, set.global=FALSE) {
    data.dir = "/data/public_labshare/resource/Meta_scATAC/cluster/"
    sample.file.list = list(BICCN2=file.path(data.dir, "output_BICCN2_gene_global_index__all_scanpy_obj_clust_ave_cluster.bed"), 
                            GSE100033=file.path(data.dir, "output_GSE100033_gene_global_index_5000__all_scanpy_obj_clust_ave_cluster.bed"),
                            GSE123576=file.path(data.dir, "output_GSE123576_gene_global_index_1000__all_scanpy_obj_clust_ave_cluster.bed"),
                            GSE126074=file.path(data.dir, "output_GSE126074_gene_global_index__all_scanpy_obj_clust_ave_cluster.bed"),
                            GSE111586=file.path(data.dir, "output_GSE111586_gene_global_index_5000__all_scanpy_obj_clust_ave_cluster.bed"),
                            GSE1303990=file.path(data.dir, "output_GSE1303990_gene_global_index_5000__all_scanpy_obj_clust_ave_cluster.bed"),
                            DeepLearning=file.path(data.dir, "deep_prediction_210910.bed"))
    if (length(which(names(sample.file.list) == dataset)) == 0) return(list())
    file.list = list(sample.file.list[[which(names(sample.file.list) == dataset)]])
    names(file.list) = dataset
    results <- setATACSignalCluster(file.list)
    if (set.global) {
        cgr <<- results[[1]]
        cdata.label <<- results[[2]]
    } else {
        return(results)
    }
}

readDataFiles <- function(marker) {
    results <- readRDS(paste0('./www/signal_object/atac.signal.marker_', marker, '.rds'))
    gr <<- results[[1]]
    data.label <<- results[[2]]
}

readDataFilesCluster <- function(dataset) {
    cat(file=stderr(), paste0('./www/signal_object/atac.signal.dataset_', dataset, '.rds'))
    results <- readRDS(paste0('./www/signal_object/atac.signal.dataset_', dataset, '.rds'))
    if (length(results) == 2) {
        cgr <<- results[[1]]
        cdata.label <<- results[[2]]
    }
}

updateTrack <- function(input) {
    return(plotTargetRegion(sequence=input$region, gene.name=input$gene, cluster=c(3), dtrack.type=input$dtrack.type, offset=input$offset))
}

# preprocess the datasets (Run not on the server)

if (FALSE) {
    setDbObjects()
    extractTssDict()
    writeDbObjects()
    for (marker in celltype.list) {
        results <- setDataFiles(marker)
        saveRDS(results, paste0('atac.signal.marker_', marker, '.rds'))
    }
    cluster.label.list <<- list()
    for (dataset in dataid.list) {
        results <- setDataFilesCluster(dataset)
        if (length(results) == 2) {
            cluster.label.list[[dataset]] <- results[[2]]
            saveRDS(results, paste0('./www/atac.signal.dataset_', dataset, '.rds'))
        }
    }
    saveRDS(cluster.label.list, './www/cluster_names.rds')
}


readDbObjects()
# extractOverlappedGenes(sequence=input$region, gene.name=input$gene, offset=input$offset)
id = 'BICCN2'
readDataFilesCluster(id)
clusters = as.vector(unlist(cluster.label.list[[dataid.list[[dataset.list[1]]]]]))
enhancer <- read.table("/data/rkawaguc/data/211104_enhancer/mouse_super_enh_from_endb.tsv", header=F, sep="\t", quote="", fill=NA)
for (i in 1:dim(enhancer)[1]) {
    region = paste0(toupper(enhancer[i,5]), ":", enhancer[i,6], '-', enhancer[i,7])
    print(region)
    plotTargetRegionCluster(id, clusters, sequence=region, dtrack.type='line', offset=5000, header=paste0('enhancer_', i, '_', region))
    plotTargetRegionCluster(id, clusters, sequence=region, dtrack.type='line', offset=10000, header=paste0('enhancer_', i, '_', region, '_10000'))
    # break
}


enhancer_mat <- enhancer[,c(5, 6, 7)]
colnames(enhancer_mat) <- c('chr', 'start', 'end')
er <- makeGRangesFromDataFrame(enhancer_mat, keep.extra.columns=TRUE, ignore.strand=TRUE)

    draw(g)

enhancer_mat
colnames(enhancer_mat) = c('chr', 'start', 'end')
MAX <- max(mat[[key]][,])
MIN <- min(mat[[key]][,])
# MIN <- 0
print(mat[[key]][,])
print(c(MAX,MIN))
if (MIN == MAX) next
col_fun = colorRamp2(c(MIN, (MAX-MIN)/2, MAX), c("blue", "white", "yellow"))
if (marker) {
    top_ann = HeatmapAnnotation(enriched=anno_enriched(ylim=c(0, MAX)))
    tg <- EnrichedHeatmap(mat[[key]], name=key, col=col_fun, row_order=1:length(er_annotation), top_annotation=top_ann, axis_name=axis_name)
} else {
    top_ann = HeatmapAnnotation(enriched=anno_enriched(ylim=c(0, MAX), gp=gpar(col=c('blue', 'red'))))
    tg <- EnrichedHeatmap(mat[[key]], name=key, col=col_fun, row_order=1:length(er_annotation), row_split=er_annotation, top_annotation=top_ann, axis_name=axis_name)
}
row_order <- row_order(tg)
g <- g+tg
count <- count+1
}
if (count == 0) return()
if (!marker) {
cluster_ann = Heatmap(er_annotation, col = structure(c('red', 'blue'), names = c('valid', 'invalid')), name = "enhancer",
    show_row_names = FALSE, width = unit(3, "mm"), row_split=er_annotation, row_order=1:length(er_annotation))
g <- cluster_ann+g
}
pdf(paste0('heatmap_', subtype, '_', gse, '_', dset, '.pdf'), width=(200*count+100)/100*1.5)
draw(g)
dev.off()
row_order_list = row_order(g)    
    er <- makeGRangesFromDataFrame(enhancer_mat, keep.extra.columns=TRUE, ignore.strand=TRUE)
