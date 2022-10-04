library(shiny)
library(shiny,quietly=T)
library(shinyWidgets,quietly=T)
library(shinydashboard,quietly=T)
library(shinydashboardPlus,quietly=T)
library(shinyBS,quietly=T)  
library(shinyjs,quietly=T)
library(corrplot)
library(rhdf5)
library(networkD3)
library(ggplot2)
library(stringr)
library(gridExtra)
library(DT)

get_celltype_labels <- function(){
  cls_label_by_sub <- c('Oligo|9', 'Astro|1', 'OPC|12', 'OPC|13', 'Micro-PVM|7', 'Endo|5', 'VLMC|15',
        'L4 IT|10', 'L2/3 L4 IT|4', 'L4 IT|9', 'L2/3 L4 IT|1', 'L2/3 IT|7',
        'L6 IT Car3|16', 'L5 IT|13', 'L5 IT|12', 'L6 IT|19', 'L6 IT Car3|17',
        'L5 ET|1', 'L5 ET|2', 'L5/6 NP|3', 'L5/6 NP|4', 'L6 CT|7', 'L6 CT|8', 'L6b|9', 'L6b|10', 'L6b|13',
        'Pvalb|10', 'Pvalb|12', 'Chandelier|1', 'Pvalb|3', 'Pvalb|6',           
        'Sst|32', 'Sst|30', 'Sst|15', 'Sst|18', 'Sst|22', 'Sst|20', 'Sst chodl|13', 'Sst|25', 'Sst|26', 'Sst|29',
        'Lamp5 Sncg|3', 'Lamp5_Lhx6|4', 'Lamp5|2', 'Pax6|5', 'Pax6 Sncg|6', 'Pax6 Sncg|7',
        'Sncg|11', 'Sncg Vip|9', 'Sncg|10', 'Vip|15', 'Vip|27', 'Vip|26', 'Vip|25', 'Vip|22', 'Vip|24', 'Vip|17')
  cls_label_by_sub2 = cls_label_by_sub[c(8:26, 27:57, 1:7)]
  return(cls_label_by_sub2)
}

get_exp_profiles <- function(id){
  file0 = './data/primate_mtg_coexp_data.h5'
  h_exp = h5read(file0, name = 'human_exp', index = list(id,1:57))
  nhp_exp = h5read(file0, name = 'nhp_exp', index = list(id,1:57))
  return(rbind(h_exp, nhp_exp))
}

get_lineage_scores <- function(id){
  file0 = './data/primate_mtg_coexp_data.h5'
  scores = h5read(file0, name = 'cell_lineage_scores_hnhp', index = list(id,1:20))
  return(scores)
}

get_exp_heatmap <- function(id){
  file0 = './data/primate_mtg_coexp_data.h5'
  scores = h5read(file0, name = 'expressolog_score', index = list(id,1:10))
  return(scores)
}

get_orthologs <- function(id){
  file0 = './data/primate_mtg_coexp_data.h5'
  g1 = h5read(file0, name = 'gene', index = list(id))
  g2 = h5read(file0, name = 'chimp_gene', index = list(id))
  g3 = h5read(file0, name = 'gorilla_gene', index = list(id))
  g4 = h5read(file0, name = 'macaque_gene', index = list(id))
  g5 = h5read(file0, name = 'marmoset_gene', index = list(id))
  return(c(g1,g2,g3,g4,g5))
}

get_score_dist <- function(id){
  file0 = './data/primate_mtg_coexp_data.h5'
  scores = h5read(file0, name = 'h_nhp_all', index = list(c(id,(id+14131)),1:4))
  return(as.vector(unlist(scores)))
}

get_bulk_pvals <- function(id){
  file0 = './data/primate_mtg_coexp_data.h5'
  pvals = h5read(file0, name = 'bulk_pvals', index = list(id,1:2))
  return(as.vector(unlist(pvals)))
}
    
get_coexp_cons_mat <- function(id){
  scorelist = unlist(data[id, 2:172])
  scores = matrix(NA, nrow = 19, ncol = 19)
  scores[lower.tri(scores, diag = F)] <- scorelist
  scores = t(scores)
  scores[lower.tri(scores, diag = F)] <- scorelist
  rownames(scores) = spe
  colnames(scores) = spe
  return(scores)
}

function(input, output) {

  output$human_nhp_exp = renderPlot({
    expmat = get_exp_profiles(match(input$gene, data$gene))
    cls_label = get_celltype_labels()
    minexp = max(0, min(expmat, na.rm = T)-1)
    maxexp = max(expmat, na.rm = T)+1
    
    plot(0, type = "n", xlab = "", main = "", ylab = "log CPM", xlim = c(1, 57), ylim = c(minexp, maxexp), xaxt = 'n')
    lines(1:57, expmat[1,], col = 'red', type = "b", pch = 20)
    lines(1:57, expmat[2,], col = '#555555', type = "b", pch = 20)
  
  
    abline(v=19.5, col="grey", lty = 5)
    abline(v=50.5, col="grey", lty = 5)
    axis(side = 1, at = 1:57, labels = cls_label, las = 2, cex.axis = 0.7)
    #legend(41, 10, legend=c('human','non-human primate'), col=c('red', '#555555'), lty=1, cex=1)
    
  },height = 300
  )
  
  output$class_explog_table = DT::renderDataTable({
      
    score_vec = get_lineage_scores(match(input$gene, data$gene))
    df_table = data.frame(Class = c('Excitatory neurons', 'Inhibitory neurons', 'Non-neurons'),
    AUROC = round(100*score_vec[c(3,4,2)])/100)
    DT::datatable(df_table, options = list(dom = 't'))
  })
  
  output$expressolog = renderPlot({
    cols = c('#3B9AB2','#47A0B6','#54A6BA','#61ACBE','#6EB2C2','#7EB8BC','#96BC9C','#AEC07B','#C6C55A','#DEC93A',
      '#E9C825','#E7C21C','#E5BC13','#E3B60B','#E1B002','#E39700','#E77800','#EA5800','#EE3900','#F21A00',
      '#3B9AB2','#47A0B6','#54A6BA','#61ACBE','#6EB2C2','#7EB8BC','#96BC9C','#AEC07B','#C6C55A','#DEC93A',
      '#E9C825','#E7C21C','#E5BC13','#E3B60B','#E1B002','#E39700','#E77800','#EA5800','#EE3900','#F21A00')
    expressolog_list = get_exp_heatmap(match(input$gene, data$gene))
    
    scores = matrix(NA, nrow = 5, ncol = 5)
    scores[lower.tri(scores, diag = F)] <- expressolog_list
    scores = t(scores)
    scores[lower.tri(scores, diag = F)] <- expressolog_list
    rownames(scores) = c('human', 'chimp', 'gorilla', 'macaque', 'marmoset')
    colnames(scores) = rownames(scores)    
    
    #  # pdf(NULL)
    corrplot(scores, method = 'color', is.corr = T, col = cols, col.lim = c(0,1),
    tl.col = 'black', tl.cex = 0.8, tl.srt = 45, na.label = ' ', na.label.col = 'white',
    cl.length = 5, cl.align.text = 'c', cl.ratio = 0.2, cl.cex = 0.75)    
  })
  
  output$exp_ortho_table = DT::renderDataTable({
  
    id = match(input$gene, data$gene)
    # get overall expressolog scores
    expressolog_list = as.vector(unlist(get_exp_heatmap(id)))
    expressolog_list <- round(100*expressolog_list)/100
    
    # get ortholog list
    orthologs = get_orthologs(id)
    primates = c('human', 'chimp', 'gorilla', 'macaque', 'marmoset')
    combos = combn(1:5, 2)
    Gene1 = paste0(primates[combos[1,]], '|', orthologs[combos[1,]])
    Gene2 = paste0(primates[combos[2,]], '|', orthologs[combos[2,]])
    
    df_table = data.frame(Gene1, Gene2, Expressolog_score = expressolog_list)
    DT::datatable(df_table, options = list(pageLength = 5))
  })
  
  output$lineage_aurocs = renderDiagonalNetwork({
    #cols = c(colorRampPalette(c('red','#bebebe'))(10), colorRampPalette(c('#bebebe','blue'))(10))
    #cols = colorRampPalette(c('white','blue'))(20)
    cols = c('#3B9AB2','#47A0B6','#54A6BA','#61ACBE','#6EB2C2','#7EB8BC','#96BC9C','#AEC07B','#C6C55A','#DEC93A',
    '#E9C825','#E7C21C','#E5BC13','#E3B60B','#E1B002','#E39700','#E77800','#EA5800','#EE3900','#F21A00')
    
    # Create fake data
    multisc <- list(name = "MTG", children = list(list(name = "Non-neurons",
                                                       children = list(list(name = "Astro, Oligo, OPC"), list(name = "MicroPVM, Endo, VLMC"))),             list(name = "Excitatory neurons", children = list(list(name = "Deep types",
                                                                                                                                                                                                                   children = list(list(name = "L5 ET, L5/6 NP"), list(name = "L6 CT, L6b"))),
                                                                                                                                                                                                              list(name = "IT types", children = list(list(name = "L2/3, L4 IT"), list(name = "L5, L6 IT"))))),
                                                  list(name = "Inhibitory neurons", children = list(list(name = "CGE-derived",
                                                                                                         children = list(list(name = "Lamp5"), list(name = "Pax6, Sncg"), list(name = "Vip"))),
                                                                                                    list(name = "MGE-derived", children = list(list(name = "Chandelier, Pvalb"),
                                                                                                                                               list(name = "Sst"), list(name = "Sst Chodl")))))))
    
    # plot for input gene    
    score_vec = get_lineage_scores(match(input$gene, data$gene))
    cids = unlist(lapply(1:20, function(jj) (findInterval(score_vec[jj], seq(0, 1, 0.05)))))
    cids[cids==21] = 20   # assign last itnerval for genes with expressolog score 1          
    colorVector <- cols[cids]      
    
    jsarray <- paste0('["', paste(colorVector, collapse = '", "'), '"]')
    nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
        
    diagonalNetwork(List = multisc, height = 800, width = 1000, fontSize = 14, linkColour = '#666666',
                  nodeColour = nodeStrokeJS, nodeStroke = nodeStrokeJS, opacity = 1)    
  }
  )  
  
  output$human_nhp_boxplots = renderPlot({
  
    gene_score_vec = get_score_dist(match(input$gene, data$gene))
    
    # get full profile over all genes        
    df_box1 = data.frame(lower_whisker = c(0.565, 0.351), q1 = c(0.821, 0.73), med = c(0.942, 0.894),
    q3 = c(0.992, 0.983), upper_whisker = c(1, 1), species = c('human', 'nhp'), scores = gene_score_vec[1:2])    
    
    df_box2 = data.frame(lower_whisker = c(0.263, 0.192), q1 = c(0.664, 0.574), med = c(0.808, 0.718),
    q3 = c(0.931, 0.881), upper_whisker = c(1, 1), species = c('human', 'nhp'), scores = gene_score_vec[3:4])
    
    df_box3 = data.frame(lower_whisker = c(0.268, 0.169), q1 = c(0.675, 0.586), med = c(0.828, 0.741),
    q3 = c(0.947, 0.903), upper_whisker = c(0.999, 0.998), species = c('human', 'nhp'), scores = gene_score_vec[5:6])
    
    df_box4 = data.frame(lower_whisker = c(0.195, 0.226), q1 = c(0.624, 0.55), med = c(0.78, 0.704),
    q3 = c(0.912, 0.874), upper_whisker = c(0.996, 0.994), species = c('human', 'nhp'), scores = gene_score_vec[7:8])
    
    p1 <- ggplot(df_box1, aes(x = as.factor(species))) +
    geom_boxplot(aes(lower = q1, upper = q3, middle = med, ymin = lower_whisker, ymax = upper_whisker),
    stat = "identity", fill = 'white', outlier.shape = NA, width = 0.7) +
    geom_hline(yintercept = 0.55, linetype = "dashed", color = "#999999") +
    geom_point(data = df_box1, aes(x = as.factor(species), y = scores), color = 'red', size = 3) + theme_bw() +
    ylim(0, 1) + ggtitle('Overall') + xlab('') + ylab('Expressolog score') +
    scale_x_discrete(labels=c("human" = "Human - NHP", "nhp" = "Across NHP"))

    p2 <- ggplot(df_box2, aes(x = as.factor(species))) +
    geom_boxplot(aes(lower = q1, upper = q3, middle = med, ymin = lower_whisker, ymax = upper_whisker),
    stat = "identity", fill = 'white', outlier.shape = NA, width = 0.7) +
    geom_hline(yintercept = 0.55, linetype = "dashed", color = "#999999") +
    geom_point(data = df_box2, aes(x = as.factor(species), y = scores), color = 'red', size = 3) + theme_bw() +
    ylim(0, 1) + ggtitle('Excitatory neurons') + xlab('') + ylab('Expressolog score') +
    scale_x_discrete(labels=c("human" = "Human - NHP", "nhp" = "Across NHP"))
    
    p3 <- ggplot(df_box3, aes(x = as.factor(species))) +
    geom_boxplot(aes(lower = q1, upper = q3, middle = med, ymin = lower_whisker, ymax = upper_whisker),
    stat = "identity", fill = 'white', outlier.shape = NA, width = 0.7) +
    geom_hline(yintercept = 0.55, linetype = "dashed", color = "#999999") +
    geom_point(data = df_box3, aes(x = as.factor(species), y = scores), color = 'red', size = 3) + theme_bw() +
    ylim(0, 1) + ggtitle('Inhibitory neurons') + xlab('') + ylab('Expressolog score') +
    scale_x_discrete(labels=c("human" = "Human - NHP", "nhp" = "Across NHP"))
    
    p4 <- ggplot(df_box4, aes(x = as.factor(species))) + 
    geom_boxplot(aes(lower = q1, upper = q3, middle = med, ymin = lower_whisker, ymax = upper_whisker), 
    stat = "identity", fill = 'white', outlier.shape = NA, width = 0.7) +
    geom_hline(yintercept = 0.55, linetype = "dashed", color = "#999999") +
    geom_point(data = df_box4, aes(x = as.factor(species), y = scores), color = 'red', size = 3) + theme_bw() +
    ylim(0, 1) + ggtitle('Non-neurons') + xlab('') + ylab('Expressolog score') +
    scale_x_discrete(labels=c("human" = "Human - NHP", "nhp" = "Across NHP"))
  
    grid.arrange(p1, p2, p3, p4, nrow = 1)
  }
  )
  
  
  output$pval_table = DT::renderDataTable({
    pval_vec = get_bulk_pvals(match(input$gene, data$gene))
    pvals = data.frame(Group1 = rep('Human - Non-human mammals',2), 
    Group2 = c('Within non-human mammals', 'Vertebrates - Non-human mammals'),
    adjusted_pval = c(formatC(pval_vec[1], format = "e", digits = 2),
    formatC(pval_vec[2], format = "e", digits = 2)))    
    DT::datatable(pvals, options = list(lengthChange = FALSE, dom = 't'), rownames = F)
  })
  
  output$vertBoxplots = renderPlot({
    mat1 = get_coexp_cons_mat(match(input$gene, data$gene))
    
    # shortlist upto zfish    
    h_m = mat1[1,2:13]  # humans with mammals
    mat2 = mat1[2:13,2:13]
    m_m = mat2[lower.tri(mat2, diag = F)]  # mammals within themselves
    mat3 = mat1
    mat3[upper.tri(mat3, diag = T)] = NA
    a_m = mat3[14:15,-1]   # animals (chicken, zfish) with non-human mammals
    
    # data.frame for ggplot2
    df_vert = data.frame(score = c(h_m, m_m, a_m), species = c(rep('human - non-human mammals', length(h_m)), rep('within non-human mammals', 66), rep('non-human mammals - vertebrates', length(a_m))))
    df_vert$species <- factor(df_vert$species, levels = c('human - non-human mammals', 'within non-human mammals', 'non-human mammals - vertebrates'), ordered = TRUE) 
    
    if(sum(df_vert$score, na.rm = T)){
      ggplot(df_vert, aes(x = species, y = score, fill = species)) +
      stat_boxplot(geom = "errorbar", width = 0.25) +
      geom_boxplot(fill = '#d3d3d3', outlier.shape = NA, width = 0.7)  +
      theme(axis.text.x=element_text(angle=45,size = 7)) + 
      geom_jitter(size = 2.5, colour = c(rep("#900c3f", length(h_m)), rep("#ffc30f", 66), rep("#ff5733", length(a_m)))) +     theme_bw() + xlab('') + ylab('Coexpression conservation') + 
      theme(text = element_text(size=12), axis.text.x = element_text(angle=0, hjust=0.5), legend.position = "none") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
    }else{
      ggplot() + theme_void() +
      geom_text(aes(0, 0, label = str_wrap("Coexpression conservation scores not available for this gene", 40))) +
      xlab(NULL)
    }
  }
  )
  
  output$aucHeat = renderPlot({
    scorelist = unlist(data[data$gene==input$gene, 2:172])
    scores = matrix(NA, nrow = 19, ncol = 19)
    scores[lower.tri(scores, diag = F)] <- scorelist
    scores = t(scores)
    scores[lower.tri(scores, diag = F)] <- scorelist
    rownames(scores) = spe
    colnames(scores) = spe
   
    if(sum(scores, na.rm = T)){
      corrplot(scores, method = 'color', is.corr = F, col = my_palette, #col.lim = c(0,1),
      tl.col = 'black', tl.cex = 0.8, tl.srt = 45, na.label = ' ', na.label.col = 'white')
    }else{
      ggplot() + theme_void() +
      geom_text(aes(0, 0, label = str_wrap("Coexpression conservation scores not available for this gene", 100))) +
      xlab(NULL)
    }
    
  }
  )
    
}
  
  


