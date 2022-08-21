.libPaths("/usr/lib64/R/library")

library(shiny,quietly=T)
library(shinyalert,quietly=T)

library(shinyWidgets,quietly=T) 
library(shinydashboard,quietly=T) 
library(shinydashboardPlus,quietly=T)
library(shinyBS,quietly=T)  # bootstrap

shinyalert::useShinyalert()
shiny::enableBookmarking(store = "url")

DATADIR = '/home/ftp/data/FAIRit'
AVAILABLE_PROJECTS = dir(DATADIR,pattern = 'h5ad')
AVAILABLE_PROJECTS = sub('.h5ad','',AVAILABLE_PROJECTS)

SEARCH_DF = read.delim(paste0(DATADIR,'/metadata/search_table.tsv'),sep = "\t",row.names = 1)
SEARCH_DF$pdate = as.Date(SEARCH_DF$pdate, format = '%B %d, %Y')
# SEARCH_DF$pdate[is.na(SEARCH_DF$pdate)] = as.Date('2021-01-01')
# as.Date(format(SEARCH_DF$pdate ,'%Y %b %d'))
# check global auroc for SRP126648, SRP172768

# total counts auroc > .6 implies smartseq

# default candidates - SRP179101

datasets2drop = c(
    'SRP278583', # multi tissue T cells
    'SRP066963', # too few cells
    'SRP308387',# too few cells
    'SRP239491',# too few cells
    'SRP150863',# too few cells
    'SRP142629',# too few cells
    'SRP150630',# too few cells
    'SRP071876',# too few cells
    'SRP308826',# too few cells
    'SRP066314',# too few cells
    'SRP066317',
    'SRP228572',# too few cells
    'SRP108034',# too few cells
    'SRP303200'# too few cells
)

AVAILABLE_PROJECTS = setdiff(AVAILABLE_PROJECTS, datasets2drop)

BIOPROJ_URL = 'https://www.ncbi.nlm.nih.gov/bioproject/'
GEO_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
LABSHARE_URL = 'https://labshare.cshl.edu/shares/gillislab/resource/FAIRit/'

# obs column title to display

# discrete 
DISCRETE_LABS=c(
    # 'Batch' = "batch", 
    'Experiment ID' = "exp_id",
    'Run ID' = "run_id", 
    'Sample ID' = "samp_id",
    'Technique' ="tech"
)


CTS_LABS =c(
    'log(total counts)' = "total_counts",
    'log(# unique genes)'="n_genes_by_counts",
    'Correlation to mean' = "corr_to_mean",
    'Gini coefficient (depth)' = "gini",
    # 'Max correlation with others' = "max_corr_with_others",
    'Ribosomal genes' = "ribo",
    'Mitochondria genes' = "mt",
    'Essential genes' = "essential",
    # 'Cell cycle genes' = "cell_cycle",
    'Housekeeping genes' = "housekeeping",  
    'Female marker - Xist' = "female",
    'Male marker - Ddx3y' = "male",  
    'Highly variable genes' = "highly_variable",
    # 'Top 50 genes' = "top_50_gene" ,
    'Top 100 genes' = "top_100_gene" 
    # 'Top 200 genes' = "top_200_gene" ,
    # 'Top 500 genes' = "top_500_gene" 
)

# "min_corr_with_others",
# "total_counts_cell_cycle"      "total_counts_essential"      
# "total_counts_female"          "total_counts_highly_variable"
# "total_counts_housekeeping"    "total_counts_male"           
# "total_counts_mt"              "total_counts_ribo"  


# "log1p_corr_to_mean"  ,
# "log1p_gini",
# "pct_counts_cc_cluster_1"      "pct_counts_cc_cluster_2"     
# "pct_counts_cc_cluster_3"      "pct_counts_cc_cluster_4"     
# "pct_counts_cc_cluster_5"      "pct_counts_cc_cluster_6"     
# "pct_counts_cc_cluster_7"      "pct_counts_cc_cluster_8"    
# "total_counts_cc_cluster_1"    "total_counts_cc_cluster_2"   
# "total_counts_cc_cluster_3"    "total_counts_cc_cluster_4"   
# "total_counts_cc_cluster_5"    "total_counts_cc_cluster_6"   
# "total_counts_cc_cluster_7"    "total_counts_cc_cluster_8"  


# deprecated
OBSNAMES = c(
    'Total Counts'="total_counts",
    'Number of Genes'="n_genes_by_counts",
    'Correlation to mean' = "corr_to_mean", 
    'Gini Coefficient' = "gini",
    'Cell Cycle genes'= "cell_cycle"       ,
    'Essential Genes'="essential"         ,
    'Xist'="female"            ,
    'Ddx3y'="male"             ,
    'Mitochondia Genes'="mt"                ,
    'Ribosomal Genes'="ribo"              
    # "cc_cluster_1" ,     
    # "cc_cluster_2"  ,    
    # "cc_cluster_3"   ,  
    # "cc_cluster_4"    ,  
    # "cc_cluster_5"     , 
    # "cc_cluster_6"     ,
    # "cc_cluster_7"      ,
    # "cc_cluster_8"      ,
    # "in_top_100_genes" ,
    # "in_top_200_genes"  ,
    # "in_top_500_genes"  ,
    # "in_top_50_genes"  
)

# deprecated
OBSUNITS = c(
    'Total Counts' = 'total_counts_',
    'Log1p' = 'log1p_total_counts_',
    'Percent' = 'pct_counts_'
)

# deprectated
CLASSNAMES = c('GABAergic','Glutamatergic','Non-neuronal','unassigned')

CLASS_COLORS = data.frame(
    label = c(
        'GABAergic','Glutamatergic','Non-neuronal','unassigned'), 
    color = c(
        '#F05A28','#00ADEE','#000000','#D3D3D3')
)

SUBCLASS_COLORS = data.frame( 
    label = c(
        "Lamp5","Sncg","Vip","Sst Chodl","Sst","Pvalb" ,"L2-3 IT",
        "L4-5 IT","L5 IT","L6 IT","L6 IT Car3","L5 PT","L5-6 NP","L6 CT",
        "L6b",'Meis2','Oligo' , 'Astro','Endo', 'VLMC','SMC',  'Peri'  ,'Micro-PVM','unassigned','Other'),
    color =c(
        '#DA808C','#D633FF','#B864CC','#ECD704','#FF9900','#D93137','#C4EC04','#09CCC6','#50B2AD','#A19922', '#5100FF',
        '#0D5B78','#3E9E64','#2D8CB8','#53377D','#C60C0F','#2E3E39','#665C47', '#8D6C62','#697255', '#807059', '#665547','#94AF97','#D3D3D3','#000000'
    )
)
# pageLength

# assortment of 20 distinct colors
# DISTINCT_COLORS = c(
#     '#ff4500','#0000cd','#2e8b57', '#8b0000',
#     '#ffa500','#ffff00','#191970','#7fff00',
#     '#ba55d3','#00ffff','#98fb98','#ff00ff','#1e90ff',
#     '#87cefa','#dda0dd','#f08080','#ff1493','#808000','#2f4f4f','#ffe4b5')

# assortment of 30 distinct colors
DISTINCT_COLORS = c(
    '#0000cd',
    '#7f0000',
    '#006400',
    # '#483d8b',
    '#7f007f',
    '#008b8b',
    '#9acd32',
    '#00008b',
    '#808000',
    '#8fbc8f',
    '#696969',
    '#b03060',
    '#ff0000',
    '#ff8c00',
    '#ffd700',
    '#00ff00',
    '#9400d3',
    '#00ff7f',
    '#dc143c',
    '#00ffff',
    '#00bfff',
    '#f4a460',
    '#0000ff',
    '#f08080',
    '#da70d6',
    '#d8bfd8',
    '#1e90ff',
    '#f0e68c',
    '#90ee90',
    '#ff1493',
    '#7b68ee'
)


N_COLORS  = length(DISTINCT_COLORS)

ALL_AUROCS = read.delim(paste0(DATADIR,'/metadata/auc_wrt_global.tsv'),sep="\t", row.names=1)
# ALL_AUROCS[,c('ribo','mt','essential','cell_cycle','housekeeping',)]
# BOOKMARK LINKS 

# https://gillisweb.cshl.edu/MetaQC/?_inputs_&projid=%22SRP167086%22


#    store: Either ‘"url"’, which encodes all of the relevant values in a
#           URL, ‘"server"’, which saves to disk on the server, or
#           ‘"disable"’, which disables any previously-enabled
#           bookmarking.



# In general, an AUC of 0.5 suggests no discrimination 
# 0.7 to 0.8 is considered acceptable, 
# 0.8 to 0.9 is considered excellent, 
# 0.9 is considered outstanding


# to use: 
#   paste(TEXT_OPTIONS[[list_name]] ,collapse = '')
# TEXT_OPTIONS = list(
#     total_counts = list('Comparing total counts for each cell, with the global distribution, we get an AUROC of ' ,auc, ' suggesting that this dataset', adj ,'likely used ' ,tech_used, '.' ),
#     n_genes_by_counts = list('Comparing the number of unique genes for each cell, with the global distribution, we get an AUROC of ' ,auc, ' which is')
# ) 