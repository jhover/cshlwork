library(shiny)

# load coexp cons scores across 19 species, 14,131 genes
datafile = './data/coexp_cons_across_19_species.csv' 
data = read.table(file = datafile, sep = ',', header = TRUE)
genes = data$gene

# labels of plots
spe = c('human', 'chimp', 'rhesus-macaque', 'crab-eating-macaque', 'mouse', 'rat', 'rabbit', 'pig',
        'cow', 'dog', 'horse', 'goat', 'sheep', 'chicken', 'zebrafish', 'salmon', 'trout', 'fruitfly', 'roundworm')

# color palettes for plots
my_palette = c('#000004FF','#180F3EFF','#451077FF','#721F81FF','#9F2F7FFF',
               '#CD4071FF','#F1605DFF','#FD9567FF','#FEC98DFF','#FCFDBFFF')



