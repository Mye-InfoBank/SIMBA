#!/usr/bin/env Rscript
"
compute_rogue_entropy.R

Usage:
    compute_rogue_entropy.R --atlas_path=<atlas_path> --variable_entropy=<variable_entropy> [options]

Mandatory options:
    --atlas_path=<atlas_path>               The path to the intermediate atlas object in h5ad format
    --variable_entropy=<variable_entropy>   The variable that defines the groups for which the entropy will be computed          


Optional options:
    --variable_grouping=<variable_grouping> The variable to group the data by. If absent, no groups will be considered
" -> doc

suppressMessages(library(conflicted))
suppressMessages(library(docopt))
suppressMessages(library(ROGUE))
suppressMessages(library(tidyverse))
suppressMessages(library(anndata))
suppressMessages(library(Matrix))
suppressMessages(library(MatrixExtra))

arguments <- docopt(doc, version = "0.1")
print(arguments)

h5ad.file.path <- arguments$atlas_path
variable.entropy <- arguments$variable_entropy
variable.grouping <- ''

if(!is.null(arguments$variable_grouping)){
  variable.grouping <- arguments$variable_grouping
} else {
  variable.grouping <- 'identity'
}

ad <- read_h5ad(h5ad.file.path)
matrix <- MatrixExtra::as.csc.matrix(ad$T$X)
metadata <- ad$obs
metadata$identity <- 'integrated_dataset'

# Entropy of genes
entropy.genes <- SE_fun(matrix) 

SEplot(entropy.genes) + ggtitle('Entropy of genes vs log mean expression')


# Entropy of genes, across conditions
for(c in unique(metadata[[variable.entropy]])){
  
  cur.cells <- rownames(metadata[metadata[[variable.entropy]] == c, ])
  
  cur.entropy.genes <- entropy.genes <- SE_fun(matrix[, cells])
  
  SEplot(entropy.genes) + ggtitle('Entropy of genes vs log mean expression', 
                                  subtitle = paste0('Subset ', variable_entropy, ' = ', c))
}


# Entropy across the whole dataset
rogue.value <- CalculateRogue(entropy.genes, platform = "UMI")

cat('Entropy across the whole dataset: ', rogue.value)

# Entropy by variable 
rogue.res <- rogue(matrix_right, labels = metadata[[variable.entropy]], 
                   samples = metadata[[variable.grouping]], platform = "UMI", span = 0.5)

rogue.boxplot(rogue.res)

if(is.null(variable.grouping)){
  saveRDS(rogue.res, file=paste0('entropy_', variable.entropy, '.rds'))
} else {
  saveRDS(rogue.res, file=paste0('entropy_', variable.entropy, ' ', variable.grouping, '.rds'))
}
