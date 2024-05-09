#!/usr/bin/env Rscript

# Data preparation for PLIER model:
#   * z-score (PLIER::rowNorm)
#   * Determine svdres (svd)


# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 2) {
  stop("Usage: script.R <exprs_path> <output_file>", call. = FALSE)
}

# Assign arguments to variables 
expression_dataset_path <- args[1]
output_file <- args[2]

# Load libraries
`%>%` <- dplyr::`%>%`
library(PLIER)
library(dplyr)
source(here::here('scripts/plier_util.R'))

# Load PLIER pathway and cell type data
expression_dataset_path <- here::here('output/gtex/GTEx_v8_gene_median_tpm.rds')
output_file <- here::here('output/gtex/gtex_plier_prep_chr21.rds')
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)
data(immunePathways)
data(oncogenicPathways)
chr21_pathway=readRDS(here::here('output/pathways/chr21_pathway.rds'))

# Load data
expression_dataset <- readRDS(expression_dataset_path)

# Prepare output directory
output_file_path=dirname(output_file)
dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)

# Remove gene ens id column and duplicate genes
expression_dataset <- subset(expression_dataset, select = -c(gene_ens_id))
expression_dataset <- expression_dataset[!duplicated(expression_dataset["gene_symbol"]),]

# Rename rows with gene symbols
rownames(expression_dataset) <- expression_dataset[,"gene_symbol"]

# Remove gene symbol column
expression_dataset <- subset(expression_dataset, select = -c(gene_symbol))

# Remove NA
expression_dataset = na.omit(expression_dataset)

# Convert to matrix
expression_matrix <- as.matrix(expression_dataset)

# Combine the pathway data from PLIER
all_paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways,
                                immunePathways, oncogenicPathways, chr21_pathway)

output_combine_allPaths_expressionMatrix = combine_allPaths_expressionMatrix(expression_matrix, all_paths)

expression_matrix_cm=output_combine_allPaths_expressionMatrix$expression_matrix
all_paths_cm=output_combine_allPaths_expressionMatrix$all_paths

# compute rsvd/svd
set.seed(123456)
ns=ncol(expression_matrix_cm)
message("Computing SVD")
if(ns>500){
  message("Using rsvd")
  set.seed(123456);svdres=rsvd(expression_matrix_cm, k=min(ns, max(200, ns/4)), q=3)
}else{
  svdres=svd(expression_matrix_cm)
}
message("Done")

# save z-scored expression data, the prior information matrix and svdres to be supplied to PLIER::PLIER and the number of PCs

plier_data_list <- list("expression_matrix_cm" = expression_matrix_cm,
                        "all_paths_cm" = all_paths_cm,
                        "svdres" = svdres)

saveRDS(plier_data_list, file = output_file)
