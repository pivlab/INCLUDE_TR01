#!/usr/bin/env Rscript

# Processes GTEx gene expression data using PLIER for pathway analysis, including data Z-score normalization

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 2) {
  stop("Usage: script.R <exprs_path> <output_file>", call. = FALSE)
}

# Assign arguments to variables 
gtex_expression_path <- args[1]
output_file_path <- args[2]

# Load libraries
library(PLIER)
library(dplyr)

# Load PLIER pathway and cell type data
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)

gtex_expression_path='output/gtex/GTEx_v8_gene_median_tpm.rds'

# Load GTEx data
gtex_expression_data <- readRDS(gtex_expression_path)

# Prepare output directory
dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)

# Remove gene ens id column and duplicate genes
gtex_expression_data <- subset(gtex_expression_data, select = -c(gene_ens_id))
gtex_expression_data <- gtex_expression_data[!duplicated(gtex_expression_data["gene_symbol"]),]

# Rename rows with gene symbols
rownames(gtex_expression_data) <- gtex_expression_data[,"gene_symbol"]

# Remove gene symbol column
gtex_expression_data <- subset(gtex_expression_data, select = -c(gene_symbol))

# Remove NA
gtex_expression_data = na.omit(gtex_expression_data)

# Convert to matrix
gtex_expression_matrix <- as.matrix(gtex_expression_data)

# Z-score normalization
gtex_expression_matrix <- PLIER::rowNorm(gtex_expression_matrix)

# Combine the pathway data from PLIER
all_paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways)

# Run PLIER
plier_result <- PLIER(gtex_expression_matrix, all_paths, allGenes = FALSE)

# Save results
saveRDS(plier_result, file = output_file_path)

