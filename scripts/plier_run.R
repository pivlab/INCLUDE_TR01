#!/usr/bin/env Rscript

# Data run PLIER model:
#   * log2 transformation
#   * z-score (PLIER::rowNorm)
#   * Determine svdres (svd)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 4) {
  stop("Usage: script.R <exprs_path> <output_file> <parameter_k> <frac>", call. = FALSE)
}

# Assign arguments to variables 
plier_data_list_path <- args[1]
output_file <- args[2]
parameter_k <- as.numeric(args[3])
frac <- as.numeric(args[4])

plier_data_list_path <- 'output/gtex/gtex_plier_prep_chr21.rds'
output_file <- 'output/gtex/plier_result_k1_frac0.7_chr21.rds'
parameter_k <- 1
frac <- 0.7

# Load libraries
`%>%` <- dplyr::`%>%`
library(PLIER)
library(dplyr)

# Load data
plier_data_list=readRDS(plier_data_list_path)
expression_matrix_cm=plier_data_list$expression_matrix_cm
all_paths_cm=plier_data_list$all_paths_cm
svdres=plier_data_list$svdres

# compute k
k=num.pc(svdres)*2
k <- min(k, floor(ncol(expression_matrix_cm)*0.9))
k = k * parameter_k
k = round(k, 0)
message("k is set to ", k)

# Run PLIER (with common genes)
plier_result=PLIER::PLIER(data=expression_matrix_cm, priorMat=all_paths_cm , svdres=svdres, k=k, frac=frac, scale=FALSE)

# Prepare output directory
output_file_path=dirname(output_file)
dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)

# Save results
saveRDS(plier_result, file = output_file)
