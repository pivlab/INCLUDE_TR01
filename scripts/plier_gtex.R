#!/usr/bin/env Rscript

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 2) {
  stop("Usage: script.R <exprs_path> <output_file>", call. = FALSE)
}

# Assign arguments to variables for clarity
gtex_data_p_path <- args[1]
output_file <- args[2]

library(PLIER)
library(arrow)
library(biomaRt)
library(dplyr)

# load PLIER pathway and cell type data
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)

gtex_data_p=read_feather(gtex_data_p_path)

dir.create(dirname(output_file), showWarnings=F)

# set matrix rownames
exprs_mat=as.data.frame(gtex_data_p)

# to set the rownames as a gene symbols is better to annotate using biomaRt deals better with the isoforms 
# genes with a frequency bigger than 1: 5S_rRNA, 5_8S_rRNA and 7SK

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id_version", "external_gene_name")
filters <- "ensembl_gene_id_version"
results <- getBM(attributes = attributes, filters = filters, values = exprs_mat$gene_ens_id, mart = mart)

exprs_mat=exprs_mat %>%
  dplyr::rename(ensembl_gene_id_version=gene_ens_id)

exprs_mat = dplyr::inner_join(exprs_mat, results)

# remove duplicates
exprs_mat=exprs_mat %>%
  dplyr::count(external_gene_name, name = "frequency") %>%
  dplyr::left_join(exprs_mat) %>%
  dplyr::filter(frequency <= 1)

rownames(exprs_mat)=exprs_mat$external_gene_name

exprs_mat <- as.matrix(dplyr::select(exprs_mat,
                                     -external_gene_name,
                                     -frequency,
                                     -ensembl_gene_id_version,
                                     -gene_symbol))


# remove na
exprs_mat_norm=na.omit(exprs_mat_norm)

# combine the pathway data from PLIER
all_paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, 
                                 canonicalPathways)

# what genes are common to the pathway data and the expression matrix
cm_genes <- PLIER::commonRows(all_paths, exprs_mat)

# row normalize
exprs_mat_norm <- PLIER::rowNorm(exprs_mat)

# what should we set the minimum k parameter to in PLIER? estimate the number 
# of PC for the SVD decomposition 
set_k <- PLIER::num.pc(exprs_mat_norm[cm_genes, ])

# PLIER main function + return results
plier_res <- PLIER::PLIER(exprs_mat_norm[cm_genes, ], all_paths[cm_genes, ])

saveRDS(plier_res, file = output_file)

