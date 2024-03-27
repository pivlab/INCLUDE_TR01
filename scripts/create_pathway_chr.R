#!/usr/bin/env Rscript

# Create chr21 pathway for PLIER

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 1) {
  stop("Usage: script.R <output_file>", call. = FALSE)
}

# Assign arguments to variables 
output_file <- args[1]

# Prepare output directory
output_file_path=dirname(output_file)
dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)

# Load libraries
library(biomaRt)
library(dplyr)

# Initialize the connection to the Ensembl database for human genes
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve a dataframe of gene information including symbols and positions
genes <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'band', 'start_position', 'end_position'), mart = mart)

# Define a vector of canonical chromosome names
canonical_chr <- paste0('chr', c(1:22, 'X', 'Y'))

# Filter genes specifically for chromosome 21
chr21_genes <- genes %>%
  mutate(chromosome_name=paste0('chr', chromosome_name)) %>%
  filter(hgnc_symbol != '', chromosome_name == 'chr21') %>%
  select(hgnc_symbol, chromosome_name)

# List of gene symbols from chr21_genes
unique_chr21_genes <- unique(chr21_genes$hgnc_symbol)

# Create a matrix to indicate gene presence on chromosome 21
m_chr21_genes <- matrix(0, nrow = length(unique_chr21_genes), ncol = 1, dimnames = list(unique_chr21_genes, "chr21"))

# Update the matrix for genes present on chromosome 21
m_chr21_genes[chr21_genes$hgnc_symbol, "chr21"] <- 1

# Further filter genes on chromosome 21 for band information
chr21_bands_genes <- genes %>%
  mutate(chromosome_name=paste0('chr', chromosome_name), chr_band=paste0(chromosome_name, '_', band)) %>%
  filter(hgnc_symbol != '', chromosome_name == 'chr21') %>%
  select(hgnc_symbol, chr_band)

# Create a matrix to represent gene-band pairings, initialized to 0
m_chr21_bands_genes <- matrix(0, nrow = length(unique(chr21_bands_genes$hgnc_symbol)), 
                              ncol = length(unique(chr21_bands_genes$chr_band)), 
                              dimnames = list(unique(chr21_bands_genes$hgnc_symbol), unique(chr21_bands_genes$chr_band)))

# Populate the matrix with 1 for existing gene-band pairs
for (i in 1:nrow(chr21_bands_genes)) {
  m_chr21_bands_genes[chr21_bands_genes$hgnc_symbol[i], chr21_bands_genes$chr_band[i]] <- 1
}

# Convert matrices to data frames and prepare for merging
df_chr21_bands_genes <- as.data.frame(m_chr21_bands_genes)
df_chr21_genes <- as.data.frame(m_chr21_genes)
df_chr21_bands_genes$Gene <- rownames(df_chr21_bands_genes)
df_chr21_genes$Gene <- rownames(df_chr21_genes)

# Merge data frames on gene names and replace NA with 0
merged_df <- merge(df_chr21_bands_genes, df_chr21_genes, by = "Gene", all = TRUE)
merged_df[is.na(merged_df)] <- 0

# Set genes as row names and remove the now redundant Gene column
rownames(merged_df) <- merged_df$Gene
merged_df <- merged_df[, -which(names(merged_df) == "Gene")]

# Optionally convert back to a matrix for further analysis
chr21_pathway <- as.matrix(merged_df)

saveRDS(chr21_pathway, file=output_file)
