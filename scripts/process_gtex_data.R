#!/usr/bin/env Rscript

# Processes a gene expression TSV file, renames columns, and saves it in RDS format

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 2) {
  stop("Usage: script.R <exprs_path> <output_file>", call. = FALSE)
}

# Assign arguments to variables for clarity
exprs_path <- args[1]
output_file <- args[2]
dir.create(dirname(output_file), showWarnings=FALSE)

# Load required libraries
library(dplyr)

# Load expression data
exprs_data = read.table(exprs_path, header = TRUE, sep = "\t", skip=2, check.names = FALSE) %>%
  dplyr::rename(
    gene_ens_id = Name,    # Change 'Name' to 'gene_ens_id'
    gene_symbol = Description # Change 'Description' to 'gene_symbol'
  )

# Save the processed data to a RDS file
saveRDS(exprs_data, output_file)


# Add a message at the end to indicate success
cat("File successfully written to", output_file, "\n")
