#!/usr/bin/env Rscript

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 2) {
  stop("Usage: script.R <exprs_path> <output_file>", call. = FALSE)
}

# Assign arguments to variables for clarity
exprs_path <- args[1]
output_file <- args[2]

# Load required libraries
library(dplyr)
library(arrow)
library(readr) # Explicitly load readr for read_tsv

# Load expression data
# Adding comment: Skip the first two lines of the file and read the rest
exprs_data <- read_tsv(exprs_path, skip = 2) %>%
  # Rename columns for clarity and consistency
  rename(
    gene_ens_id = Name,    # Change 'Name' to 'gene_ens_id'
    gene_symbol = Description # Change 'Description' to 'gene_symbol'
  )

# Save the processed data to a Feather file
# Comment: Write the modified data frame to a Feather file for efficient storage
write_feather(exprs_data, output_file)

# Add a message at the end to indicate success
cat("File successfully written to", output_file, "\n")



