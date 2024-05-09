#!/usr/bin/env Rscript

# PLIER to Pickle

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that we have the correct number of arguments
if (length(args) != 1) {
  stop("Usage: script.R <input_rds_plier_model>", call. = FALSE)
}

# Assign arguments to variables for clarity
input_rds_plier_model <- args[1]
input_rds_plier_model=here::here(input_rds_plier_model)

# Load libraries
library(reticulate)

# Define functions
save_as_pickle <- function(object, filename, save_directory) {
  full_path <- file.path(save_directory, filename)
  py_save_object(r_to_py(object), full_path)
}

PLIER_model_to_pickle = function(PLIER_model, save_directory){
    
    # Check if the directory exists, create if it does not
    if (!dir.exists(save_directory)) {
      dir.create(save_directory, recursive = TRUE)
    }
    
    # Assuming gtex_tmp_1 is a list with various data types
    names_list <- names(PLIER_model)
    
    for (name in names_list) {
      element <- PLIER_model[[name]]
      if (is.matrix(element) || is.array(element)) {
        # Convert matrices/arrays to data frames before saving
        df <- as.data.frame(element)
        save_as_pickle(df, paste0(name, ".pkl"), save_directory)
      } else {
        # Save other data types directly
        save_as_pickle(element, paste0(name, ".pkl"), save_directory)
      }
    }  
}

# Define your function
PLIER_model <- readRDS(input_rds_plier_model)
save_directory <- gsub("\\.rds$", "", input_rds_plier_model)
PLIER_model_to_pickle(PLIER_model, save_directory)

