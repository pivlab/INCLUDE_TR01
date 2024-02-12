library(targets)
library(tarchetypes)
library(here)
library(googlesheets4)
library(dplyr, warn.conflicts=FALSE)
options(tidyverse.quiet=TRUE)

source(here("R", "mutational_signatures.R"))
source(here("R", "genomics_regulation.R"))
source(here("R", "motifs.R"))
source(here("R", "h3k27ac.R"))
source(here("R", "global.R"))
source(here('R','sv_functions.R'))

tar_option_set(
  packages = c("RColorBrewer", "gplots", "pheatmap", "glmnet", "rsvd", "qvalue", "devtools")
  )

set.seed(20191218)

list(
  ##..............................................................................
  ## PLIER in GTEx                                                            ####
  ##..............................................................................
  tar_target()
  
)