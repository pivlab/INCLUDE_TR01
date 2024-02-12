library(biomaRt)

gtex_path='data/plier_rnaseq/gene_tpm_adrenal_gland.gct'
exprs_mat <- read.delim(file=gtex_path, skip=2)

# only leave the matrix with gene symbol as rownames
exprs_mat <- dplyr::select(exprs_mat, -id, -Description)
exprs_mat=as.data.frame(exprs_mat)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id_version", "external_gene_name")
filters <- "ensembl_gene_id_version"

results <- getBM(attributes = attributes, filters = filters, values = exprs_mat$Name, mart = mart)

exprs_mat=exprs_mat %>% 
  dplyr::rename(ensembl_gene_id_version=Name) 

exprs_mat=dplyr::inner_join(exprs_mat, results)

exprs_mat=exprs_mat %>% 
  dplyr::count(external_gene_name, name = "frequency") %>% 
  dplyr::left_join(exprs_mat) %>% 
  dplyr::filter(frequency <= 1) 
  
rownames(exprs_mat)=exprs_mat$external_gene_name

exprs_mat <- as.matrix(dplyr::select(exprs_mat,
                                     -external_gene_name,
                                     -frequency,
                                     -ensembl_gene_id_version))


# load PLIER pathway and cell type data
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)

# combine the pathway data from PLIER
all_paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, 
                                 canonicalPathways)

# what genes are common to the pathway data and the expression matrix
cm_genes <- PLIER::commonRows(all_paths, exprs_mat)

# row normalize
exprs_norm <- PLIER::rowNorm(exprs_mat)

# what should we set the minimum k parameter to in PLIER? estimate the number 
# of PC for the SVD decomposition 
set.k <- PLIER::num.pc(exprs.norm[cm.genes, ])

# PLIER main function + return results
plier.res <- PLIER::PLIER(exprs.norm[cm.genes, ], all.paths[cm.genes, ], trace = TRUE)


