library(PLIER)
library(arrow)

# load PLIER pathway and cell type data
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)

gtex_data_p_path='output/gtex/GTEx_v8_gene_median_feather'
gtex_data_p=read_feather(gtex_data_p_path)

# set matrix rownames
exprs_mat=as.data.frame(gtex_data_p)

# to set the rownames as a gene symbols is better to annotate using biomaRt deals better with the isoforms 
# genes with a frequency bigger than 1: 5S_rRNA, 5_8S_rRNA and 7SK

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id_version", "external_gene_name")
filters <- "ensembl_gene_id_version"
results <- getBM(attributes = attributes, filters = filters, values = exprs_mat$gene_ens_id, mart = mart)

exprs_mat=exprs_mat %>%
  dplyr::rename(ensembl_gene_id_version=gene_ens_id)

exprs_mat = dplyr::inner_join(exprs_mat, results)

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

# exprs_mat_norm_na=na.omit(exprs_mat_norm)

PLIER(exprs_mat_norm, all_paths)

data(dataWholeBlood)
PLIER(dataWholeBlood, all_paths)

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



####

#data stands for the RNA-seq matrix and priorMat stands for the pathway matrix
#Before everything starts, make sure there is no NA values in the data
sum(is.na(exprs_mat_norm_na))==0 

cm=intersect(rownames(data), rownames(priorMat)) 
data=data[cm,]
priorMat=priorMat[cm,]


#filter our pathways with too few genes
priorMat = priorMat[,which(apply(priorMat,2,sum)>=10)]
cm = which(apply(priorMat,1,sum)>0)
data = data[cm,]
priorMat = priorMat[cm,]

#Before applying scaling and svd, make sure there is no 0-variance row
cm=which(apply(exprs_mat_norm_na,1,sd)>0)
data = data[cm,]
priorMat = priorMat[cm,]


#start rsvd/svd steps
set.seed(123456)
ns=ncol(data)
data = tscale(data)

#make sure there is no NA values after applying tscale
sum(is.na(data))==0
svdres=rsvd(data, k=min(ns, max(200, ns/4)), q=3) #You can also try svdres=svd(data)


#if there is no problem getting svdres, then set PLIER to run
PLIER.res <- PLIER(data = data, priorMat = priorMat, svdres = svdres, scale = F)