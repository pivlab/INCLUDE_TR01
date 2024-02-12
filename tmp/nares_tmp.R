exprs.mat <- readr::read_tsv('data/NARES_SCANfast_ComBat_with_GeneSymbol.pcl',
                              progress = FALSE)

exprs.mat=as.data.frame(exprs.mat)
rownames(exprs.mat) <- exprs.mat$GeneSymbol
exprs.mat <- as.matrix(dplyr::select(exprs.mat, -GeneSymbol))

require(PLIER)

# load PLIER pathway and cell type data
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)

# combine the pathway data from PLIER
all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, 
                                 canonicalPathways)

# what genes are common to the pathway data and the expression matrix
cm.genes <- PLIER::commonRows(all.paths, exprs.mat)

# row normalize
exprs.norm <- PLIER::rowNorm(exprs.mat)

# what should we set the minimum k parameter to in PLIER? estimate the number 
# of PC for the SVD decomposition 
set.k <- PLIER::num.pc(exprs.norm[cm.genes, ])

# PLIER main function + return results
plier.res <- PLIER::PLIER(exprs.norm[cm.genes, ], all.paths[cm.genes, ], trace = TRUE)
# 