#Loading needed packages----
library(scDD)
library(SingleCellExperiment)
library(tidyverse)

#Loading data----
load("all_tumor_cells_fpkm_subset.RData")
load("all_nm_cells_fpkm_subset.RData")

#Doing pre-processing to get my data ready to go into the SingleCellExperiment data----
#container format
condition <- c(rep(1, ncol(all_tumor_cells_fpkm_subset)), rep(2, ncol(all_nm_cells_fpkm_subset)))

#Making the SingleCellExperiment object of my data---- 
sce <- SingleCellExperiment(assays=list(normcounts=cbind(all_tumor_cells_fpkm_subset,
                                                         all_nm_cells_fpkm_subset)),
                            colData=data.frame(condition))

#Filtering the sce object
sce_filtered <- preprocess(sce, zero.thresh=0.9)


#Trying to run the code in serial to determine what is going on with scDD function----
#BiocParallel::register(BiocParallel::SerialParam())

#For multiple cores----
BiocParallel::register(BiocParallel::MulticoreParam())

#scDD function and results----
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
sce_significance_test<- scDD(sce_filtered, prior_param=prior_param, testZeroes=FALSE, categorize = FALSE)
res <- results(sce_significance_test)
res <- res[with(res, order(nonzero.pvalue.adj)), ]
res <- filter(res, nonzero.pvalue.adj<0.05)