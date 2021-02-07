#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Code to test several scRNA-seq methods to my method

#Loading needed packages----
library(MAST);packageVersion("MAST")
library(scDD);packageVersion("scDD")
library(SingleCellExperiment);packageVersion("SingleCellExperiment")

#Loading our single cell data----
all_tumor_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_tumor_all_cells_FPKM.csv")
all_tumor_cells_fpkm <- subset(all_tumor_cells_fpkm, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
all_nm_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_NM_all_cells_FPKM.csv")
all_nm_cells_fpkm <- subset(all_nm_cells_fpkm, select=c(RHC3934__Bcell__.7DEA7B:RHC6187__Macrophage__.FFFF55))

#Doing pre-processing to get my data ready to go into the SingleCellExperiment data----
#container format
condition <- c(rep(1, ncol(all_tumor_cells_fpkm)), rep(2, ncol(all_nm_cells_fpkm)))
tumor_colnames <- paste("CC",seq(1:length(colnames(all_tumor_cells_fpkm))),sep = ".")
nm_colnames <- paste("NM",seq(1:length(colnames(all_nm_cells_fpkm))),sep = ".")
colnames(all_tumor_cells_fpkm) <- tumor_colnames
colnames(all_nm_cells_fpkm) <- nm_colnames
names(condition) <- c(colnames(all_tumor_cells_fpkm), colnames(all_nm_cells_fpkm))

#Making the SingleCellExperiment object of my data---- 
sce <- SingleCellExperiment(assays=list(normcounts=cbind(all_tumor_cells_fpkm,
                                                         all_nm_cells_fpkm)),
                            colData=data.frame(condition))

sce_filtered <- preprocess(sce, zero.thresh=0.9)

prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
sce_significance_test<- scDD(sce_filtered, prior_param=prior_param, testZeroes=FALSE, categorize = FALSE)
res <- results(sce_significance_test)
res <- res[with(res, order(nonzero.pvalue.adj)), ]
res <- filter(res, nonzero.pvalue.adj<0.05)
