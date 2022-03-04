#DEsingle server code
#Author: Andrew Willems <awillems@vols.utk.edu
#Purpose: To evaluate the DESingle method on my scRNA-seq data

#Loading needed packages----
library(BiocParallel)
library(DEsingle)
library(SingleCellExperiment)
library(tidyverse)

#Setting the directory for where the files are----
#setwd("/home/awillems/Projects/CC_Singlecell")

#Loading in the CRC raw files----
tumor_cells_des <- read.csv(file = "GSE81861_CRC_tumor_all_cells_COUNT.csv")
nm_cells_des <- read.csv(file = "GSE81861_CRC_NM_all_cells_COUNT.csv")

#Gene name cleaner function loaded----
gene_name_cleaner <- function(data.to.clean=all_tumor_cells_fpkm_denoised_df){
  data.to.clean <-t(data.to.clean)
  current_colname_split <- strsplit(colnames(data.to.clean), "_")
  finished_gene_list <- c()
  current_list <- current_colname_split
  for (y in seq(1:length(current_list))){
    finished_gene_list <- c(finished_gene_list, current_list[[y]][2])
  }
  colnames(data.to.clean) <- finished_gene_list
  return(data.to.clean)
}

#CRC patient pre-processing----
rownames(tumor_cells_des) <- tumor_cells_des$X
tumor_cells_des <- gene_name_cleaner(data.to.clean = tumor_cells_des)
tumor_cells_des <- t(tumor_cells_des)
tumor_cells_des <- subset(tumor_cells_des, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))

rownames(nm_cells_des) <- nm_cells_des$X
nm_cells_des <- gene_name_cleaner(data.to.clean = nm_cells_des)
nm_cells_des <- t(nm_cells_des)
nm_cells_des <- subset(nm_cells_des, select=c(RHC3934__Bcell__.7DEA7B:RHC6187__Macrophage__.FFFF55))

#Making the condition vector
condition <- c(rep(1, ncol(tumor_cells_des)), rep(2, ncol(nm_cells_des)))
condition <- factor(condition)

tumor_colnames <- paste("C1",seq(1:length(colnames(tumor_cells_des))),sep = ".")
nm_colnames <- paste("C2",seq(1:length(colnames(nm_cells_des))),sep = ".")
colnames(tumor_cells_des) <- tumor_colnames
colnames(nm_cells_des) <- nm_colnames
names(condition) <- c(colnames(tumor_cells_des), colnames(nm_cells_des))

i <- c(1:ncol(nm_cells_des))
z <- c(1:ncol(tumor_cells_des))

nm_cells_des <- apply(nm_cells_des, c(1,2),
                      function(x) as.numeric(as.character(x)))

tumor_cells_des <- apply(tumor_cells_des, c(1,2),
                         function(x) as.numeric(as.character(x)))

sce <- SingleCellExperiment(assays=list(counts=cbind(tumor_cells_des,
                                                     nm_cells_des)))

#Setting up MacOS/Linux specific multi-core parameters for this method
param <- MulticoreParam(workers = 12, progressbar = TRUE)
register(param)

#CRC patients DESingle----
des_results <- DEsingle(counts = sce, group = condition, parallel = TRUE, BPPARAM = param)
save(des_results, file = "des_results.RData")

# #Cancer cell lines (H1 vs GM)----
# #Reading in my cell-line data----
# cell_lines <- read.csv("Data/GSE81861_Cell_Line_COUNT.csv")
# 
# #Pre-processing for the cell lines----
# rownames(cell_lines) <- cell_lines$X
# cell_lines <- gene_name_cleaner(data.to.clean = cell_lines)
# cell_lines <- t(cell_lines)
# cell_lines<- subset(cell_lines, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
# cell_lines <- as.data.frame(cell_lines)
# h1_and_gm_cells <- select(cell_lines, contains("H1") | contains("GM"))
# h1_and_gm_cells <- select(h1_and_gm_cells, contains("RHG") | contains("RHC"))
# 
# h1_cells <- select(h1_and_gm_cells, contains("H1"))
# gm_cells <- select(h1_and_gm_cells, contains("GM"))
# 
# h1_colnames <- paste("C1",seq(1:length(colnames(h1_cells))),sep = ".") 
# gm_colnames <- paste("C2",seq(1:length(colnames(gm_cells))),sep = ".")
# colnames(h1_cells) <- h1_colnames
# colnames(gm_cells) <- gm_colnames
# 
# unique_h1_names <- unique(rownames(h1_cells))
# h1_cells <- h1_cells[unique_h1_names,]
# 
# unique_gm_names <- unique(rownames(gm_cells))
# gm_cells <- gm_cells[unique_gm_names,]
# 
# condition <- c(rep(1, ncol(h1_cells)), rep(2, ncol(gm_cells)))
# condition <- factor(condition)
# 
# names(condition) <- c(colnames(h1_cells), colnames(gm_cells))
# 
# 
# i <- c(1:ncol(h1_cells))  
# z <- c(1:ncol(gm_cells))
# 
# h1_cells <- apply(h1_cells, c(1,2),           
#                   function(x) as.numeric(as.character(x)))
# 
# gm_cells <- apply(gm_cells, c(1,2),            
#                   function(x) as.numeric(as.character(x)))
# 
# sce <- SingleCellExperiment(assays=list(counts=cbind(h1_cells,
#                                                     gm_cells)))
# 
# #Setting up MacOS/Linux specific multi-core parameters for this method
# param <- MulticoreParam(workers = 12, progressbar = TRUE)
# register(param)
# 
# #DESingle for H1 vs. GM cell lines----
# des_results_h1_vs_gm_lines <- DEsingle(counts = sce, group = condition, parallel = TRUE, BPPARAM = param)
# save(des_results_h1_vs_gm_lines, file = "Data/des_results_h1_vs_gm_cell_lines.RData")
# results.sig.h1.gm <- des_results_h1_vs_gm_lines[des_results_h1_vs_gm_lines$pvalue.adj.FDR < 0.05, ]
# save(results.sig.h1.gm, file = "Data/des_results_h1_vs_gm_cell_lines_sig.RData")