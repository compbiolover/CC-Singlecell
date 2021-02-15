#DEsingle server code----
library(BiocParallel)
library(DEsingle)
library(SingleCellExperiment)

#Setting the directory for where the files are
setwd("/home/awillems/Projects/CC_Singlecell")
tumor_cells_des <- read.csv(file = "Data/GSE81861_CRC_tumor_all_cells_COUNT.csv")
nm_cells_des <- read.csv(file = "Data/GSE81861_CRC_NM_all_cells_COUNT.csv")

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


des_results <- DEsingle(counts = sce, group = condition, parallel = TRUE, BPPARAM = param)
save(des_results, "/home/awillems/Projects/CC_Singlecell/Data/Data-from-pipeline/des_results.RData")
