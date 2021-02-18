#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Code to test several scRNA-seq methods to my method

#Loading needed packages----
library(BiocParallel)
library(DESeq2)
library(DEsingle)
library(edgeR)
library(EMDomics)
library(magrittr)
library(scDD)
library(scde)
library(SINCERA)
library(SingleCellExperiment)
library(tidyverse)

#Loading needed functions----
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
#scDD method----
#Loading our single cell data for scdd
load("Data/Exported-data/R-objects/all_tumor_cells_fpkm_for_scdd.RData")
load("Data/Exported-data/R-objects/all_nm_cells_fpkm_for_scdd.RData")
# all_tumor_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_tumor_all_cells_FPKM.csv")
# rownames(all_tumor_cells_fpkm) <- all_tumor_cells_fpkm$X
# all_tumor_cells_fpkm <- gene_name_cleaner(data.to.clean = all_tumor_cells_fpkm)
# all_tumor_cells_fpkm <- t(all_tumor_cells_fpkm)
# all_tumor_cells_fpkm <- subset(all_tumor_cells_fpkm, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
#save(all_tumor_cells_fpkm, file = "Data/Exported-data/R-objects/all_tumor_cells_fpkm_for_scdd.RData")

# all_nm_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_NM_all_cells_FPKM.csv")
# rownames(all_nm_cells_fpkm) <- all_nm_cells_fpkm$X
# all_nm_cells_fpkm <- gene_name_cleaner(data.to.clean = all_nm_cells_fpkm)
# all_nm_cells_fpkm <- t(all_nm_cells_fpkm)
# all_nm_cells_fpkm <- subset(all_nm_cells_fpkm, select=c(RHC3934__Bcell__.7DEA7B:RHC6187__Macrophage__.FFFF55))
#save(all_nm_cells_fpkm, file = "Data/Exported-data/R-objects/all_nm_cells_fpkm_for_scdd.RData")

#Reading in my cell-line data----
cell_lines <- read.csv("Data/Single-cell-data/FPKM/GSE81861_Cell_Line_FPKM.csv")
rownames(cell_lines) <- cell_lines$X
cell_lines <- gene_name_cleaner(data.to.clean = cell_lines)
cell_lines <- t(cell_lines)
cell_lines<- subset(cell_lines, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cell_lines <- as.data.frame(cell_lines)
h1_and_gm_cells <- select(cell_lines, contains("H1") | contains("GM"))
h1_and_gm_cells <- select(h1_and_gm_cells, contains("RHG") | contains("RHC"))

h1_cells <- select(h1_and_gm_cells, contains("H1"))
gm_cells <- select(h1_and_gm_cells, contains("GM"))

h1_colnames <- paste("C1",seq(1:length(colnames(h1_cells))),sep = ".") 
gm_colnames <- paste("C2",seq(1:length(colnames(gm_cells))),sep = ".")
colnames(h1_cells) <- h1_colnames
colnames(gm_cells) <- gm_colnames

unique_h1_names <- unique(rownames(h1_cells))
h1_cells <- h1_cells[unique_h1_names,]

unique_gm_names <- unique(rownames(gm_cells))
gm_cells <- gm_cells[unique_gm_names,]

condition <- c(rep(1, ncol(h1_cells)), rep(2, ncol(gm_cells)))

names(condition) <- c(colnames(h1_cells), colnames(gm_cells))


i <- c(1:ncol(h1_cells))  
z <- c(1:ncol(gm_cells))

h1_cells <- apply(h1_cells, c(1,2),           
                           function(x) as.numeric(as.character(x)))

gm_cells <- apply(gm_cells, c(1,2),            
                              function(x) as.numeric(as.character(x)))


sce <- SingleCellExperiment(assays=list(normcounts=cbind(h1_cells,
                                                         gm_cells)),
                            colData=data.frame(condition))

#Filtering the sce object
sce_filtered <- preprocess(sce, zero.thresh=0.9)




prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
sce_significance_test<- scDD(sce_filtered, prior_param=prior_param, testZeroes=FALSE, categorize = FALSE)
scdd_res_cell_line <- results(sce_significance_test)
scdd_res_cell_line <- res[with(scdd_res_cell_line, order(nonzero.pvalue.adj)), ]
scdd_res_cell_line <- filter(scdd_res_cell_line, nonzero.pvalue.adj<0.05)



#Doing pre-processing to get my data ready to go into the SingleCellExperiment data for scDD method
#container format
tumor_colnames <- paste("C1",seq(1:length(colnames(all_tumor_cells_fpkm))),sep = ".")
nm_colnames <- paste("C2",seq(1:length(colnames(all_nm_cells_fpkm))),sep = ".")
colnames(all_tumor_cells_fpkm) <- tumor_colnames
colnames(all_nm_cells_fpkm) <- nm_colnames
names(condition) <- c(colnames(all_tumor_cells_fpkm), colnames(all_nm_cells_fpkm))
i <- c(1:ncol(all_nm_cells_fpkm))  
z <- c(1:ncol(all_tumor_cells_fpkm))

all_nm_cells_fpkm <- apply(all_nm_cells_fpkm, c(1,2),           
                    function(x) as.numeric(as.character(x)))

all_tumor_cells_fpkm <- apply(all_tumor_cells_fpkm, c(1,2),            
                                 function(x) as.numeric(as.character(x)))

#Making the SingleCellExperiment object of my data---- 
unique_tumor_names <- unique(rownames(all_tumor_cells_fpkm))
all_tumor_cells_fpkm <- all_tumor_cells_fpkm[unique_tumor_names,]

unique_nm_names <- unique(rownames(all_nm_cells_fpkm))
all_nm_cells_fpkm <- all_nm_cells_fpkm[unique_nm_names,]

condition <- c(rep(1, ncol(all_tumor_cells_fpkm)), rep(2, ncol(all_nm_cells_fpkm)))


sce <- SingleCellExperiment(assays=list(normcounts=cbind(all_tumor_cells_fpkm,
                                                         all_nm_cells_fpkm)),
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
scdd_res <- results(sce_significance_test)
scdd_res <- res[with(scdd_res, order(nonzero.pvalue.adj)), ]
scdd_res <- filter(scdd_res, nonzero.pvalue.adj<0.05)

#save(scdd_res,file = "Data/Exported-data/R-objects/scdd_res.RData")
#write.csv(scdd_res, file = "Data/Exported-data/Csv-files/scdd_res.csv")

#DESeq2----
#Reading in the count files
tumor_cells_ds2 <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_tumor_all_cells_COUNT.csv")
nm_cells_ds2 <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_NM_all_cells_COUNT.csv")

#DEsingle----
tumor_cells_des <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_tumor_all_cells_COUNT.csv")
nm_cells_des <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_NM_all_cells_COUNT.csv")

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
param <- MulticoreParam(workers = 2, progressbar = TRUE)
register(param)


des_results <- DEsingle(counts = sce, group = condition, parallel = TRUE, BPPARAM = param)
#save(des_results, "/home/awillems/Projects/CC_Singlecell/Data/Data-from-pipeline/des_results.RData")

#load("Data/Exported-data/R-objects/des_results.RData")



#edgeR----

#EMDomics----

#scde----

#SINCERA----
