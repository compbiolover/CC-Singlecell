#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Code to test several scRNA-seq methods to my method

#Loading needed packages----
library(BiocParallel)
library(cowplot)
library(DESeq2)
library(DEsingle)
library(EBImage)
library(edgeR)
library(EMDomics)
library(grid)
library(gridExtra)
library(magrittr)
library(png)
library(scDD)
library(scde)
library(scran)
library(SINCERA)
library(SingleCellExperiment)
library(splatter)
library(tidyverse)
library(zinbwave)

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

#Method to clean just a vector of gene names
gene_vector_cleaner <- function(unclean.data=deseq2_rows){
  split_data <- strsplit(unclean.data, "_")
  finished_genes <- c()
  current_list <- split_data
  for (y in seq(1:length(current_list))){
    finished_genes <- c(finished_genes, current_list[[y]][2])
  }
  unclean.data <- finished_genes
  return(unclean.data)
  
  
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

a549_cells <- select(cell_lines, contains("A549"))
gm_cells <- select(h1_and_gm_cells, contains("GM"))

a549_colnames <- paste("C1",seq(1:length(colnames(a549_cells))),sep = ".") 
gm_colnames <- paste("C2",seq(1:length(colnames(gm_cells))),sep = ".")
colnames(a549_cells) <- a549_colnames
colnames(gm_cells) <- gm_colnames

unique_a549_names <- unique(rownames(a549_cells))
a549_cells <- a549_cells[unique_a549_names,]

unique_gm_names <- unique(rownames(gm_cells))
gm_cells <- gm_cells[unique_gm_names,]

condition <- c(rep(1, ncol(a549_cells)), rep(2, ncol(gm_cells)))

names(condition) <- c(colnames(a549_cells), colnames(gm_cells))


i <- c(1:ncol(a549_cells))  
z <- c(1:ncol(gm_cells))

a549_cells <- apply(a549_cells, c(1,2),           
                           function(x) as.numeric(as.character(x)))

gm_cells <- apply(gm_cells, c(1,2),            
                              function(x) as.numeric(as.character(x)))


sce <- SingleCellExperiment(assays=list(normcounts=cbind(a549_cells,
                                                         gm_cells)),
                            colData=data.frame(condition))

#Filtering the sce object
sce_filtered <- preprocess(sce, zero.thresh=0.9)




prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
sce_significance_test<- scDD(sce_filtered, prior_param=prior_param, testZeroes=FALSE, categorize = FALSE)
scdd_res_cell_line <- results(sce_significance_test)
scdd_res_cell_line <- scdd_res_cell_line[with(scdd_res_cell_line, order(nonzero.pvalue.adj)), ]
scdd_res_cell_line <- filter(scdd_res_cell_line, nonzero.pvalue.adj<0.05)



#Doing pre-processing to get my data ready to go into the SingleCellExperiment data for scDD method
#container format
tumor_colnames <- paste("C1",seq(1:length(colnames(all_tumor_cells_fpkm))),sep = ".")
nm_colnames <- paste("C2",seq(1:length(colnames(all_nm_cells_fpkm))),sep = ".")
colnames(all_tumor_cells_fpkm) <- tumor_colnames
colnames(all_nm_cells_fpkm) <- nm_colnames

condition <- c(rep(1, ncol(all_tumor_cells_fpkm)), rep(2, ncol(all_nm_cells_fpkm)))

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
scdd_res <- scdd_res[with(scdd_res, order(nonzero.pvalue.adj)), ]
scdd_res <- filter(scdd_res, nonzero.pvalue.adj<0.05)

save(scdd_res,file = "Data/Exported-data/R-objects/scdd_res.RData")
write.csv(scdd_res, file = "Data/Exported-data/Csv-files/scdd_res.csv")


#For Glio dataset----
#Doing pre-processing to get my data ready to go into the SingleCellExperiment data for scDD method
#container format
tumor_colnames <- paste("C1",seq(1:293),sep = ".")
nm_colnames <- paste("C2",seq(1:2),sep = ".")
colnames(glio_dg)[1:293] <- tumor_colnames
colnames(glio_dg)[294:295] <- nm_colnames

condition <- c(rep(1, 293), rep(2, 2))

names(condition) <- c(colnames(glio_dg[1:293]), colnames(glio_dg[294:295]))
i <- c(1:293)  
z <- c(1:2)

all_nm_cells_fpkm <- apply(glio_dg[294:295], c(1,2),           
                           function(x) as.numeric(as.character(x)))

all_tumor_cells_fpkm <- apply(glio_dg[1:293], c(1,2),            
                              function(x) as.numeric(as.character(x)))

#Making the SingleCellExperiment object of my data---- 
unique_tumor_names <- unique(rownames(all_tumor_cells_fpkm))
all_tumor_cells_fpkm <- all_tumor_cells_fpkm[unique_tumor_names,]

unique_nm_names <- unique(rownames(all_nm_cells_fpkm))
all_nm_cells_fpkm <- all_nm_cells_fpkm[unique_nm_names,]

condition <- c(rep(1, ncol(all_tumor_cells_fpkm)), rep(2, 2))


sce <- SingleCellExperiment(assays=list(normcounts=cbind(all_tumor_cells_fpkm,
                                                         all_nm_cells_fpkm)),
                            colData=data.frame(condition))

#Filtering the sce object
sce_filtered <- preprocess(sce, zero.thresh=0.9)


#For multiple cores----
BiocParallel::register(BiocParallel::MulticoreParam())

#scDD function and results----
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
sce_significance_test<- scDD(sce_filtered, prior_param=prior_param, testZeroes=FALSE, categorize = FALSE, min.size = 3)
scdd_res <- results(sce_significance_test)
scdd_res <- scdd_res[with(scdd_res, order(nonzero.pvalue.adj)), ]
scdd_res <- filter(scdd_res, nonzero.pvalue.adj<0.05)

save(scdd_res,file = "Data/Exported-data/R-objects/scdd_res_glio.RData")
write.csv(scdd_res, file = "Data/Exported-data/Csv-files/scdd_res_glio.csv")


#DESeq2----
#Reading in the count files
lung_cells_ds2 <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_tumor_all_cells_COUNT.csv")
nm_cells_ds2 <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_NM_all_cells_COUNT.csv")
normal_num <- rep("normal", 266)
tumor_num <- rep("tumor", 375)
all_nums <- c(tumor_num, normal_num)
# nm_names <- rep(paste("NM",1:ncol(nm_cells_ds2)))
# tumor_names <- rep(paste("TM",1:ncol(tumor_cells_ds2)))
# all_names <- c(tumor_names, nm_names)

#Merging them together into one large dataframe
deseq2_df <- merge(tumor_cells_ds2, nm_cells_ds2)
rownames(deseq2_df) <- deseq2_df$X
deseq2_df <- subset(deseq2_df, select=c(RHC3546__Tcell__.C6E879:RHC6187__Macrophage__.FFFF55))
deseq2_rows <- rownames(deseq2_df)
deseq2_rows <- gene_vector_cleaner(unclean.data = deseq2_rows)
deseq2_cols <- colnames(deseq2_df)
deseq2_df <- deseq2_df[complete.cases(deseq2_df), ]
deseq2_df <- apply(deseq2_df, c(1,2), as.integer)
#rownames(deseq2_df) <- seq(1, nrow(deseq2_df), by=1)
#colnames(deseq2_df) <- seq(1, ncol(deseq2_df), by=1)
deseq2_mat <- as.matrix(deseq2_df)

#Making a summarizedExperiment object from the combined dataframe
deseq2_se <- SummarizedExperiment(assays=list(counts=deseq2_mat),
                                  colData=DataFrame(label=deseq2_cols),
                                  rowData=DataFrame(length=deseq2_rows),)



keep <-  rowSums(assay(deseq2_se)>5)>10
table(keep)

zinb_data <- deseq2_se[keep,]
zinb_data$cell_state <- all_nums


zinb_data <- zinb_data[names(zinb_data)[1:5000],]

#For multicore
BiocParallel::register(BiocParallel::MulticoreParam())

zinb_data <- zinbwave(zinb_data, K=0, BPPARAM = BiocParallel::bpparam(), epsilon=1e12, normalizedValues=FALSE, observationalWeights = TRUE, verbose=TRUE)

dds <- DESeqDataSet(zinb_data, design=~cell_state)

dds <- estimateSizeFactors(dds, type="poscounts")

scr <- scran::calculateSumFactors(dds)
sizeFactors(dds) <- scr

dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6, minRep=Inf) 


dds_res <- DESeq2::results(dds)

resOrdered <- dds_res[order(dds_res$padj),]
save(resOrdered, file = "Data/Data-from-Cleaner-code/deseq2_top5000_genes_cc_patients.RData")

resOrdered_subset <- resOrdered[1:1800,]
save(resOrdered_subset, file = "Data/Data-from-Cleaner-code/deseq2_top1800_genes_cc_patients.RData")

cleaned_gene_names <- gene_vector_cleaner(unclean.data = rownames(resOrdered_subset))
resOrdered_subset_finished <- data.frame(gene=cleaned_gene_names, baseMean=resOrdered_subset$baseMean, log2FoldChange=resOrdered_subset$log2FoldChange, lfcSE=resOrdered_subset$lfcSE, stat=resOrdered_subset$stat, pvalue=resOrdered_subset$pvalue, padj=resOrdered_subset$padj)
write.csv(resOrdered_subset_finished, file = "Data/Data-from-Cleaner-code/deseq2_top1800_genes_cc_patients.csv")


#For TCGA-LUAD
#Reading in the count files
tumor_cells_ds2 <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_tumor_all_cells_COUNT.csv")
nm_cells_ds2 <- read.csv(file = "Data/Single-cell-data/Counts/GSE81861_CRC_NM_all_cells_COUNT.csv")
normal_num <- rep("normal", 266)
tumor_num <- rep("tumor", 375)
all_nums <- c(tumor_num, normal_num)
# nm_names <- rep(paste("NM",1:ncol(nm_cells_ds2)))
# tumor_names <- rep(paste("TM",1:ncol(tumor_cells_ds2)))
# all_names <- c(tumor_names, nm_names)

#Merging them together into one large dataframe
deseq2_df <- merge(tumor_cells_ds2, nm_cells_ds2)
rownames(deseq2_df) <- deseq2_df$X
deseq2_df <- subset(deseq2_df, select=c(RHC3546__Tcell__.C6E879:RHC6187__Macrophage__.FFFF55))
deseq2_rows <- rownames(deseq2_df)
deseq2_rows <- gene_vector_cleaner(unclean.data = deseq2_rows)
deseq2_cols <- colnames(deseq2_df)
deseq2_df <- deseq2_df[complete.cases(deseq2_df), ]
deseq2_df <- apply(deseq2_df, c(1,2), as.integer)
#rownames(deseq2_df) <- seq(1, nrow(deseq2_df), by=1)
#colnames(deseq2_df) <- seq(1, ncol(deseq2_df), by=1)
deseq2_mat <- as.matrix(deseq2_df)

#Making a summarizedExperiment object from the combined dataframe
deseq2_se <- SummarizedExperiment(assays=list(counts=deseq2_mat),
                                  colData=DataFrame(label=deseq2_cols),
                                  rowData=DataFrame(length=deseq2_rows),)



keep <-  rowSums(assay(deseq2_se)>5)>10
table(keep)

zinb_data <- deseq2_se[keep,]
zinb_data$cell_state <- all_nums


zinb_data <- zinb_data[names(zinb_data)[1:5000],]

#For multicore
BiocParallel::register(BiocParallel::MulticoreParam())

zinb_data <- zinbwave(zinb_data, K=0, BPPARAM = BiocParallel::bpparam(), epsilon=1e12, normalizedValues=FALSE, observationalWeights = TRUE, verbose=TRUE)

dds <- DESeqDataSet(zinb_data, design=~cell_state)

dds <- estimateSizeFactors(dds, type="poscounts")

scr <- scran::calculateSumFactors(dds)
sizeFactors(dds) <- scr

dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6, minRep=Inf) 


dds_res <- DESeq2::results(dds)

resOrdered <- dds_res[order(dds_res$padj),]
save(resOrdered, file = "Data/Data-from-Cleaner-code/deseq2_top5000_genes_cc_patients.RData")

resOrdered_subset <- resOrdered[1:1800,]
save(resOrdered_subset, file = "Data/Data-from-Cleaner-code/deseq2_top1800_genes_cc_patients.RData")

cleaned_gene_names <- gene_vector_cleaner(unclean.data = rownames(resOrdered_subset))
resOrdered_subset_finished <- data.frame(gene=cleaned_gene_names, baseMean=resOrdered_subset$baseMean, log2FoldChange=resOrdered_subset$log2FoldChange, lfcSE=resOrdered_subset$lfcSE, stat=resOrdered_subset$stat, pvalue=resOrdered_subset$pvalue, padj=resOrdered_subset$padj)
write.csv(resOrdered_subset_finished, file = "Data/Data-from-Cleaner-code/deseq2_top1800_genes_cc_patients.csv")



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
param <- MulticoreParam(workers = 4, progressbar = TRUE)
register(param)


des_results <- DEsingle(counts = sce, group = condition, parallel = TRUE, BPPARAM = param)
save(des_results, file = "Data/Data-from-Cleaner-code/des_results_cc_patients.RData")

#load("Data/Exported-data/R-objects/des_results.RData", verbose = TRUE)



#edgeR----
dge <- DGEList(assay(zinb_data))
dge <- calcNormFactors(dge)

weights <- assay(zinb_data, "weights")


design <- model.matrix(~cell_state, data = colData(zinb_data))
dge$weights <- weights
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

lrt <- glmWeightedF(fit, coef = 1:2)
topTags(lrt)

finished_edger <- head(lrt$table, n=1800)
finished_edger <- finished_edger[order(finished_edger$padjFilter), ]
write.csv(finished_edger, file = "Data/Data-from-Cleaner-code/finished_edgeR_genes.csv")
clean_edger_names <- gene_vector_cleaner(edger_df_done$X)





#EMDomics----
#Reading in the files
all_tumor_cells_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_tumor_all_cells_FPKM.csv")
all_nm_cells_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_NM_all_cells_FPKM.csv")

rownames(all_tumor_cells_fpkm) <- all_tumor_cells_fpkm$X
all_tumor_cells_fpkm <- gene_name_cleaner(data.to.clean = all_tumor_cells_fpkm)
all_tumor_cells_fpkm <- t(all_tumor_cells_fpkm)
all_tumor_cells_fpkm <- subset(all_tumor_cells_fpkm, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
save(all_tumor_cells_fpkm, file = "Data/Exported-data/R-objects/all_tumor_cells_fpkm_for_emdomics.RData")
write.csv(all_tumor_cells_fpkm, file = "Data/Exported-data/Csv-files/all_tumor_cells_fpkm_for_emdomics.csv")

all_nm_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_NM_all_cells_FPKM.csv")
rownames(all_nm_cells_fpkm) <- all_nm_cells_fpkm$X
all_nm_cells_fpkm <- gene_name_cleaner(data.to.clean = all_nm_cells_fpkm)
all_nm_cells_fpkm <- t(all_nm_cells_fpkm)
all_nm_cells_fpkm <- subset(all_nm_cells_fpkm, select=c(RHC3934__Bcell__.7DEA7B:RHC6187__Macrophage__.FFFF55))
save(all_nm_cells_fpkm, file = "Data/Exported-data/R-objects/all_nm_cells_fpkm_for_emdomics.RData")
write.csv(all_nm_cells_fpkm, file = "Data/Exported-data/Csv-files/all_nm_cells_fpkm_for_emdomics.csv")

unique_tumor_names <- unique(rownames(all_tumor_cells_fpkm))
all_tumor_cells_fpkm <- all_tumor_cells_fpkm[unique_tumor_names,]
save(all_tumor_cells_fpkm, file = "Data/Exported-data/R-objects/all_tumor_cells_fpkm_for_emdomics_formatted.RData")
write.csv(all_tumor_cells_fpkm, file = "Data/Exported-data/Csv-files/all_tumor_cells_fpkm_for_emdomics_formatted.csv")



unique_nm_names <- unique(rownames(all_nm_cells_fpkm))
all_nm_cells_fpkm <- all_nm_cells_fpkm[unique_nm_names,]
save(all_nm_cells_fpkm, file = "Data/Exported-data/R-objects/all_nm_cells_fpkm_for_emdomics_formatted.RData")
write.csv(all_nm_cells_fpkm, file = "Data/Exported-data/Csv-files/all_nm_cells_fpkm_for_emdomics_formatted.csv")


condition <- c(rep("A", ncol(all_tumor_cells_fpkm)), rep("B", ncol(all_nm_cells_fpkm)))

all_crc_data_together <- merge(all_tumor_cells_fpkm, all_nm_cells_fpkm, by=0)
rownames(all_crc_data_together) <- all_crc_data_together$Row.names
all_crc_data_together <- subset(all_crc_data_together, select=c(RHC3546__Tcell__.C6E879:RHC6187__Macrophage__.FFFF55))

names(condition) <- colnames(all_crc_data_together)


all_crc_data_together <- as.matrix(all_crc_data_together)

counts <- all_crc_data_together

#Seurat filtering
seurat<-CreateSeuratObject(raw.data = counts, min.cells = 3, min.genes = 350, project = "My_project")



#For multi-core support
BiocParallel::register(BiocParallel::MulticoreParam())

#Doing the actual calculation
emdomics_results <- calculate_emd(all_crc_data_together, condition, nperm=100, parallel=TRUE)

#scde----

#SINCERA----



#Plotting of the different metrics----
#Image loading function----
my_img_loader <- function(use.ebimage=TRUE,
                          img.to.load=NULL,
                          show.img   =TRUE, 
                          my.method  ="raster",
                          png.info   =TRUE){
  require(EBImage)
  require(png)
  require(grid)
  require(gridExtra)
  
  
  if(use.ebimage==TRUE){
    loaded_img <- EBImage::readImage(img.to.load)
    
    if(show.img==TRUE){
      display(loaded_img, method = my.method)
    }
    
  }else{
    loaded_img <- png::readPNG(source = img.to.load, info = png.info)
    
  }
  
  
  return(loaded_img)
}

#Using it to load images----
setwd("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/")
cc_names <- list.files("Completed_KM_Curves/CC_Singlecell/")
deseq2_names <- list.files("Completed_KM_Curves/DESeq2/")
desingle_names <- list.files("Completed_KM_Curves/Desingle/")
edger_names <- list.files("Completed_KM_Curves/edgeR/")
scdd_names <- list.files("Completed_KM_Curves/scDD/")
cc_genes <- cc_names[cc_names %in% c("mad_sdes_cindex_combo_that_seperates_rc_patients.png", "mirna_sdes_cc_patients_cindex_that_seperates.png", "mirna_sdes_top_cindex_combo_rc_patients.png", "mirna_sdes_top_cindex_combo_cc_and_rc_patients.png", "mad_sdes_top_cindex_combo_cc_patients.png", "mad_sdes_top_cindex_combo_cc_and_rc_patients.png")]
deseq2_genes <- deseq2_names[deseq2_names %in% c("cc_patients.png","rc_patients.png", "cc_and_rc_patients.png")]
desingle_genes <- desingle_names
edger_genes <- edger_names[edger_names %in% c("cc_patients.png", "rc_patients.png", "cc_and_rc_patients.png")]
scdd_genes <- scdd_names[scdd_names %in% c("genes_only_cc_patients.png", "genes_only_rc_patients.png", "genes_only_cc_and_rc_patients.png")]

all_genes <- c(cc_genes, deseq2_genes, edger_genes, desingle_genes, scdd_genes)

cindex_boxplot <- readPNG("Completed_KM_Curves/genes_only_even_more_progress.png")
cindex_boxplot <- rasterGrob(cindex_boxplot)

#Looping through and loading them----
setwd("Completed_KM_Curves/scDD/")
#my_imgs <- list()
#my_deseq2_imgs <- list()
#my_edger_imgs <- list()
#my_desingle_imgs <- list()
my_scdd_imgs <- list()
counter <- 1
for (f in scdd_genes) {
  current_img <- my_img_loader(img.to.load = f, show.img = FALSE, use.ebimage = FALSE, png.info = TRUE)
  my_scdd_imgs[[as.character(counter)]] <- current_img
  counter <- counter + 1
}


#Putting all the images together----
km_plots <- grid.arrange(rasterGrob(my_imgs$`4`),rasterGrob(my_imgs$`6`),rasterGrob(my_imgs$`3`), rasterGrob(my_imgs$`1`), rasterGrob(my_imgs$`5`), rasterGrob(my_imgs$`2`), rasterGrob(my_deseq2_imgs$`2`), rasterGrob(my_deseq2_imgs$`3`), rasterGrob(my_deseq2_imgs$`1`), rasterGrob(my_edger_imgs$`2`), rasterGrob(my_edger_imgs$`3`), rasterGrob(my_edger_imgs$`1`),rasterGrob(my_desingle_imgs$`2`),rasterGrob(my_desingle_imgs$`3`), rasterGrob(my_desingle_imgs$`1`),rasterGrob(my_scdd_imgs$`2`), rasterGrob(my_scdd_imgs$`3`), rasterGrob(my_scdd_imgs$`1`), ncol=3, nrow=6)
plot_grid(cindex_boxplot, km_plots, rel_heights = c(1/3,2/3),nrow = 2, ncol = 1, labels = c("A.", "B."), hjust = -8.0)
