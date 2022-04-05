#Name: reproducible.R
#Purpose: For reproducing all figures & results in the manuscript
#Author: Andrew Willems <awillems@vols.utk.edu>

#Loading needed packages----
library(boot)
library(data.table)
library(ggplot2)
library(phateR)
library(selectiveInference)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(svglite)
library(TCGAbiolinks)
library(tidyverse)
library(viridis)

#Loading the needed functions from their respective files----
source("Code/cell_dataset_builder.R")
source("Code/cox_model.R")
source("Code/mad_calculator.R")
source("Code/magic_denoiser.R")
source("Code/mirna_calculator_original.R")
source("Code/model_optimizer.R")
source("Code/switchde_calculator.R")

#Loading single-cell data that is used for both colon and rectal cancer data sets----
cc_tumor_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_tumor_all_cells_FPKM.csv")

#Processing the single-cell data to clean up the gene names----
rownames(cc_tumor_fpkm) <- cc_tumor_fpkm$X
current_colname_split <- strsplit(colnames(t(cc_tumor_fpkm)), "_")
finished_gene_list <- c()
for (x in seq(1:length(current_colname_split))){
  finished_gene_list <- c(finished_gene_list, current_colname_split[[x]][2])
}

cc_tumor_fpkm <- t(cc_tumor_fpkm)
colnames(cc_tumor_fpkm) <- finished_gene_list
cc_tumor_fpkm <- t(cc_tumor_fpkm)
cc_tumor_fpkm <- subset(cc_tumor_fpkm,
                        select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
cc_tumor_fpkm <- apply(cc_tumor_fpkm, c(1,2), as.numeric)


#First pre-processing the scRNA-seq data before sending it to MAGIC
#Keeping genes expressed in at least 10 cells
keep_rows <- rowSums(cc_tumor_fpkm > 0) > 10
cc_tumor_fpkm <- cc_tumor_fpkm[keep_rows,]

#Looking at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=colSums(cc_tumor_fpkm)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

#Normalizing the library size
cc_tumor_fpkm <- library.size.normalize(cc_tumor_fpkm)
cc_tumor_fpkm <- sqrt(cc_tumor_fpkm)


#Denoising the single-cell data and saving the output----
cc_tumor_fpkm <- magic_denoiser(sc.data = cc_tumor_fpkm,
                                magic.seed = 123, magic.solver = 'approximate')
saveRDS(cc_tumor_fpkm,
        file = "Data/Reproducible-results/denoised-colon-and-rectal-single-cell-data.rds")


#Generating the VIM based pseudotime progression with Monocle3----
#Based on the expression of the VIM graph we select the bottom-rightmost point
#as the root of our pseudotime as this matches the gradient of VIM expression
#level
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP"),
                                   cell.data = cc_tumor_fpkm$denoised_sc_dataframe,
                                   cell.meta = cc_tumor_fpkm$cds_gene_names,
                                   my.root = "Y_19",
                                   my.cds.filename = "Data/Reproducible-results/Data/cds_ouput.rds",
                                   my.monocle.graph = "Data/Reproducible-results/Data/vim_pseudotime.svg",
                                   my.monocle.graph.genes = "Data/Reproducible-results/Data/vim_genes.svg",
                                   my.pt.data.filename = "Data/Reproducible-results/Data/vim_pseudotime_data.csv")


#MAD metric for colon and rectal cancer----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
saveRDS(mad.genes,
        file = "Data/Reproducible-results/Data/mad_colon_and_rectal_cancer.rds")

#Switchde metric for colon and rectal cancer----
sde.genes <- switchde_calculator(cc_tumor_fpkm$denoised_sc_dataframe,
                                 pseudo.time = cds_output$Pseudotime)
saveRDS(sde.genes,
        file = "Data/Reproducible-results/Data/sde_colon_and_rectal_cancer.rds")

#Now getting our Colon cancer bulk data set----
#TCGA-COAD
# coad_query <- GDCquery(project       = "TCGA-COAD",
#                        data.category = "Transcriptome Profiling",
#                        data.type     = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - FPKM")
# 
# #Downloading the files
# GDCdownload(query           = coad_query,
#             method          = "api",
#             files.per.chunk = 10,
#             directory       = "Data/Reproducible-results/Data/Bulk-data/")
# 
# 
# #Making the SummarizedExperiment object
# coad_data_se <- GDCprepare(coad_query, summarizedExperiment = TRUE,
#                            directory = "Data/Reproducible-results/Data/Bulk-data/")
# coad_data_df <- as.data.frame(colData(coad_data_se))
# coad_data_df$vital_status <- factor(coad_data_df$vital_status,
#                                     levels = c("Alive", "Dead"),
#                                     labels = c(0,1))
# coad_data_df$vital_status <- as.numeric(as.character(coad_data_df$vital_status))
# 
# 
# bulk_rna_df <- coad_data_se@assays@data@listData[["HTSeq - FPKM"]]
# colnames(bulk_rna_df) <- coad_data_se@colData@rownames
# rownames(bulk_rna_df) <- coad_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
# bulk_rna_df <- t(bulk_rna_df)
# bulk_rna_df <- as.data.frame(bulk_rna_df)
# bulk_rownames <- rownames(bulk_rna_df)
# bulk_rna_df$barcode <- bulk_rownames
# 
# bulk_rna_df_unique <- subset(bulk_rna_df,
#                              select = unique(colnames(bulk_rna_df)))
# coad_data_df_unique <- subset(coad_data_df,
#                               select = unique(colnames(coad_data_df)))
# merged_df <- merge(bulk_rna_df_unique, coad_data_df_unique, by = 'barcode')
# rownames(merged_df) <- merged_df$barcode
# merged_df <- merged_df[,2:length(colnames(merged_df))]
# 
# 
# 
# merged_df$days_to_last_follow_up <- ifelse(merged_df$vital_status==1,
#                                            merged_df$days_to_death, 
#                                            merged_df$days_to_last_follow_up)
# 
# merged_df <- filter(merged_df, days_to_last_follow_up != "NA")
# 
# 
# cox_time <- merged_df$days_to_last_follow_up
# cox_event <- merged_df$vital_status
# cox_tumor <- merged_df$ajcc_pathologic_stage
# cox_tumor_n <- merged_df$ajcc_pathologic_n
# cox_tumor_m <- merged_df$ajcc_pathologic_m
# cox_gender <- merged_df$gender
# cox_eth <- merged_df$ethnicity
# cox_race <- merged_df$race
# cox_type <- merged_df$definition
# cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
# cox_df$days.to.last.follow.up <- cox_time
# cox_df$vital.status <- cox_event
# cox_df$tumor.stage <- cox_tumor
# cox_df$ajcc.m <- cox_tumor_m
# cox_df$ajcc.n <- cox_tumor_n
# cox_df$race <- cox_race
# cox_df$ethnicity <- cox_eth
# cox_df$gender <- cox_gender
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage IV", replacement = 4)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage III", replacement = 3)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage II", replacement = 2)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage I", replacement = 1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
# cox_df$sample.type <- cox_type
# cox_df <- filter(cox_df, !tumor.stage=="not reported")
# cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]
# cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up < 1, 1,
#                                         cox_df$days.to.last.follow.up)
# saveRDS(cox_df, "Data/TCGA-COAD/coad_df_finished.rds")

cox_df <- readRDS("Data/Reproducible-results/Data/coad_df_finished.rds")

#Getting ideal gene number for MAD metric on TCGA-COAD----
gene_sizes <- seq(100, 3000, 50)
mad_cindices <- rep(0, length(gene_sizes))

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = mad.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/mad_coad_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mad_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mad_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mad_cindices_coad_df <- as.data.frame(cbind(gene_sizes, mad_cindices))
mad_cindices_coad_df$method <- rep("MAD", nrow(mad_cindices_coad_df))
mad_cindices_coad_df <- mad_cindices_coad_df[,2:4]
colnames(mad_cindices_coad_df)[2] <- "c_index"
write.csv(mad_cindices_coad_df,
          "Data/Reproducible-results/Data/mad_cindices_coad_across_gene_size.csv")




#Getting ideal gene number for SDE metric on TCGA-COAD----
gene_sizes <- seq(100, 3000, 50)
sde_cindices <- rep(0, length(gene_sizes))

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = sde.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/sde_coad_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the sde_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  sde_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
sde_cindices_coad_df <- as.data.frame(cbind(gene_sizes, sde_cindices))
sde_cindices_coad_df$method <- rep("SDE", nrow(sde_cindices_coad_df))
sde_cindices_coad_df <- sde_cindices_coad_df[,2:4]
colnames(sde_cindices_coad_df)[2] <- "c_index"
write.csv(sde_cindices_coad_df,
          "Data/Reproducible-results/Data/sde_cindices_coad_across_gene_size.csv")


#Getting ideal gene number for miRNA metric on TCGA-COAD----
#For this metric we don't know which combination of miRNA and miRNA targets
#will yield the best result in advance so we are trying 3 different miRNA-
#miRNA target combinations that encompass specific areas of our grid search.
#They are low miRNA number and low miRNA target number, medium miRNA number and
#medium miRNA target number, and high miRNA number and high miRNA target number.
#Once the ideal miRNA-miRNA target pair is known from the grid search we will
#also include a gene size search for it here
gene_sizes <- seq(100, 3000, 50)
mirna_high_cindices <- rep(0, length(gene_sizes))
mirna_medium_cindices <- rep(0, length(gene_sizes))
mirna_low_cindices <- rep(0, length(gene_sizes))

#Loading the high miRNA-miRNA target file
load("Data/Reproducible-results/Data/800_1010_targets.RData", verbose = TRUE)

high.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = high.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/high_mirna_coad_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_high_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_high_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_high_cindices_coad_df <- as.data.frame(cbind(gene_sizes,
                                                   mirna_high_cindices))

colnames(mirna_high_cindices_coad_df) [2] <- "c_index"
mirna_high_cindices_coad_df$mirna_type <- rep("high",
                                              nrow(mirna_high_cindices_coad_df))

write.csv(mirna_high_cindices_coad_df,
          "Data/Reproducible-results/Data/mirna_high_cindices_coad_across_gene_size.csv")


#Medium miRNA-miRNA target number
#Loading the medium miRNA-miRNA target file
load("Data/Reproducible-results/Data/400_510_targets.RData", verbose = TRUE)

medium.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = medium.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/medium_mirna_coad_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_medium_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_medium_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_medium_cindices_coad_df <- as.data.frame(cbind(gene_sizes,
                                                   mirna_medium_cindices))

colnames(mirna_medium_cindices_coad_df) [2] <- "c_index"
mirna_medium_cindices_coad_df$mirna_type <- rep("medium",
                                                nrow(mirna_medium_cindices_coad_df))

write.csv(mirna_medium_cindices_coad_df,
          "Data/Reproducible-results/Data/mirna_medium_cindices_coad_across_gene_size.csv")

#low miRNA-miRNA target number
#Loading the low miRNA-miRNA target file
load("Data/Reproducible-results/Data/100_10_targets.RData", verbose = TRUE)

low.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = low.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/low_mirna_coad_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_medium_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_low_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_low_cindices_coad_df <- as.data.frame(cbind(gene_sizes,
                                                     mirna_low_cindices))

colnames(mirna_low_cindices_coad_df)[2] <- "c_index"
mirna_low_cindices_coad_df$mirna_type <- rep("low",
                                             nrow(mirna_low_cindices_coad_df))

write.csv(mirna_low_cindices_coad_df,
          "Data/Reproducible-results/Data/mirna_low_cindices_coad_across_gene_size.csv")


#Binding all the rows of the 3 different miRNA data frames together into 1 big
#data frame
all_mirna_metrics_coad_df <- bind_rows(mirna_low_cindices_coad_df,
                                            mirna_medium_cindices_coad_df,
                                            mirna_high_cindices_coad_df)

all_mirna_metrics_coad_df <- all_mirna_metrics_coad_df[,2:4]

#Reordering the labels to make them look nicer on the plot
all_mirna_metrics_coad_df$mirna_type <- factor(all_mirna_metrics_coad_df$mirna_type,
                                               levels = c("low", "medium", "high"))

#Basic miRNA-miRNA target plot
all_mirna_metrics_coad_plot <- ggplot(data = all_mirna_metrics_coad_df,
                                      aes(x=gene_sizes, y=c_index, color=mirna_type))

#Making the plot much nicer
all_mirna_metrics_coad_plot_finished <-all_mirna_metrics_coad_plot +
  geom_line(size=2.5) + geom_point(size=3.0)+
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", 
        legend.text = element_text(size = 35, face = "plain"),
        axis.text=element_text(size=25, face="bold"))+
  xlab("Gene Size")+ ylab("Mean Concordance Index")+
  ggtitle("miRNA-miRNA Target Number")+ 
  labs(color="miRNA-miRNA Target Number")+
  scale_color_viridis_d()


#Saving the finished graph in .svg format
ggsave(filename = "Data/Reproducible-results/Figures/mirna_mirna_target_num_across_gene_size.svg",
       plot     = print(all_mirna_metrics_coad_plot_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")


#Optimal gene size for MAD, SDE, and miRNA high metric----
#First bind the data frames together
methods_optimal_gene_coad_df <- bind_rows(mad_cindices_coad_df, 
                                        sde_cindices_coad_df,
                                        mirna_high_cindices_coad_df)

methods_optimal_gene_coad_df <- methods_optimal_gene_coad_df[,1:3]

methods_optimal_gene_coad_df$method[119:177] <- rep("miRNA", 58)

#Saving the data frame
write.csv(methods_optimal_gene_coad_df,
          "Data/Reproducible-results/Data/individual_optimal_gene_point_coad_data.csv")

#Basic miRNA-miRNA target plot
methods_optimal_gene_coad_plot <- ggplot(data = methods_optimal_gene_coad_df,
                                      aes(x=gene_sizes, y=c_index, color=method))

#Making the plot much nicer
methods_optimal_gene_coad_plot_finished <-methods_optimal_gene_coad_plot +
  geom_line(size=2.5) + geom_point(size=3.0)+
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", 
        legend.text = element_text(size = 35, face = "plain"),
        axis.text=element_text(size=25, face="bold"))+
  xlab("Gene Size")+ ylab("Mean Concordance Index")+
  ggtitle("Optimal Gene Number")+ 
  labs(color="Method")+
  scale_color_manual(values = c("#21908CFF","#FDE725FF", "#440154FF"))


#Saving the finished graph in .svg format
ggsave(filename = "Data/Reproducible-results/Figures/individual_metrics_optimal_gene_size_coad.svg",
       plot     = print(methods_optimal_gene_coad_plot_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")




#Now doing all of these steps for the rectal cancer data set from TCGA-READ
#Now getting our rectal cancer bulk data set----
#TCGA-READ
# read_query <- GDCquery(project       = "TCGA-READ",
#                        data.category = "Transcriptome Profiling",
#                        data.type     = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - FPKM")
# 
# #Downloading the files
# GDCdownload(query           = read_query,
#             method          = "api",
#             files.per.chunk = 10,
#             directory       = "Data/Reproducible-results/Data/Bulk-data/")
# 

#Making the SummarizedExperiment object
# read_data_se <- GDCprepare(read_query, summarizedExperiment = TRUE,
#                            directory = "Data/Reproducible-results/Data/Bulk-data/")
# read_data_df <- as.data.frame(colData(read_data_se))
# read_data_df$vital_status <- factor(read_data_df$vital_status,
#                                     levels = c("Alive", "Dead"),
#                                     labels = c(0,1))
# read_data_df$vital_status <- as.numeric(as.character(read_data_df$vital_status))
# 
# 
# bulk_rna_df <- read_data_se@assays@data@listData[["HTSeq - FPKM"]]
# colnames(bulk_rna_df) <- read_data_se@colData@rownames
# rownames(bulk_rna_df) <- read_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
# bulk_rna_df <- t(bulk_rna_df)
# bulk_rna_df <- as.data.frame(bulk_rna_df)
# bulk_rownames <- rownames(bulk_rna_df)
# bulk_rna_df$barcode <- bulk_rownames
# 
# bulk_rna_df_unique <- subset(bulk_rna_df,
#                              select = unique(colnames(bulk_rna_df)))
# read_data_df_unique <- subset(read_data_df,
#                               select = unique(colnames(read_data_df)))
# merged_df <- merge(bulk_rna_df_unique, read_data_df_unique, by = 'barcode')
# rownames(merged_df) <- merged_df$barcode
# merged_df <- merged_df[,2:length(colnames(merged_df))]
# 
# 
# 
# merged_df$days_to_last_follow_up <- ifelse(merged_df$vital_status==1,
#                                            merged_df$days_to_death,
#                                            merged_df$days_to_last_follow_up)
# 
# merged_df <- filter(merged_df, days_to_last_follow_up != "NA")
# 
# 
# cox_time <- merged_df$days_to_last_follow_up
# cox_event <- merged_df$vital_status
# cox_tumor <- merged_df$ajcc_pathologic_stage
# cox_tumor_n <- merged_df$ajcc_pathologic_n
# cox_tumor_m <- merged_df$ajcc_pathologic_m
# cox_gender <- merged_df$gender
# cox_eth <- merged_df$ethnicity
# cox_race <- merged_df$race
# cox_type <- merged_df$definition
# cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
# cox_df$days.to.last.follow.up <- cox_time
# cox_df$vital.status <- cox_event
# cox_df$tumor.stage <- cox_tumor
# cox_df$ajcc.m <- cox_tumor_m
# cox_df$ajcc.n <- cox_tumor_n
# cox_df$race <- cox_race
# cox_df$ethnicity <- cox_eth
# cox_df$gender <- cox_gender
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage IV", replacement = 4)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage III", replacement = 3)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage II", replacement = 2)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage I", replacement = 1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
# cox_df$sample.type <- cox_type
# cox_df <- filter(cox_df, !tumor.stage=="not reported")
# cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]
# cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up < 1, 1,
#                                         cox_df$days.to.last.follow.up)
# saveRDS(cox_df, "Data/TCGA-READ/read_df_finished.rds")

cox_df <- readRDS("Data/Reproducible-results/Data/read_df_finished.rds")

#Getting ideal gene number for MAD metric on TCGA-READ----
gene_sizes <- seq(100, 3000, 50)
mad_cindices <- rep(0, length(gene_sizes))

#Loading in the MAD and SDE files if they aren't already loaded into the 
#environment from earlier
mad.genes <-readRDS("Data/Reproducible-results/Data/mad_colon_and_rectal_cancer.rds")
sde.genes <- readRDS("Data/Reproducible-results/Data/sde_colon_and_rectal_cancer.rds")

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = mad.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/mad_read_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mad_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mad_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mad_cindices_read_df <- as.data.frame(cbind(gene_sizes, mad_cindices))
mad_cindices_read_df$method <- rep("MAD", nrow(mad_cindices_read_df))
colnames(mad_cindices_read_df)[2] <- "c_index"
write.csv(mad_cindices_read_df,
          "Data/Reproducible-results/Data/mad_cindices_read_across_gene_size.csv")




#Getting ideal gene number for SDE metric on TCGA-READ----
gene_sizes <- seq(100, 3000, 50)
sde_cindices <- rep(0, length(gene_sizes))

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = sde.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/sde_read_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the sde_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  sde_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
sde_cindices_read_df <- as.data.frame(cbind(gene_sizes, sde_cindices))
sde_cindices_read_df$method <- rep("SDE", nrow(sde_cindices_read_df))
colnames(sde_cindices_read_df)[2] <- "c_index"
write.csv(sde_cindices_read_df,
          "Data/Reproducible-results/Data/sde_cindices_read_across_gene_size.csv")


#Getting ideal gene number for miRNA metric on TCGA-read----
#For this metric we don't know which combination of miRNA and miRNA targets
#will yield the best result in advance so we are trying 3 different miRNA-
#miRNA target combinations that encompass specific areas of our grid search.
#They are low miRNA number and low miRNA target number, medium miRNA number and
#medium miRNA target number, and high miRNA number and high miRNA target number.
#Once the ideal miRNA-miRNA target pair is known from the grid search we will
#also include a gene size search for it here
gene_sizes <- seq(100, 3000, 50)
mirna_high_cindices <- rep(0, length(gene_sizes))
mirna_medium_cindices <- rep(0, length(gene_sizes))
mirna_low_cindices <- rep(0, length(gene_sizes))

#Loading the high miRNA-miRNA target file
load("Data/Reproducible-results/Data/800_1010_targets.RData", verbose = TRUE)
load("Data/Reproducible-results/Data/1000_110_targets.RData", verbose = TRUE)
load("Data/Reproducible-results/Data/1000_1010_targets.RData", verbose = TRUE)

high.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = high.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/high_mirna3_read_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_high_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_high_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_high_cindices_read_df <- as.data.frame(cbind(gene_sizes,
                                                   mirna_high_cindices))

colnames(mirna_high_cindices_read_df) [2] <- "c_index"
mirna_high_cindices_read_df$mirna_type <- rep("high",
                                              nrow(mirna_high_cindices_read_df))

write.csv(mirna_high_cindices_read_df,
          "Data/Reproducible-results/Data/mirna_high_cindices_read_across_gene_size3.csv")


#Medium miRNA-miRNA target number
#Loading the medium miRNA-miRNA target file
load("Data/Reproducible-results/Data/400_510_targets.RData", verbose = TRUE)

medium.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = medium.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/medium_mirna_read_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_medium_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_medium_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_medium_cindices_read_df <- as.data.frame(cbind(gene_sizes,
                                                     mirna_medium_cindices))

colnames(mirna_medium_cindices_read_df) [2] <- "c_index"
mirna_medium_cindices_read_df$mirna_type <- rep("medium",
                                                nrow(mirna_medium_cindices_read_df))

write.csv(mirna_medium_cindices_read_df,
          "Data/Reproducible-results/Data/mirna_medium_cindices_read_across_gene_size.csv")

#low miRNA-miRNA target number
#Loading the low miRNA-miRNA target file
load("Data/Reproducible-results/Data/100_10_targets.RData", verbose = TRUE)

low.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = low.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/low_mirna_read_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_medium_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_low_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_low_cindices_read_df <- as.data.frame(cbind(gene_sizes,
                                                  mirna_low_cindices))

colnames(mirna_low_cindices_read_df)[2] <- "c_index"
mirna_low_cindices_read_df$mirna_type <- rep("low",
                                             nrow(mirna_low_cindices_read_df))

write.csv(mirna_low_cindices_read_df,
          "Data/Reproducible-results/Data/mirna_low_cindices_read_across_gene_size.csv")


#Binding all the rows of the 3 different miRNA data frames together into 1 big
#data frame
all_mirna_metrics_read_df <- bind_rows(mirna_low_cindices_read_df,
                                       mirna_medium_cindices_read_df,
                                       mirna_high_cindices_read_df)

all_mirna_metrics_read_df <- all_mirna_metrics_read_df[,2:4]

#Reordering the labels to make them look nicer on the plot
all_mirna_metrics_read_df$mirna_type <- factor(all_mirna_metrics_read_df$mirna_type,
                                               levels = c("low", "medium", "high"))

#Basic miRNA-miRNA target plot
all_mirna_metrics_read_plot <- ggplot(data = all_mirna_metrics_read_df,
                                      aes(x=gene_sizes, y=c_index, color=mirna_type))

#Making the plot much nicer
all_mirna_metrics_read_plot_finished <-all_mirna_metrics_read_plot +
  geom_line(size=2.5) + geom_point(size=3.0)+
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", 
        legend.text = element_text(size = 35, face = "plain"),
        axis.text=element_text(size=25, face="bold"))+
  xlab("Gene Size")+ ylab("Mean Concordance Index")+
  ggtitle("miRNA-miRNA Target Number")+ 
  labs(color="miRNA-miRNA Target Number")+
  scale_color_viridis_d()


#Saving the finished graph in .svg format
ggsave(filename = "Data/Reproducible-results/Figures/mirna_mirna_target_num_across_gene_size_read.svg",
       plot     = print(all_mirna_metrics_read_plot_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")


#Optimal gene size for MAD, SDE, and miRNA high metric----
#First bind the data frames together
methods_optimal_gene_read_df <- bind_rows(mad_cindices_read_df, 
                                          sde_cindices_read_df,
                                          mirna_high_cindices_read_df)

methods_optimal_gene_read_df <- methods_optimal_gene_read_df[,1:3]

methods_optimal_gene_read_df$method[119:177] <- rep("miRNA", 58)

#Saving the data frame
write.csv(methods_optimal_gene_read_df,
          "Data/Reproducible-results/Data/individual_optimal_gene_point_read_data.csv")

methods_optimal_gene_read_df <- methods_optimal_gene_read_df[,2:4]

#Basic miRNA-miRNA target plot
methods_optimal_gene_read_plot <- ggplot(data = methods_optimal_gene_read_df,
                                         aes(x=gene_sizes, y=c_index, color=method))

#Making the plot much nicer
methods_optimal_gene_read_plot_finished <-methods_optimal_gene_read_plot +
  geom_line(size=2.5) + geom_point(size=3.0)+
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", 
        legend.text = element_text(size = 35, face = "plain"),
        axis.text=element_text(size=25, face="bold"))+
  xlab("Gene Size")+ ylab("Mean Concordance Index")+
  ggtitle("Optimal Gene Number")+ 
  labs(color="Method")+
  scale_color_manual(values = c("#21908CFF","#FDE725FF", "#440154FF"))


#Saving the finished graph in .svg format
ggsave(filename = "Data/Reproducible-results/Figures/individual_metrics_optimal_gene_size_read.svg",
       plot     = print(methods_optimal_gene_read_plot_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")


#Random sample of genes (10 random samples taken) for READ.
#We will take the average of these 10 samplings. We will select the same number
#of genes as the number of genes that give our miRNA metric the best performance
#(350 genes)
random_read1 <- sample(colnames(cox_df), size = 350)
random_read2 <- sample(colnames(cox_df), size = 350)
random_read3 <- sample(colnames(cox_df), size = 350)
random_read4 <- sample(colnames(cox_df), size = 350)
random_read5 <- sample(colnames(cox_df), size = 350)
random_read6 <- sample(colnames(cox_df), size = 350)
random_read7 <- sample(colnames(cox_df), size = 350)
random_read8 <- sample(colnames(cox_df), size = 350)
random_read9 <- sample(colnames(cox_df), size = 350)
random_read10 <- sample(colnames(cox_df), size = 350)

#Saving the random results for reproducibility
write.csv(random_read1, "Data/Reproducible-results/Data/random_read_genes1.csv")
write.csv(random_read2, "Data/Reproducible-results/Data/random_read_genes2.csv")
write.csv(random_read3, "Data/Reproducible-results/Data/random_read_genes3.csv")
write.csv(random_read4, "Data/Reproducible-results/Data/random_read_genes4.csv")
write.csv(random_read5, "Data/Reproducible-results/Data/random_read_genes5.csv")
write.csv(random_read6, "Data/Reproducible-results/Data/random_read_genes6.csv")
write.csv(random_read7, "Data/Reproducible-results/Data/random_read_genes7.csv")
write.csv(random_read8, "Data/Reproducible-results/Data/random_read_genes8.csv")
write.csv(random_read9, "Data/Reproducible-results/Data/random_read_genes9.csv")
write.csv(random_read10, "Data/Reproducible-results/Data/random_read_genes10.csv")

#Putting all the gene vectors together in a list to loop over for cox model
random_read_gene_list <- list(random_read1, random_read2, random_read3,
                              random_read4, random_read5, random_read6,
                              random_read7, random_read8, random_read9,
                              random_read10)


all_random_gene_cindices <- seq(1,10, 1)

for(rg in all_random_gene_cindices){
  current_genes <- random_read_gene_list[[rg]]
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = current_genes,
                                cox.df = cox_df,
                                gene.num = 350,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/random_read_coefs_for_random_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the all_random_gene_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  all_random_gene_cindices[rg] <- top_cindex
  
}

write.csv(all_random_gene_cindices,
          "Data/Reproducible-results/Data/read_random_genes.csv")

mean_random <- mean(all_random_gene_cindices)


write.csv(read_methods_comp_df, "Data/Reproducible-results/Data/read_method_comparison_df.csv")



#MAD + SDE metric for READ
#We first do alpha weight optimization and then fit the elastic-net penalized
#cox model to every combination of alpha and assess its performance through 
#10-fold cross-validation

#Weight optimization
mad_sde_read_optimized <- two_weight_optimizer(first.metric = mad.genes,
                                               second.metric = sde.genes,
                                               my.filename = "Data/Reproducible-results/Data/mad_sde_read_optimized.rds")





#Penalized cox model----
c_index <- rep(0, 59)
gene_sizes <- seq(100, 3000, 50)

for(gs in gene_sizes){
  for(ms in mad_sde_read_optimized[1:11]){
    cox_model <- cox_model_fitter(my.seed = 1,
                                  cox.predictors = ms,
                                  cox.df = cox_df,
                                  gene.num = gs,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  my.filename = paste0("Data/Reproducible-results/Data/mad_sde_read_coefs_",gs,"_genes_",ms,"_index.csv"))
    
    #Getting the top concordance index from the cross validation and then rounding
    #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
    #the c_index list with the result
    top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
    c_index[which(gene_sizes==gs)] <- top_cindex
    
  }
  
  
}

mad_sde_finished_df <- as.data.frame(cbind(gene_sizes, c_index))

write.csv(mad_sde_finished_df, "Data/Reproducible-results/Data/mad_sde_read_df.csv")
  
#Now constructing a data frame of all the methods
read_methods_comp_df <- data.frame(Method=c("MAD", "SDE", 
                                            "miRNA", "MAD + SDE" ,
                                            "Random Genes"),
                                   c_index=c(0.6932, 0.7224, 
                                             0.7425,0.6932, mean_random))

#Factoring the levels to make the plot nicer
read_methods_comp_df$Method <- factor(read_methods_comp_df$Method,
                                      levels = c("miRNA", "MAD", "SDE",
                                                 "MAD + SDE",
                                                 "Random Genes"))

#Plotting the comparison of methods for READ
individual_graph_read <-ggplot(data = read_methods_comp_df,
                               aes(x=Method, y=c_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-READ",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40,
                                  family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey", size = 2.5, lineend = "round"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 34, family = "sans", angle = 45, vjust = 0.50),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none")

individual_graph_read + coord_cartesian(ylim = c(0.5,0.75))+ 
  scale_fill_manual(values = c("#FDE725FF","#404788FF", "#404788FF",
                               "#404788FF", "#404788FF")) 

individual_graph_read <- individual_graph_read + coord_cartesian(ylim = c(0.5,0.75))+ 
  scale_fill_manual(values = c("#FDE725FF","#404788FF", "#404788FF",
                               "#404788FF", "#404788FF")) 

#Saving the result
ggsave(filename = "Data/Reproducible-results/Figures/methods_comparison_read.svg",
       plot     = print(individual_graph_read, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")




#MAD + SDE metric for COAD
#We first do alpha weight optimization and then fit the elastic-net penalized
#cox model to every combination of alpha and assess its performance through 
#10-fold cross-validation

#Weight optimization
mad_sde_coad_optimized <- two_weight_optimizer(first.metric = mad.genes,
                                               second.metric = sde.genes,
                                               my.filename = "Data/Reproducible-results/Data/mad_sde_coad_optimized.rds")



#Penalized cox model----
c_index <- rep(0, 59)
gene_sizes <- seq(100, 3000, 50)

for(gs in gene_sizes){
  for(ms in mad_sde_coad_optimized[1:11]){
    cox_model <- cox_model_fitter(my.seed = 1,
                                  cox.predictors = ms,
                                  cox.df = cox_df,
                                  gene.num = gs,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  my.filename = paste0("Data/Reproducible-results/Data/mad_sde_coad_coefs_",gs,"_genes_",ms,"_index.csv"))
    
    #Getting the top concordance index from the cross validation and then rounding
    #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
    #the c_index list with the result
    top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
    c_index[which(gene_sizes==gs)] <- top_cindex
    
  }
  
  
}

mad_sde_finished_df <- as.data.frame(cbind(gene_sizes, c_index))

write.csv(mad_sde_finished_df, "Data/Reproducible-results/Data/mad_sde_coad_df.csv")

#Random sample of genes (10 random samples taken) for COAD.
#We will take the average of these 10 samplings. We will select the same number
#of genes as the number of genes that give our miRNA metric the best performance
#(1,900 genes)
random_coad1 <- sample(colnames(cox_df), size = 1900)
random_coad2 <- sample(colnames(cox_df), size = 1900)
random_coad3 <- sample(colnames(cox_df), size = 1900)
random_coad4 <- sample(colnames(cox_df), size = 1900)
random_coad5 <- sample(colnames(cox_df), size = 1900)
random_coad6 <- sample(colnames(cox_df), size = 1900)
random_coad7 <- sample(colnames(cox_df), size = 1900)
random_coad8 <- sample(colnames(cox_df), size = 1900)
random_coad9 <- sample(colnames(cox_df), size = 1900)
random_coad10 <- sample(colnames(cox_df), size = 1900)

#Saving the random results for reproducibility
write.csv(random_coad1, "Data/Reproducible-results/Data/random_coad_genes1.csv")
write.csv(random_coad2, "Data/Reproducible-results/Data/random_coad_genes2.csv")
write.csv(random_coad3, "Data/Reproducible-results/Data/random_coad_genes3.csv")
write.csv(random_coad4, "Data/Reproducible-results/Data/random_coad_genes4.csv")
write.csv(random_coad5, "Data/Reproducible-results/Data/random_coad_genes5.csv")
write.csv(random_coad6, "Data/Reproducible-results/Data/random_coad_genes6.csv")
write.csv(random_coad7, "Data/Reproducible-results/Data/random_coad_genes7.csv")
write.csv(random_coad8, "Data/Reproducible-results/Data/random_coad_genes8.csv")
write.csv(random_coad9, "Data/Reproducible-results/Data/random_coad_genes9.csv")
write.csv(random_coad10, "Data/Reproducible-results/Data/random_coad_genes10.csv")

#Putting all the gene vectors together in a list to loop over for cox model
random_coad_gene_list <- list(random_coad1, random_coad2, random_coad3,
                              random_coad4, random_coad5, random_coad6,
                              random_coad7, random_coad8, random_coad9,
                              random_coad10)

all_random_gene_cindices <- seq(1,10, 1)

for(rg in all_random_gene_cindices){
  current_genes <- random_coad_gene_list[[rg]]
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = current_genes,
                                cox.df = cox_df,
                                gene.num = 1900,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("Data/Reproducible-results/Data/random_coad_coefs_for_random_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the all_random_gene_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  all_random_gene_cindices[rg] <- top_cindex
  
}

write.csv(all_random_gene_cindices,
          "Data/Reproducible-results/Data/coad_random_genes.csv")

mean_random <- mean(all_random_gene_cindices)


#Now constructing a data frame of all the methods
coad_methods_comp_df <- data.frame(Method=c("MAD", "SDE", 
                                            "miRNA", "MAD + SDE" ,
                                            "Random Genes"),
                                   c_index=c(0.6395, 0.6721, 
                                             0.7076,0.6395, mean_random))

#Factoring the levels to make the plot nicer
coad_methods_comp_df$Method <- factor(coad_methods_comp_df$Method,
                                      levels = c("miRNA", "MAD", "SDE",
                                                 "MAD + SDE",
                                                 "Random Genes"))




write.csv(coad_methods_comp_df, "Data/Reproducible-results/Data/coad_method_comparison_df.csv")


#Plotting the comparison of methods for COAD
individual_graph_coad <-ggplot(data = coad_methods_comp_df,
                               aes(x=Method, y=c_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-COAD",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40,
                                  family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey", size = 2.5, lineend = "round"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 34, family = "sans", angle = 45, vjust = 0.50),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none")

individual_graph_coad + coord_cartesian(ylim = c(0.5,0.72))+ 
  scale_fill_manual(values = c("#FDE725FF","#404788FF", "#404788FF",
                               "#404788FF", "#404788FF")) 

individual_graph_coad <- individual_graph_coad + coord_cartesian(ylim = c(0.5,0.72))+
  scale_fill_manual(values = c("#FDE725FF","#404788FF", "#404788FF",
                               "#404788FF", "#404788FF"))

#Saving the result
ggsave(filename = "Data/Reproducible-results/Figures/methods_comparison_coad.svg",
       plot     = print(individual_graph_coad, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")




#CC Singlecell MS grid search COAD
#See server_speedup_coad.R for the code


#CC Singlecell MS grid search READ
#See server_speedup_read.R for the code


#Reading in the results files of the CC Singlecell MS COAD grid search
coad_ms_0 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0_cc_singlecell_ms_coad_used_combo_100_1010_index_5_df.csv")
colnames(coad_ms_0) <- c("number", "c_index")
coad_ms_0$mirna_num <- rep(seq(800,100,-100), each = 11)
coad_ms_0$mirna_target <- rep(seq(10, 1010, by=100), times = 8)

coad_ms_05 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.5_cc_singlecell_ms_coad_used_combo_100_1010_index_5_df.csv")
colnames(coad_ms_05) <- c("number", "c_index")
coad_ms_05 <- coad_ms_05[89:176,]
coad_ms_05$mirna_num <- rep(seq(800,100,-100), each = 11)
coad_ms_05$mirna_target <- rep(seq(10, 1010, by=100), times = 8)

coad_ms_1 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_1_cc_singlecell_ms_coad_used_combo_100_1010_index_5_df.csv")
colnames(coad_ms_1) <- c("number", "c_index")
coad_ms_1 <- coad_ms_1[177:264,]
coad_ms_1$mirna_num <- rep(seq(800,100,-100),each = 11)
coad_ms_1$mirna_target <- rep(seq(10, 1010, by=100), times = 8)


#Plotting the results of the CC Singlecell MS COAD grid search
#Alpha 0
heatmap_coad_ccs_ms_0 <- ggplot(data = coad_ms_0, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#Alpha 0.5
heatmap_coad_ccs_ms_05 <- ggplot(data = coad_ms_05, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#Alpha 1
heatmap_coad_ccs_ms_1 <- ggplot(data = coad_ms_1, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_1_test.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")



#Reading in the results files of the CC Singlecell MS READ grid search
read_ms_0 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_0_cc_singlecell_ms_read_used_combo_100_1010_index_8_df.csv")
colnames(read_ms_0) <- c("number", "c_index")
read_ms_0$mirna_num <- rep(seq(1000,100,-100), 11)
read_ms_0$mirna_target <- rep(seq(10, 1010, by=100), 10)


read_ms_05 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_0.5_cc_singlecell_ms_read_used_combo_100_1010_index_8_df.csv")
colnames(read_ms_05) <- c("number", "c_index")
read_ms_05 <- read_ms_05[111:220,]
read_ms_05$mirna_num <- rep(seq(1000,100,-100), 11)
read_ms_05$mirna_target <- rep(seq(10, 1010, by=100), 10)

read_ms_1 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_1_cc_singlecell_ms_read_used_combo_100_1010_index_11_df.csv")
colnames(read_ms_1) <- c("number", "c_index")
read_ms_1 <- read_ms_1[221:330,]
read_ms_1$mirna_num <- rep(seq(1000,100,-100), 11)
read_ms_1$mirna_target <- rep(seq(10, 1010, by=100), 10)



#Plotting the results of the CC Singlecell MS READ grid search
#READ MS Alpha 0
heatmap_read_ccs_ms_0 <- ggplot(data = read_ms_0, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS READ Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ms_0_ccs_ms_finished <- heatmap_read_ccs_ms_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_read_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_read_ms_0_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#READ MS Alpha 0.5
heatmap_read_ccs_ms_05 <- ggplot(data = read_ms_05, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS READ Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ms_05_ccs_ms_finished <- heatmap_read_ccs_ms_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_read_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_read_ms_05_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")



#READ MS Alpha 1
heatmap_read_ccs_ms_1 <- ggplot(data = read_ms_1, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS READ Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ms_1_ccs_ms_finished <- heatmap_read_ccs_ms_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_read_grid_search_heatmap_alpha_1.svg",
       plot     = print(heatmap_read_ms_1_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#CC Singlecell MM COAD
coad_mm_0 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0_cc_singlecell_mm_coad_used_combo_100_1010_index_5_df.csv")
colnames(coad_mm_0) <- c("number", "c_index")
coad_mm_0$mirna_num <- rep(seq(800,100,-100), 11)
coad_mm_0$mirna_target <- rep(seq(10, 1010, by=100), 8)

coad_mm_05 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.5_cc_singlecell_mm_coad_used_combo_100_1010_index_10_df.csv")
colnames(coad_mm_05) <- c("number", "c_index")
coad_mm_05$mirna_num <- rep(seq(800,100,-100), 11)
coad_mm_05$mirna_target <- rep(seq(10, 1010, by=100), 8)


coad_mm_1 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_1_cc_singlecell_mm_coad_used_combo_100_1010_index_10_df.csv")
colnames(coad_mm_1) <- c("number", "c_index")
coad_mm_1$mirna_num <- rep(seq(800,100,-100), 11)
coad_mm_1$mirna_target <- rep(seq(10, 1010, by=100), 8)



#Alpha 0
heatmap_coad_ccs_mm_0 <- ggplot(data = coad_mm_0, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MM COAD Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_mm_finished <- heatmap_coad_ccs_mm_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mm_coad_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_coad_ccs_mm_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.5
heatmap_coad_ccs_mm_05 <- ggplot(data = coad_mm_05, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MM COAD Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans", margin=margin(0,20,0,0)),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_mm_finished <- heatmap_coad_ccs_mm_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mm_coad_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_coad_ccs_mm_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")




#Alpha 1
heatmap_coad_ccs_mm_1 <- ggplot(data = coad_mm_1, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MM COAD Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans", margin=margin(0,20,0,0)),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_mm_finished <- heatmap_coad_ccs_mm_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mm_coad_grid_search_heatmap_alpha_1.svg",
       plot     = print(heatmap_coad_ccs_mm_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")






#CC Singlecell MM READ
read_mm_0 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_0_cc_singlecell_mm_read_used_combo_100_1010_index_1_df.csv")
colnames(read_mm_0) <- c("number", "c_index")
read_mm_0$mirna_num <- rep(seq(1000,100,-100), 11)
read_mm_0$mirna_target <- rep(seq(10, 1010, by=100), 10)

read_mm_05 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_0.5_cc_singlecell_mm_read_used_combo_100_1010_index_7_df.csv")
colnames(read_mm_05) <- c("number", "c_index")
read_mm_05$mirna_num <- rep(seq(1000,100,-100), 11)
read_mm_05$mirna_target <- rep(seq(10, 1010, by=100), 10)


read_mm_1 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_1_cc_singlecell_mm_read_used_combo_100_1010_index_10_df.csv")
colnames(read_mm_1) <- c("number", "c_index")
read_mm_1$mirna_num <- rep(seq(1000,100,-100), 11)
read_mm_1$mirna_target <- rep(seq(10, 1010, by=100), 10)



#Alpha 0
heatmap_read_ccs_mm_0 <- ggplot(data = read_mm_0, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MM READ Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans", margin = margin(0,20,0,0)),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ccs_mm_finished <- heatmap_read_ccs_mm_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mm_read_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_read_ccs_mm_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.5
heatmap_read_ccs_mm_05 <- ggplot(data = read_mm_05, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MM READ Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ccs_mm_finished <- heatmap_read_ccs_mm_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mm_read_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_read_ccs_mm_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#Alpha 1
heatmap_read_ccs_mm_1 <- ggplot(data = read_mm_1, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MM READ Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ccs_mm_finished <- heatmap_read_ccs_mm_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mm_read_grid_search_heatmap_alpha_1.svg",
       plot     = print(heatmap_read_ccs_mm_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#CC Singlecell MMS COAD
coad_mms_0 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0_cc_singlecell_mms_coad_used_combo_100_1010_index_86_df.csv")
colnames(coad_mms_0) <- c("number", "c_index")
coad_mms_0$mirna_num <- rep(seq(800,100,-100), each=11)
coad_mms_0$mirna_target <- rep(seq(10, 1010, by=100),times =8)

coad_mms_05 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.5_cc_singlecell_mms_coad_used_combo_100_1010_index_17_df.csv")
colnames(coad_mms_05) <- c("number", "c_index")
coad_mms_05$mirna_num <- rep(seq(800,100,-100), each = 11)
coad_mms_05$mirna_target <- rep(seq(10, 1010, by=100), times = 8)


coad_mms_1 <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_1_cc_singlecell_mms_coad_used_combo_100_1010_index_17_df.csv")
colnames(coad_mms_1) <- c("number", "c_index")
coad_mms_1$mirna_num <- rep(seq(800,100,-100), each = 11)
coad_mms_1$mirna_target <- rep(seq(10, 1010, by=100), times = 8)



#Alpha 0
heatmap_coad_ccs_mms_0 <- ggplot(data = coad_mms_0, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MMS COAD Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_mms_finished <- heatmap_coad_ccs_mms_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mms_coad_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_coad_ccs_mms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.5
heatmap_coad_ccs_mms_05 <- ggplot(data = coad_mms_05, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MMS COAD Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans", margin=margin(0,20,0,0)),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_mms_finished <- heatmap_coad_ccs_mms_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mms_coad_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_coad_ccs_mms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")




#Alpha 1
heatmap_coad_ccs_mms_1 <- ggplot(data = coad_mms_1, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MMS COAD Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans", margin=margin(0,20,0,0)),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_mms_finished <- heatmap_coad_ccs_mms_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mms_coad_grid_search_heatmap_alpha_1.svg",
       plot     = print(heatmap_coad_ccs_mms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")


#CC Singlecell MMS READ
read_mms_0 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_0_cc_singlecell_mms_read_used_combo_900_1010_index_56_df.csv")
colnames(read_mms_0) <- c("number", "c_index")
read_mms_0$mirna_num <- rep(seq(1000,100,-100), each=11)
read_mms_0$mirna_target <- rep(seq(10, 1010, by=100), times=10)

read_mms_05 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_0.5_cc_singlecell_mms_read_used_combo_100_1010_index_7_df.csv")
colnames(read_mms_05) <- c("number", "c_index")
read_mms_05$mirna_num <- rep(seq(1000,100,100), 11)
read_mms_05$mirna_target <- rep(seq(10, 1010, by=100), 10)


read_mms_1 <- read.csv("Data/Reproducible-results/Data/Outputs/READ/top_cindices_alpha_1_cc_singlecell_mms_read_used_combo_100_1010_index_10_df.csv")
colnames(read_mms_1) <- c("number", "c_index")
read_mms_1$mirna_num <- rep(seq(1000,100,100), 11)
read_mms_1$mirna_target <- rep(seq(10, 1010, by=100), 10)



#Alpha 0
heatmap_read_ccs_mms_0 <- ggplot(data = read_mms_0, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MMS READ Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans", margin = margin(0,20,0,0)),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ccs_mms_finished <- heatmap_read_ccs_mms_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mms_read_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_read_ccs_mms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.5
heatmap_read_ccs_mms_05 <- ggplot(data = read_mms_05, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MMS READ Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ccs_mms_finished <- heatmap_read_ccs_mms_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mms_read_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_read_ccs_mms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#Alpha 1
heatmap_read_ccs_mms_1 <- ggplot(data = read_mms_1, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MMS READ Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_read_ccs_mms_finished <- heatmap_read_ccs_mms_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_mms_read_grid_search_heatmap_alpha_1.svg",
       plot     = print(heatmap_read_ccs_mms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")




#Additional alpha values for glmnet on COAD for CC Singlecell MS
alpha_01_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.1_cc_singlecell_ms_coad_df.csv")
alpha_02_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.2_cc_singlecell_ms_coad_df.csv")
alpha_03_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.3_cc_singlecell_ms_coad_df.csv")
alpha_04_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.4_cc_singlecell_ms_coad_df.csv")
alpha_06_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.6_cc_singlecell_ms_coad_df.csv")
alpha_07_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.7_cc_singlecell_ms_coad_df.csv")
alpha_08_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.8_cc_singlecell_ms_coad_df.csv")
alpha_09_df <- read.csv("Data/Reproducible-results/Data/Outputs/COAD/top_cindices_alpha_0.9_cc_singlecell_ms_coad_df.csv")


colnames(alpha_01_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_02_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_03_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_04_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_06_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_07_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_08_df) <- c("X", "mirna_num", "mirna_target", "c_index")
colnames(alpha_09_df) <- c("X", "mirna_num", "mirna_target", "c_index")

alpha_01_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_02_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_03_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_04_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_06_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_07_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_08_df$mirna_num <- rep(seq(800,100,100), each=11)
alpha_09_df$mirna_num <- rep(seq(800,100,100), each=11)



#Alpha 0.1
heatmap_coad_ccs_ms_01 <- ggplot(data = alpha_01_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_01 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_01.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.2
heatmap_coad_ccs_ms_02 <- ggplot(data = alpha_02_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.2",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_02 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_02.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.3
heatmap_coad_ccs_ms_03 <- ggplot(data = alpha_03_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.3",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_03 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_03.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.4
heatmap_coad_ccs_ms_04 <- ggplot(data = alpha_04_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.4",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_04 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_04.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")



#Alpha 0.6
heatmap_coad_ccs_ms_06 <- ggplot(data = alpha_06_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.6",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_06 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_06.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.7
heatmap_coad_ccs_ms_07 <- ggplot(data = alpha_07_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.7",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_07 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_07.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha 0.8
heatmap_coad_ccs_ms_08 <- ggplot(data = alpha_08_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.8",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_08 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_08.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")

#Alpha 0.9
heatmap_coad_ccs_ms_09 <- ggplot(data = alpha_09_df, aes(x=mirna_num, y=mirna_target, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "CC Singlecell MS COAD Alpha = 0.9",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_ccs_ms_09 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_grid_search_heatmap_alpha_09.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")





#KM risk calculation for COAD----
#Now loading the top performing result (alpha = 1, 700 miRNA and 1010 miRNA targets)
mirna_sde_optimized <- readRDS("~/Desktop/Optimization_700_1010_targets_cc_singlecell_ms_1_alpha_coad.rds")
data_set <- "coad"
a <- 0.5
counter <- 1
my_cindices <-c()

for(ms in mirna_sde_optimized[1:11]){
  cox_model <- cox_model_fitter(my.seed = 1,
                                my.alpha = 0.5,
                                my.dataset = "COAD",
                                cox.predictors = ms,
                                cox.df = cox_df,
                                gene.num = 1900,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("~/Desktop/cc_singlecell_ms_",data_set,"_alpha_",a,"_coefs__index.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the c_index list with the result
  current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  my_cindices[counter] <- current_cindex
  counter <- counter + 1
  
}

top_cindex <-max(my_cindices)
top_index <- which(my_cindices==top_cindex)
print(top_index)
print(top_cindex)



patient_risk <- risk_score_calculator(my.file = "~/Desktop/cc_singlecell_ms_coad_alpha_1_coefs__index.csv",
                                      tumor.data = FALSE, n.data = FALSE,
                                      set.ci = TRUE, cox.df = cox_df, 
                                      set.test = TRUE, my.km.plot = "test_coad_km.svg",
                                      plot.title = "TCGA-COAD")
patient_risk$`KM Plot`


#Saving the KM plot to .svg format
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_km_plot.svg",
       plot     = print(patient_risk$`KM Plot`, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")


#KM risk calculation for READ----
patient_risk <- risk_score_calculator(my.file = "~/Desktop/cc_singlecell_ms_read_alpha_1_coefs__index.csv",
                                      tumor.data = FALSE, n.data = FALSE,
                                      set.ci = TRUE,
                                      cox.df = cox_df, 
                                      plot.title = "TCGA-READ")
patient_risk$`KM Plot`

#Saving the KM plot to .svg format
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_read_km_plot.svg",
       plot     = print(patient_risk$`KM Plot`$plot, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 34, height = 34,
       units    = "cm")


#For plotting the COAD coefficients----
coad_coef_df <- read.csv("~/Desktop/cc_singlecell_ms_coad_alpha_1_coefs__index.csv")
coad_coef_df <- coad_coef_df[,2:3]
colnames(coad_coef_df) <- c("Gene", "coefs")
coad_coef_df$hazard_ratio <- exp(coad_coef_df$coefs)
coad_coef_df$effect_size <- ifelse(coad_coef_df$hazard_ratio>1,
                                   abs(1- coad_coef_df$hazard_ratio)*100,
                                   (1 - coad_coef_df$hazard_ratio)*100)

coad_coef_df_sub <- filter(coad_coef_df, effect_size>80.0)

#These top 10 active genes are associated with an FDR corrected p-value pathway
#of NOTCH3 which has been implicated in colon cancer here from reactome: 
#1. A third Notch in colorectal cancer progression and metastasis (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7537388/)
#2. AKT-dependent NOTCH3 activation drives tumor progression in a model of mesenchymal colorectal cancer (https://pubmed.ncbi.nlm.nih.gov/32749453/)
#3. The miR-1-NOTCH3-Asef pathway is important for colorectal tumor cell migration (https://pubmed.ncbi.nlm.nih.gov/24244701/)
#4. Role of Notch signaling in colorectal cancer (https://pubmed.ncbi.nlm.nih.gov/19793799/)
#The FDR corrected p-values are seen in the report. 

coad_coef_df_sub <- filter(coad_coef_df_sub, Gene != "ZNF705D")

coad_coef_plot <- ggplot(data = coad_coef_df_sub, aes(x=Gene, y=effect_size,
                                                      color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("COAD Coefficients")+
  ylab("Effect Size (%)")+
  xlab("Genes")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40,
                                  family = "sans"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 30, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"))+
  scale_color_viridis_d(direction = -1)+
  scale_fill_viridis_d(direction = -1)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coad_coef_plot

#Saving the KM plot to .svg format
ggsave(filename = "Data/Reproducible-results/Figures/cc_singlecell_ms_coad_coef_plot_no_znf705.svg",
       plot     = print(coad_coef_plot, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 40, height = 40,
       units    = "cm")


#Citation for ZNF705D at: Authentication of differential gene expression in oral squamous cell carcinoma using machine learning applications
#(https://bmcoralhealth.biomedcentral.com/articles/10.1186/s12903-021-01642-9)

#Citation for TAL2 at: TAL2, a helix-loop-helix gene activated by the (7;9)(q34;q32) translocation in human T-cell leukemia
# (https://www.pnas.org/doi/abs/10.1073/pnas.88.24.11416)

#Citation for ST6GALNAC3 at:Comprehensive analysis of coexpressed long noncoding RNAs and genes in breast cancer
# (https://obgyn.onlinelibrary.wiley.com/doi/abs/10.1111/jog.13840?casa_token=3sf94uCg6UUAAAAA:50KHN5Q3La4JCZkKTugJBA8YLia2mJEFnGjhJKQw5ZEshIwXSh_pIMPKOTHfeOevB6V8R_Eysd_E2Ks)
# and at: Biomarker potential of ST6GALNAC3 and ZNF660 promoter hypermethylation in prostate cancer tissue and liquid biopsies
# (https://febs.onlinelibrary.wiley.com/doi/full/10.1002/1878-0261.12183)

#Citation for HTR2C at: A zebrafish embryo screen utilizing gastrulation identifies the HTR2C inhibitor pizotifen as a suppressor of EMT-mediated metastasis
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8824480/)

#Citation for HS6ST3 at: Silencing HS6ST3 inhibits growth and progression of breast cancer cells through suppressing IGF1R and inducing XAF1
# (https://www.sciencedirect.com/science/article/pii/S0014482716304281?casa_token=zKk2Dvoylr8AAAAA:VZVqZE-BKry_Mg6fZXidN-ZMjdLKL2QfBWKKuH126dhlnYDOz0-Mkru0WRu8Ud1UHuvHUiAVZuY)

#Citation for GRIK3 at: CircASXL1 knockdown represses the progression of colorectal cancer by downregulating GRIK3 expression by sponging miR-1205
# (https://link.springer.com/article/10.1186/s12957-021-02275-6)

#Citation for FABP7 at: FABP7 promotes cell proliferation and survival in colon cancer through MEK/ERK signaling pathway
# (https://www.sciencedirect.com/science/article/pii/S0753332218327021)

#Citation for AJAP1 at: Tumor-associated methylation of the putative tumor suppressor AJAP1 gene and association between decreased AJAP1 expression and shorter survival in patients with glioma
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4013351/)

#Citation for AC079612.1 at: Mechanism of lnc-AC079612.1.1-1∶1 in hepatic metastasis of colon cancer (this looks like a really low quality journal)
# (https://pesquisa.bvsalud.org/portal/resource/pt/wpr-693273)

#Now that we see the active coefficients and their hazard ratios we will attempt
#to use a relatively new technique to get p-values for the alpha value of our
#best fit model and attempt inference for both p-values of individual active
#genes and 95% confidence intervals
betas <- as.numeric(coad_alpha1$Active.Coefficients)
my_x <- read.csv("~/Desktop/x_matrix.csv")
rownames(my_x) <- my_x$X
my_x$X <- NULL
my_x$X.Intercept. <- NULL
my_x <- as.matrix(my_x)
my_y <- read.csv("~/Desktop/y_matrix.csv")
time <- my_y$time
my_status <- my_y$status

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}


#Test Cox data to see if I can get the p-values from it
data("CoxExample")
x <- CoxExample$x
y <- CoxExample$y

x = scale(x, TRUE, TRUE)

fit <- glmnet(x, y, family = "cox", standardize = FALSE)


set.seed(1)
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C")
lambda <- cvfit$lambda.min
beta <- as.numeric(coef(fit, x=x, y=y, s=lambda/1000, exact = TRUE))
time <- y[,"time"]
my_status <- y[,"status"]

lass_coefs_coad <- fixedLassoInf(x = x, time, beta = betas, status = my_status,
              lambda = 0.01073, family = "cox", verbose = TRUE)

lass_coefs_ex <- fixedLassoInf(x = x, time, beta = beta, status = my_status,
                            lambda = lambda, family = "cox", verbose = TRUE)

#Got really side-tracked with this because I was trying to see if I could use significance to improve the interpretability of the model.
#Still don't have this working yet. It might be a really good way to make the model more parsimonious than the 156 active genes we have
#currently. This paper by Taylor et al. "Statistical learning and selective inference" (https://www.pnas.org/doi/full/10.1073/pnas.1507583112)
#details the problem of selective inference after using other methods to select predictors and the challenges that are currently facing statisticians
#as they do this. "Tractable Post-Selection Maximum Likelihood Inference for the Lasso" by Meir et al. also details that in essence
#"Applying standard statistical methods after model selection may yield inefficient estimators and hypothesis tests that fail to achieve nominal type-I error rates.
#The main issue is the fact that the post-selection distribution of the data differs from the original distribution. In particular, the observed data is constrained
#to lie in a subset of the original sample space that is determined by the selected model.
#This often makes the post-selection likelihood of the observed data intractable and maximum likelihood inference difficult." (https://arxiv.org/pdf/1705.09417.pdf)
#(https://stats.stackexchange.com/questions/241082/testing-for-coefficients-significance-in-lasso-logistic-regression)
#This other paper by Lockhart et al. "A SIGNIFICANCE TEST FOR THE LASSO" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4285373/)
#attempts to implement significance testing for lasso. If we can present confidence intervals and significance results for our method I think it would
#make our model more rigorous as currently we can't say anything about the individual predictor's significance 

#miRNA grid search at different alpha values
mirna_alpha0 <- read.csv("~/Desktop/alpha_0_cindices.csv")
mirna_alpha05 <- read.csv("~/Desktop/alpha_0.5_cindices.csv")
mirna_alpha1 <- read.csv("~/Desktop/alpha_1_cindices.csv")


colnames(mirna_alpha0) <- c("mirna_num", "c_index")
mirna_alpha0$mirna_num <- rep(seq(800,100,-100), each=11)
mirna_alpha0$mirna_targets <- rep(seq(10,1010,100), times=8)


colnames(mirna_alpha05) <- c("mirna_num", "c_index")
mirna_alpha05$mirna_num <- rep(seq(800,100,-100), each=11)
mirna_alpha05$mirna_targets <- rep(seq(10,1010,100), times=8)

colnames(mirna_alpha1) <- c("mirna_num", "c_index")
mirna_alpha1$mirna_num <- rep(seq(800,100,-100), each=11)
mirna_alpha1$mirna_targets <- rep(seq(10,1010,100), times=8)



#miRNA only heat map
#Alpha = 0
heatmap_coad_mirna_0 <- ggplot(data = mirna_alpha0, aes(x=mirna_num, y=mirna_targets, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "miRNA COAD Alpha = 0",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_mirna_0 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/mirna_coad_grid_search_heatmap_alpha_0.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")



#Alpha = 0.5
heatmap_coad_mirna_05 <- ggplot(data = mirna_alpha05, aes(x=mirna_num, y=mirna_targets, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "miRNA COAD Alpha = 0.5",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_mirna_05 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/mirna_coad_grid_search_heatmap_alpha_05.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")


#Alpha = 1
heatmap_coad_mirna_1 <- ggplot(data = mirna_alpha1, aes(x=mirna_num, y=mirna_targets, fill=c_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(c_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of miRNAs",
       y = "# of miRNA Targets",
       title = "miRNA COAD Alpha = 1",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.0,"cm"))

#Changing to color-blind friendly palette
heatmap_coad_ccs_ms_finished <- heatmap_coad_mirna_1 + scale_fill_viridis_c()

#Now saving the heat map
ggsave(filename = "Data/Reproducible-results/Figures/mirna_coad_grid_search_heatmap_alpha_1.svg",
       plot     = print(heatmap_coad_ccs_ms_finished, newpage = FALSE),
       device   = "svg", dpi=300,
       width    = 32, height = 32,
       units    = "cm")








#For plotting the READ coefficients-----
read_coef_df <- read.csv("~/Desktop/read_active_genes.csv")
read_coef_df <- read_coef_df[,2:3]
colnames(read_coef_df) <- c("Gene", "coefs")
read_coef_df_sub <- filter(read_coef_df, abs(coefs)>0.5)


read_coef_plot <- ggplot(data = read_coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("READ Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"))+
  scale_color_viridis_d(direction = -1)+
  scale_fill_viridis_d(direction = -1)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


read_coef_plot


#Individual metrics vs. CC Singlecell MMS COAD----
#MAD + SDES
load("Data/TCGA-COAD/MAD/mad.RData", verbose = TRUE)
load("Data/TCGA-COAD/SDE/sde.RData", verbose = TRUE)

mad_sdes_optimized <- two_weight_optimizer(first.metric = mad.genes,
                                           second.metric = sde.genes,
                                           my.filename = "Data/Reproducible-results/Data/mad_sdes_coad_optimized.RData")


cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1

for (x in mad_sdes_optimized[1:11]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = 2100,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/coad_mad_sde_active_genes.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  current_coefs <- length(current_cox$`Active Coefficients`)
  my_active_coefs <- c(my_active_coefs, current_coefs)
}






cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1

for (x in mad_sdes_optimized[1:11]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = 1800,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/mad_sdes_active_genes.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  current_coefs <- length(current_cox$`Active Coefficients`)
  my_active_coefs <- c(my_active_coefs, current_coefs)
}



individual_metrics_df <- read.csv(file = "Data/TCGA-COAD/indivdual_metrics_cindex_updated.csv")
individual_metrics_df$Method <- factor(individual_metrics_df$Method,
                                       levels = c("CCS MM",
                                                  "CCS MS",
                                                  "CCS MMS",
                                                  "MAD + SDE",
                                                  "MAD", 
                                                  "MiRNA",
                                                  "SDE", 
                                                  "Random"))
individual_graph_coad <-ggplot(data = individual_metrics_df,
                               aes(x=Method, y=C_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-COAD",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40,
                                  family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 34, family = "sans", angle = 45, vjust = 0.58),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none")

individual_graph_coad + coord_cartesian(ylim = c(0.5,0.71))+ 
  scale_fill_manual(values = c("#FDE725FF","#FDE725FF", "#FDE725FF", "#404788FF", 
                            "#404788FF", "#404788FF", "#404788FF", "#404788FF")) 



#Individual metrics vs. CC Singlecell MMS READ----
#MAD + SDES
load("Data/TCGA-READ/MAD/mad.RData", verbose = TRUE)
load("Data/TCGA-READ/SDE/sde.RData", verbose = TRUE)

mad_sdes_optimized <- two_weight_optimizer(first.metric = mad.genes,
                                           second.metric = sde.genes,
                                           my.filename = "Data/Reproducible-results/Data/mad_sdes_read_optimized.RData")


cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1

for (x in mad_sdes_optimized[1:11]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = 1450,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/read_mad_sde_active_genes.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  current_coefs <- length(current_cox$`Active Coefficients`)
  my_active_coefs <- c(my_active_coefs, current_coefs)
}






cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1

for (x in mad_sdes_optimized[1:11]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = 1800,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/mad_sdes_active_genes.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  current_coefs <- length(current_cox$`Active Coefficients`)
  my_active_coefs <- c(my_active_coefs, current_coefs)
}



individual_metrics_df <- read.csv(file = "Data/TCGA-READ/individual_vs_cc_singlecell_mms_comparison_updated.csv")
individual_metrics_df$Method <- factor(individual_metrics_df$Method,
                                       levels = c("CCS MM",
                                                  "CCS MS",
                                                  "CCS MMS",
                                                  "MAD + SDE",
                                                  "MAD", 
                                                  "MiRNA",
                                                  "SDE", 
                                                  "Random"))
individual_graph_read <-ggplot(data = individual_metrics_df,
                               aes(x=Method, y=C_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-READ",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40,
                                  family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 34, family = "sans", angle = 45, vjust = 0.58),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none")

individual_graph_read + coord_cartesian(ylim = c(0.5,0.85))+ 
  scale_fill_manual(values = c("#FDE725FF","#FDE725FF", "#FDE725FF", "#404788FF", 
                               "#404788FF", "#404788FF", "#404788FF", "#404788FF")) 


