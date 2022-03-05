#Name: reproducible.R
#Purpose: For reproducing all figures & results in the manuscript
#Author: Andrew Willems <awillems@vols.utk.edu>

#Loading needed packages----
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library(svglite)
library(TCGAbiolinks)
library(tidyverse)
library(viridis)

#Loading the needed functions from their respecitive files----
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
#cc_tumor_fpkm <- as.matrix(cc_tumor_fpkm)

#First pre-processing the scRNA-seq data before sending it to MAGIC
#Keeping genes expressed in at least 10 cells
keep_rows <- rowSums(cc_tumor_fpkm > 0) > 10
cc_tumor_fpkm <- cc_tumor_fpkm[keep_rows,]

# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=colSums(cc_tumor_fpkm)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

#Normalizing the library size
cc_tumor_fpkm_normalized <- library.size.normalize(cc_tumor_fpkm)
cc_tumor_fpkm_normalized <- sqrt(cc_tumor_fpkm_normalized)


#Denoising the single-cell data and saving the output----
cc_tumor_fpkm <- magic_denoiser(sc.data = cc_tumor_fpkm_normalized,
                                magic.seed = 123,magic.solver = 'approximate')
saveRDS(cc_tumor_fpkm,
        file = "Data/Reproducible-results/denoised-colon-and-rectal-single-cell-data_v2.rds")


#Generating the VIM based pseudotime progression with Monocle3----
#Based on the expression of the VIM graph we select the bottom-rightmost point
#as the root of our pseudotime as this matches the gradient of VIM expression
#level
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP"),
                                   cell.data = cc_tumor_fpkm$denoised_sc_dataframe,
                                   cell.meta = cc_tumor_fpkm$cds_gene_names)


#MAD metric for colon and rectal cancer----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
saveRDS(mad.genes,
        file = "Data/Reproducible-results/Data/mad_colon_and_rectal_cancer_v2.rds")

#Switchde metric for colon and rectal cancer----
sde.genes <- switchde_calculator(cc_tumor_fpkm$denoised_sc_dataframe,
                                 pseudo.time = cds_output$Pseudotime)
saveRDS(sde.genes,
        file = "Data/Reproducible-results/Data/sde_colon_and_rectal_cancer_v2.rds")

#MiRNA metric----
#Due to it taking a while I am just including the file that we get from this
#metric. I can run it overnight so that we can show that the files are the same

# mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
#                                 cancer.type2 = "colon cancer",
#                                 max.miR.targets = 810,
#                                 cancer.up = TRUE,
#                                 mirna.remove = c("hsa-miR-129-2-3p",
#                                                  "hsa-miR-129-1-3p",
#                                                  "hsa-miR-454-3p",
#                                                  "hsa-miR-365a-5p"),
#                                 max.mirnas = 500,
#                                 ts.org = "Human",
#                                 ts.version = "7.2",
#                                 print.ts.targets = TRUE,
#                                 save.mirna.genes = TRUE,
#                                 mirna.gene.filename = "Data/Reproducible-results/mirna_genes_global_search_mirna_500_810_targets.csv",
#                                 mirna.gene.rfile = "Data/Reproducible-results/mirna_genes_global_search_mirna_500_810_targets.RData")

mirna.ranking <- readRDS("Data/Reproducible-results/Data/500_mirna_810_targets_colon_and_rectal_cancer_mirna_genes.rds")

#Optimizing the weights of the three metric linear model----
#Due to this step taking a while I am including its output and will run this
#overnight and show that the outputs are the same
# mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric  = mad.genes,
#                                                    second.metric = mirna.ranking,
#                                                    third.metric  = sde.genes,
#                                                    my.filename   = "~/Desktop/Optimized-model.RData")


load("Data/Reproducible-results/Data/Optimization_500_810_targets_for_colon_cancer.RData",
     verbose = TRUE)


mad_sdes_mirna_optimized <- integrated_gene_lists



#Elastic-net penalized cox model----
cox_df <- readRDS("Data/TCGA-COAD/coad_df_finished_v2.rds")
cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1

#Just showing the top performing result from the grid search for colon cancer
#here
for (x in mad_sdes_mirna_optimized[1:121]) {
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
                                  my.filename = "~/Desktop/coad_active_genes.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build
  #the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  current_coefs <- length(current_cox$`Active Coefficients`)
  my_active_coefs <- c(my_active_coefs, current_coefs)
}



#Just showing the top performing result from the grid search for rectal cancer
#here
load("Data/Reproducible-results/Data/Optimization_800_510_targets_for_rectal_cancer.RData",
     verbose = TRUE)

mad_sdes_mirna_optimized <- integrated_gene_lists

cox_df <- readRDS("Data/TCGA-READ/read_df_finished_v2.rds")

cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1

for (x in mad_sdes_mirna_optimized[15]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = 1100,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/read_active_genes.csv") 
  
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

#KM risk calculation for COAD----
patient_risk <- risk_score_calculator(my.file = "~/Desktop/coad_active_genes.csv",
                                      my.title = "TCGA-COAD",
                                      tumor.data = FALSE,
                                      n.data = FALSE, 
                                      cox.df = cox_df,
                                      show.pval = TRUE,
                                      show.pval.method = FALSE)
patient_risk


#KM risk calculation for READ----
patient_risk <- risk_score_calculator(my.file = "~/Desktop/read_active_genes.csv",
                                      tumor.data = FALSE,
                                      n.data = FALSE, 
                                      cox.df = cox_df, 
                                      plot.title = "TCGA-READ Test")
patient_risk



#For plotting the COAD coefficients----
coad_coef_df <- read.csv("~/Desktop/coad_active_genes.csv")
coad_coef_df <- coad_coef_df[,2:3]
colnames(coad_coef_df) <- c("Gene", "coefs")
coad_coef_df_sub <- filter(coad_coef_df, abs(coefs)>0.5)


coad_coef_plot <- ggplot(data = coad_coef_df_sub, aes(x=Gene, y=coefs,
                                                      color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("COAD Coefficients")+
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


coad_coef_plot

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


