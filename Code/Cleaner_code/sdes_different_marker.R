#sdes_different_marker.R

#Load the different marker files----
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_snai1.RData", verbose = TRUE)
snai1_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_twist1.RData", verbose = TRUE)
twist1_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_twist1-redo.RData", verbose = TRUE)
twist1_redo_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_zeb1.RData", verbose = TRUE)
zeb1_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde.RData", verbose = TRUE)
vim_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_cdh1.RData", verbose = TRUE)
cdh1_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_cdh2.RData", verbose = TRUE)
cdh2_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_grhl2.RData", verbose = TRUE)
grhl2_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_ccnd1.RData", verbose = TRUE)
ccnd1_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_cdk4.RData", verbose = TRUE)
cdk4_sde <- sde.genes
load("Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_cdk6.RData", verbose = TRUE)
cdk6_sde <- sde.genes

#Loading needed packages----
library(survival);packageVersion("survival")
library(survminer);packageVersion("survminer")

#Other required packages----
require(BiocGenerics)
require(doParallel)
require(glmnet)
require(parallel)
require(survival)

#Setting the number of processors on the----
#machine to speed up the fitting
num_of_cores <- parallel::detectCores()
registerDoParallel(cores = num_of_cores)


my_gene_list <- list(snai1=snai1_sde, twist1=twist1_sde, twist1_redo=twist1_redo_sde, zeb1=zeb1_sde, vim=vim_sde, cdh1=cdh1_sde, cdh2=cdh2_sde, grhl2=grhl2_sde, ccnd1=ccnd1_sde, cdk4=cdk4_sde, cdk6=cdk6_sde)


#Getting access to cox function from other script----
if(!exists("cox_model_fitter", mode="function")) source("cox_model.R")

#Setting the gene size to a constant value----
cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1

for (y in my_gene_list) {
  print(counter)
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1200, cox.predictors = y , tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  
}

my_df <- data.frame(marker=names(my_gene_list), concordance_index=my_cindicies)
my_df$se <- sqrt(my_df$concordance_index/sqrt(length(my_df$concordance_index)))


my_finished_df <- read.csv("~/Desktop/different_markers_coad_read_and_combined.csv")
my_finished_df$marker <- factor(my_finished_df$marker, levels = c("ccnd1","cdk4","cdk6", "cdh1", "cdh2", "grhl2", "snai1", "twist1", "twist1_redo", "vim", "zeb1"))
my_finished_df$dataset <- factor(my_finished_df$dataset, levels = c("TCGA-COAD", "TCGA-READ", "TCGA-COAD + TCGA-READ"))

res_aov <- aov(concordance_index~dataset, data = my_finished_df)
aov_sum <- summary(res_aov)
pvalue_to_plot_dataset <- round(aov_sum[[1]][["Pr(>F)"]][1], digits = 5)



dataset_panel <- ggplot(my_finished_df, aes(x=dataset, y=concordance_index, fill=dataset))+
  geom_boxplot(data=my_finished_df, aes(x=dataset, y=concordance_index))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("Datasets")+
  xlab("Dataset")+
  ylab("Concordance Index")


dataset_panel <- dataset_panel + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
  annotate(geom = 'text', label = paste0("P=", pvalue_to_plot_dataset), x = 0.68, y = 0.62, hjust = 0, vjust = 1, size=3)



marker_panel <- ggplot(my_finished_df, aes(x=marker, y=concordance_index, fill=marker))+
  geom_boxplot(data=my_finished_df, aes(x=marker, y=concordance_index))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("Pseudo-time Markers")+
  xlab("Markers")+
  ylab("Concordance Index")



marker_panel <- marker_panel + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
  annotate(geom = 'text', label = paste0("P=", pvalue_to_plot_marker), x = 1, y = 0.62, hjust = 0, vjust = 1, size=3)




# p+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
#   stat_compare_means(method = "anova", label.x = 1, label.y = 0.60)+
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "CC Singlecell")


ggarrange(marker_panel, dataset_panel, nrow = 2, ncol = 1, labels = c("A.", "B."))

