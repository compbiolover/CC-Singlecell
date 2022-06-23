#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")


deseq2_genes <- readRDS(file = "deseq2_top_genes.rds")
#cox_df <-readRDS("coad_df_finished.rds")
cox_df <-readRDS("read_df_finished.rds")

#For finding the optimal number of genes
gene_sizes <- seq(100, 1930, 50)
my_cindices <- rep(0, length(gene_sizes))
counter <- 1

for(gs in gene_sizes){
  current_cox <- cox_model_fitter(my.seed = 1, my.alpha = 1, cox.df = cox_df,
                                  gene.num = gs, tumor.m = FALSE,
                                  tumor.n = FALSE, tumor.stage = FALSE,
                                  cox.predictors = deseq2_genes[1:gs],
                                  my.dataset = "COAD",
                                  use.foldids = FALSE,
                                  calc.auc = FALSE, save.coefs = TRUE,
                                  my.filename = "deseq2_coefs_read.csv")
  
  my_cindices[counter] <- round(current_cox$CV$cvm[current_cox$CV$index[1]],
                                digits = 4)
  
  counter <- counter + 1
  
}

print(max(my_cindices))