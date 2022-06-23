#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")



#For finding the optimal number of genes
gene_sizes <- seq(100, 3000, 50)
scdd_res <- readRDS("scdd_res_read.rds")
#cox_df <-readRDS("coad_df_finished.rds")
cox_df <-readRDS("read_df_finished.rds")
my_cindices <- rep(0, length(gene_sizes))
counter <- 1

for(gs in gene_sizes){
  current_cox <- cox_model_fitter(my.seed = 1, my.alpha = 1, cox.df = cox_df,
                                  gene.num = gs, tumor.m = FALSE,
                                  tumor.n = FALSE, tumor.stage = FALSE,
                                  cox.predictors = scdd_res$gene[1:gs],
                                  my.dataset = "READ",
                                  use.foldids = FALSE,
                                  calc.auc = FALSE,
                                  save.coefs = TRUE,
                                  my.filename = "scdd_coefs_read.csv")
  
  my_cindices[counter] <- round(current_cox$CV$cvm[current_cox$CV$index[1]],
                                digits = 4)
  
  counter <- counter + 1
  
}

print(max(my_cindices))