#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")

des_results <- readRDS(file = "des_results_read.rds")
des_results <- filter(des_results, pvalue.adj.FDR < 0.05)

gene_sizes <- seq(100, 3000, 50)
# cox_df <-readRDS("coad_df_finished.rds")
cox_df <-readRDS("read_df_finished.rds")
my_cindices <- rep(0, length(gene_sizes))
counter <- 1

for(gs in gene_sizes){
  current_cox <- cox_model_fitter(my.seed = 1, my.alpha = 1, cox.df = cox_df,
                                  gene.num = gs, tumor.m = FALSE,
                                  tumor.n = FALSE, tumor.stage = FALSE,
                                  cox.predictors = rownames(des_results)[1:gs],
                                  my.dataset = "READ",
                                  save.coefs = TRUE,
                                  calc.auc = FALSE,
                                  use.foldids = FALSE,
                                  my.filename = "desingle_coefs_read.csv")
  
  my_cindices[counter] <- round(current_cox$CV$cvm[current_cox$CV$index[1]],
                                digits = 4)
  
  counter <- counter + 1
  
}

print(max(my_cindices))