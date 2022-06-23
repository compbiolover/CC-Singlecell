#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")

cox_df <- readRDS("read_df_finished.rds")

all_random_gene_cindices <- seq(1,10, 1)
random_read_gene_list <- readRDS("random_read_genes.rds")

for(rg in all_random_gene_cindices){
  current_genes <- random_read_gene_list[[rg]]
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = current_genes,
                                cox.df = cox_df,
                                gene.num = 350,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("random_read_coefs_for_random_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the all_random_gene_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  all_random_gene_cindices[rg] <- top_cindex
  
}

write.csv(all_random_gene_cindices, file = "read_random_genes.csv")

mean_random <- mean(all_random_gene_cindices)

write.csv(mean_random, file = "read_mean_random_genes.csv")