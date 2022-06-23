#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")


sde.genes <- readRDS("sde_colon_and_rectal_cancer.rds")
cox_df <- readRDS("coad_df_finished.rds")

#Getting ideal gene number for MAD metric on TCGA-COAD
gene_sizes <- seq(100, 3000, 50)
mad_cindices <- rep(0, length(gene_sizes))

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = sde.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("sde_coad_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mad_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mad_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mad_cindices_coad_df <- as.data.frame(cbind(gene_sizes, mad_cindices))
mad_cindices_coad_df$method <- rep("SDE", nrow(mad_cindices_coad_df))
colnames(mad_cindices_coad_df)[2] <- "c_index"


#Saving the output
write.csv(mad_cindices_coad_df,
          "sde_cindices_coad_across_gene_size.csv")