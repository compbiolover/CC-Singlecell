#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")



gene_sizes <- seq(100, 3000, 50)
mirna_high_cindices <- rep(0, length(gene_sizes))

#Loading the high miRNA-miRNA target file
load("800_1010_targets.RData",
     verbose = TRUE)

high.mirna.genes <- mirna.ranking
cox_df <- readRDS("coad_df_finished.rds")

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = high.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("high_mirna_coad_coefs_",gs,"_genes.csv"))
  
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
          "mirna_high_cindices_coad_across_gene_size.csv")