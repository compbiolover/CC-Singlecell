#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")


combo_used <- c("1000_10", "1000_110", "1000_210", "1000_310", "1000_410", "1000_510", "1000_610", "1000_710", "1000_810", "1000_910", "1000_1010",
                "900_10", "900_110", "900_210", "900_310", "900_410", "900_510", "900_610", "900_710", "900_810", "900_910", "900_1010")

gene_sizes <- seq(100, 3000, 50)
mad.genes <- readRDS("MAD_outputs/mad_colon_and_rectal_cancer.rds")
cox_df <- readRDS("read_df_finished.rds")

alpha_value <- c(0.0,0.5,1.0)

mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
data_set <- "read"


for(a in alpha_value[3]){
  top_cindices <- c()
  for(c in combo_used[2]){
    my_cindices <- rep(0, 11)
    counter <- 1
    load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/miRNA_inputs/Inputs/",c,"_targets.RData"), verbose = TRUE)
    mirna.genes <- mirna.ranking
    mirna_mad_optimized <- two_weight_optimizer(first.metric = mirna.genes,
                                                second.metric = mad.genes,
                                                my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_mm_",a,"_alpha_",data_set,".rds"))
    
    
    
    
    
    
    
    for(mm in mirna_mad_optimized[1:11]){
      for(gs in gene_sizes){
        cox_model <- cox_model_fitter(my.seed = 1,
                                      my.alpha = a,
                                      my.dataset = "READ",
                                      cox.predictors = mm,
                                      cox.df = cox_df,
                                      gene.num = gs,
                                      tumor.stage = FALSE,
                                      tumor.n = FALSE,
                                      tumor.m = FALSE,
                                      my.filename = paste0("cc_singlecell_mm_ideal_gene_size_coefs.csv"))
        
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
      top_index_used <- top_index[1]
      top_cindices <- c(top_cindices, top_cindex)
    }
    
    write.csv(top_cindices, file = paste0("cc_singlecell_mm_read_df_ideal_gene_size.csv"))
    write.csv(my_cindices, file = paste0("cc_singlecell_mm_read_df_all_gene_size_data.csv"))
  }
  
  
  
}