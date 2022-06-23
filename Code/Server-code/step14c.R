#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")



mad.genes <- readRDS("MAD_outputs/mad_colon_and_rectal_cancer.rds")
sde.genes <- readRDS("SDE_outputs/sde_colon_and_rectal_cancer.rds")
cox_df <- readRDS("coad_df_finished.rds")


#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")


combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
                "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010")

alpha_value <- c(0.0,0.5,1.0)

gene_sizes <- seq(100, 3000, 50)
data_set <- "coad"


for(a in alpha_value[3]){
  top_cindices <- c()
  for(c in combo_used[11]){
    my_cindices <- rep(0, 11)
    counter <- 1
    mad_sde_optimized <- two_weight_optimizer(first.metric = mad.genes,
                                              second.metric = sde.genes,
                                              my.filename = paste0("Optimization_",c,"_targets_mad_sde_",data_set,".rds"))
    
    
    
    index <- 1
    
    for(md in mad_sde_optimized[1:11]){
      for(gs in gene_sizes){
        cox_model <- cox_model_fitter(my.seed = 1,
                                      my.alpha = a,
                                      cox.predictors = md,
                                      cox.df = cox_df,
                                      gene.num = gs,
                                      tumor.stage = FALSE,
                                      tumor.n = FALSE,
                                      tumor.m = FALSE,
                                      my.filename = paste0("mad_sde_coad_coefs.csv"))
        
        #Getting the top concordance index from the cross validation and then rounding
        #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
        #the c_index list with the result
        current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
        my_cindices[counter] <- current_cindex
        counter <- counter + 1
        index <- index+1
        
      }
      
      top_cindex <-max(my_cindices)
      top_index <- which(my_cindices==top_cindex)
      print(top_index)
      print(top_cindex)
      top_index_used <- top_index[1]
      top_cindices <- c(top_cindices, top_cindex)
      index <- index+1
    }
    
    write.csv(my_cindices, file = paste0("top_cindices_mad_sde_coad_df.csv"))
  }
  
}