#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files. Originally the best combo_used was 27
source("cox_model.R")
source("model_optimizer.R")



# combo_used <- c("1000_10", "1000_110", "1000_210", "1000_310", "1000_410", "1000_510", "1000_610", "1000_710", "1000_810", "1000_910", "1000_1010",
#                 "900_10", "900_110", "900_210", "900_310", "900_410", "900_510", "900_610", "900_710", "900_810", "900_910", "900_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010")

combo_used <- "300_310"

sde.genes <- readRDS("SDE_outputs/sde_colon_and_rectal_cancer.rds")
cox_df <- readRDS("read_df_finished.rds")

alpha_value <- c(0.0,0.5,1.0)

mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
data_set <- "read"


for(a in alpha_value[3]){
  top_cindices <- c()
  for(c in combo_used){
    my_cindices <- rep(0, 11)
    counter <- 1
    load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/miRNA_inputs/Inputs/",c,"_targets.RData"), verbose = TRUE)
    mirna.genes <- mirna.ranking
    mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
                                                second.metric = sde.genes,
                                                my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_ms_compared_to_mms_",a,"_alpha_",data_set,".rds"))
    
    
    
    
    
    
    
    for(ms in mirna_sde_optimized){
      cox_model <- cox_model_fitter(my.seed = 1,
                                    my.alpha = a,
                                    my.dataset = "READ",
                                    cox.predictors = ms,
                                    cox.df = cox_df,
                                    gene.num = 600,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    my.filename = paste0("cc_singlecell_ms_ideal_mirna_combo_coefs_compared_to_mms.csv"))
      
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
  
  write.csv(top_cindices, file = paste0("cc_singlecell_ms_read_df_ideal_mirna_combo_compared_to_mms.csv"))
  write.csv(my_cindices, file = paste0("cc_singlecell_ms_read_df_all_mirna_combo_performance_data_compared_to_mms.csv"))
}