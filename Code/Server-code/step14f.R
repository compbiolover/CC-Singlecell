#Loading needed packages
library(BiocParallel)
library(tidyverse)


#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")

##Ideal CC Singlecell MM combo used
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010")
# 
# mad.genes <- readRDS("MAD_outputs/mad_colon_and_rectal_cancer.rds")
# sde.genes <- readRDS("SDE_outputs/sde_colon_and_rectal_cancer.rds")
# cox_df <- readRDS("coad_df_finished.rds")
# 
# alpha_value <- c(0.0,0.5,1.0)
# 
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# 
# 
# for(a in alpha_value[3]){
#   top_cindices <- c()
#   for(c in combo_used[19]){
#     my_cindices <- rep(0, 11)
#     counter <- 1
#     load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/miRNA_inputs/Inputs/",c,"_targets.RData"), verbose = TRUE)
#     mirna.genes <- mirna.ranking
#     mirna_mad_sde_optimized <- three_weight_optimizer(first.metric = mirna.genes,
#                                                       second.metric = mad.genes,
#                                                       third.metric = sde.genes,
#                                                       my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_mm_",a,"_alpha_",data_set,".rds"))
#     
#     
#     
#     
#     
#     
#     
#     for(mms in mirna_mad_sde_optimized){
#       cox_model <- cox_model_fitter(my.seed = 1,
#                                     my.alpha = a,
#                                     my.dataset = "COAD",
#                                     cox.predictors = mms,
#                                     cox.df = cox_df,
#                                     gene.num = 2500,
#                                     tumor.stage = FALSE,
#                                     tumor.n = FALSE,
#                                     tumor.m = FALSE,
#                                     my.filename = paste0("cc_singlecell_mms_ideal_mirna_combo_coefs.csv"))
#       
#       #Getting the top concordance index from the cross validation and then rounding
#       #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#       #the c_index list with the result
#       current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#       my_cindices[counter] <- current_cindex
#       counter <- counter + 1
#       
#     }
#     
#     top_cindex <-max(my_cindices)
#     top_index <- which(my_cindices==top_cindex)
#     print(top_index)
#     print(top_cindex)
#     top_index_used <- top_index[1]
#     top_cindices <- c(top_cindices, top_cindex)
#   }
#   
#   write.csv(top_cindices, file = paste0("cc_singlecell_mms_coad_df_ideal_mirna_combo_mm_ideal_combo.csv"))
#   write.csv(my_cindices, file = paste0("cc_singlecell_mms_coad_df_all_mirna_combo_performance_data_mm_ideal_combo.csv"))
# }


#Ideal CC Singelcell MS combo used for CC Singlecell MMS
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010")
# 
# mad.genes <- readRDS("MAD_outputs/mad_colon_and_rectal_cancer.rds")
# sde.genes <- readRDS("SDE_outputs/sde_colon_and_rectal_cancer.rds")
# cox_df <- readRDS("coad_df_finished.rds")
# 
# alpha_value <- c(0.0,0.5,1.0)
# 
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# 
# 
# for(a in alpha_value[3]){
#   top_cindices <- c()
#   for(c in combo_used[11]){
#     my_cindices <- rep(0, 121)
#     counter <- 1
#     load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/miRNA_inputs/Inputs/",c,"_targets.RData"), verbose = TRUE)
#     mirna.genes <- mirna.ranking
#     mirna_mad_sde_optimized <- three_weight_optimizer(first.metric = mirna.genes,
#                                                       second.metric = mad.genes,
#                                                       third.metric = sde.genes,
#                                                       my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_mms_ms_ideal_used",a,"_alpha_",data_set,".rds"))
#     
#     
#     
#     
#     
#     
#     
#     for(mms in mirna_mad_sde_optimized[112]){
#       cox_model <- cox_model_fitter(my.seed = 1,
#                                     my.alpha = a,
#                                     my.dataset = "COAD",
#                                     cox.predictors = mms,
#                                     cox.df = cox_df,
#                                     gene.num = 2500,
#                                     tumor.stage = FALSE,
#                                     tumor.n = FALSE,
#                                     tumor.m = FALSE,
#                                     my.filename = paste0("cc_singlecell_mms_ideal_mirna_combo_coefs_of_top_performer_ms.csv"))
#       
#       #Getting the top concordance index from the cross validation and then rounding
#       #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#       #the c_index list with the result
#       current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#       my_cindices[counter] <- current_cindex
#       counter <- counter + 1
#       
#     }
#     
#     top_cindex <-max(my_cindices)
#     top_index <- which(my_cindices==top_cindex)
#     print(top_index)
#     print(top_cindex)
#     top_index_used <- top_index[1]
#     top_cindices <- c(top_cindices, top_cindex)
#   }
#   
#   #write.csv(top_cindices, file = paste0("cc_singlecell_mms_coad_df_ideal_mirna_combo_ms_combo_used.csv"))
#   #write.csv(my_cindices, file = paste0("cc_singlecell_mms_coad_df_all_mirna_combo_performance_data_ms_combo_used.csv"))
# }
# 

#Now doing the entire grid search at the ideal gene size of CC Singlecell MS for CC Singlecell MMS. It occurs at the 20th combo used with 102 weighting of the 3 metric model
combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
                "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
                "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
                "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
                "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
                "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
                "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
                "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")

mad.genes <- readRDS("MAD_outputs/mad_colon_and_rectal_cancer.rds")
sde.genes <- readRDS("SDE_outputs/sde_colon_and_rectal_cancer.rds")
cox_df <- readRDS("coad_df_finished.rds")

alpha_value <- c(0.0,0.5,1.0)

mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
data_set <- "coad"


for(a in alpha_value[3]){
  top_cindices <- c()
  for(c in combo_used[20]){
    my_cindices <- rep(0, 121)
    counter <- 1
    load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/miRNA_inputs/Inputs/",c,"_targets.RData"), verbose = TRUE)
    mirna.genes <- mirna.ranking
    # mirna_mad_sde_optimized <- three_weight_optimizer(first.metric = mirna.genes,
    #                                                   second.metric = mad.genes,
    #                                                   third.metric = sde.genes,
    #                                                   my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_mms_ms_ideal_used_",a,"_alpha_",data_set,".rds"))
    # 
    
    mirna_mad_sde_optimized <- readRDS(file = "Optimization_700_810_targets_cc_singlecell_mms_ms_ideal_used1_alpha_coad.rds")
    
    
    
    
    for(mms in mirna_mad_sde_optimized[102]){
      cox_model <- cox_model_fitter(my.seed = 1,
                                    my.alpha = a,
                                    my.dataset = "COAD",
                                    cox.predictors = mms,
                                    cox.df = cox_df,
                                    gene.num = 2500,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    my.filename = paste0("cc_singlecell_mms_ideal_mirna_fig_s3_coefs.csv"))
      
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
  
  write.csv(top_cindices, file = paste0("cc_singlecell_mms_coad_top_cindex_tested.csv"))
}