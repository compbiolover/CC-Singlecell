#Name: server_speedup_coad.R
#Purpose: To use the power of the lab servers for COAD data set

#mad_sde_coad_optimized <- readRDS(file = "mad_sde_coad_optimized.rds")
cox_df <- readRDS(file = "coad_df_finished.rds")
#load(file = "800_1010_targets.RData", verbose = TRUE)
#mirna.genes <- mirna.ranking
#sde.genes <- readRDS(file = "sde_colon_and_rectal_cancer.rds")
#mad.genes <- readRDS(file = "mad_colon_and_rectal_cancer.rds")

#Getting the needed functions
source("cox_model.R")
source("model_optimizer.R")


#Ideal miRNA + SDE gene number----
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
#                 "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
#                 "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
#                 "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
#                 "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
#                 "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")
# 
# alpha_value <- c(0.0,0.5,1.0)
# #alpha_value <- c(0.6,0.7,0.8,0.9)
# #alpha_value <- c(0.1,0.2,0.3,0.4)
# gene_sizes <- seq(100,3000, 100)
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# # top_cindices <- c()
# 
# for(a in alpha_value){
#   for(gs in gene_sizes){
#     top_cindices <- c()
#     for(c in combo_used){
#       my_cindices <- rep(0, 11)
#       counter <- 1
#       load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
#       mirna.genes <- mirna.ranking
#       mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                                   second.metric = sde.genes,
#                                                   my.filename = paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_ms_",a,"_alpha_",data_set,"_ideal_gene_num_try2_",gs,".rds"))
#       
#       
#       
#       
#       
#       
#       
#       index <- 1
#       
#       for(ms in mirna_sde_optimized){
#         cox_model <- cox_model_fitter(my.seed = 1,
#                                       my.alpha = a,
#                                       my.dataset = "COAD",
#                                       cox.predictors = ms,
#                                       cox.df = cox_df,
#                                       gene.num = 900,
#                                       tumor.stage = FALSE,
#                                       tumor.n = FALSE,
#                                       tumor.m = FALSE,
#                                       my.filename = paste0("cc_singlecell_ms_",data_set,"_alpha_",a,"_coefs_",index,"_index_ideal_mirna_sde_gene_num_try2_",gs,".csv"))
#         
#         #Getting the top concordance index from the cross validation and then rounding
#         #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#         #the c_index list with the result
#         current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#         my_cindices[counter] <- current_cindex
#         counter <- counter + 1
#         index <- index+1
#         
#       }
#       
#       #print(cox_model$CV)
#       top_cindex <-max(my_cindices)
#       top_index <- which(my_cindices==top_cindex)
#       print(top_index)
#       print(top_cindex)
#       top_index_used <- top_index[1]
#       #print(my_cindices[top_index_used])
#       #print(top_index_used)
#       top_cindices <- c(top_cindices, top_cindex)
#       index <- index+1
#     }
#     
#     write.csv(top_cindices, file = paste0("top_cidices_mirna_sde_try2_",gs,".csv"))
#     #mirna_num <- rep(seq(800,100,-100), each=11)
#     #mirna_targets <- seq(10,1010,100)
#     #top_cindices_df <- data.frame(mirna_num, mirna_targets, top_cindices)
#     #write.csv(top_cindices_df, file = paste0("top_cindices_alpha_",a,"_cc_singlecell_ms_coad_df_finished.csv"))
#     
#     
#     
#     
#   }
#   
# 
# }
# 
# 


#Grid search TCGA-COAD----
#miRNA only grid search
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
#                 "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
#                 "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
#                 "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
#                 "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
#                 "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")
# 
# alpha_value <- c(0.0,0.5,1.0)
# 
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# 
# for(a in alpha_value[3]){
#   my_cindices <- c()
#   counter <- 1
#   for(c in combo_used){
#     load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
#     mirna.genes <- mirna.ranking
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   my.alpha = a,
#                                   cox.predictors = mirna.genes,
#                                   cox.df = cox_df,
#                                   gene.num = 1900,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("cc_singlecell_ms_alpha_",a,"_coefs_",index,"_index_finished_mirna_only_grid_search_finished.csv"))
# 
#     #Getting the top concordance index from the cross validation and then rounding
#     #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#     #the c_index list with the result
#     current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]],
#                             digits = 4)
#     my_cindices[counter] <- current_cindex
#     counter <- counter + 1
# 
#   }
#   write.csv(my_cindices, file = paste0("alpha_",a,"_cindices_mirna_only_finished.csv"))
# 
# }
# 






#CC Singlecell MS grid search
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
#                 "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
#                 "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
#                 "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
#                 "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
#                 "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")
# 
# alpha_value <- c(0.0,0.5,1.0)
# #alpha_value <- c(0.6,0.7,0.8,0.9)
# #alpha_value <- c(0.1,0.2,0.3,0.4)
# gene_sizes <- seq(100, 3000, 50)
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# # top_cindices <- c()
# 
# for(a in alpha_value[3]){
#   top_cindices <- c()
#   for(c in combo_used){
#     my_cindices <- rep(0, 11)
#     counter <- 1
#     # load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
#     # mirna.genes <- mirna.ranking
#     # mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#     #                                             second.metric = sde.genes,
#     #                                             my.filename = paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_ms_",a,"_alpha_",data_set,".rds"))
# 
#     mirna_sde_optimized <- readRDS(paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_ms_1_alpha_",data_set,".rds"))
# 
# 
# 
# 
# 
# 
#     index <- 1
# 
#     for(ms in mirna_sde_optimized){
#       cox_model <- cox_model_fitter(my.seed = 1,
#                                     my.alpha = a,
#                                     cox.predictors = ms,
#                                     cox.df = cox_df,
#                                     gene.num = 1900,
#                                     tumor.stage = FALSE,
#                                     tumor.n = FALSE,
#                                     tumor.m = FALSE,
#                                     my.filename = paste0("cc_singlecell_ms_coad_alpha_",a,"_coefs_",index,"_index_finished_no_foldid_all_combos.csv"))
# 
#       #Getting the top concordance index from the cross validation and then rounding
#       #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#       #the c_index list with the result
#       current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#       my_cindices[counter] <- current_cindex
#       counter <- counter + 1
#       index <- index+1
# 
#     }
# 
#     #print(cox_model$CV)
#     top_cindex <-max(my_cindices)
#     top_index <- which(my_cindices==top_cindex)
#     print(top_index)
#     print(top_cindex)
#     top_index_used <- top_index[1]
#     #print(my_cindices[top_index_used])
#     #print(top_index_used)
#     top_cindices <- c(top_cindices, top_cindex)
#     index <- index+1
#   }
# 
#   mirna_num <- rep(seq(800,100,-100), each=11)
#   mirna_targets <- rep(seq(10,1010,100), times=8)
#   top_cindices_df <- data.frame(mirna_num, mirna_targets, top_cindices)
#   write.csv(top_cindices_df, file = paste0("top_cindices_alpha_",a,"_cc_singlecell_ms_coad_df_finished_no_foldid_finished_all_combos.csv"))
# 
# 
# 
# }
# 
# 

#For optimal gene size of CC Singlecell MS
combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
                "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
                "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
                "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
                "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
                "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
                "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
                "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")



gene_sizes <- seq(100, 3000, 50)
top_cindices <- c()
combo_progress <- 1

for(c in combo_used){
  gene_progress <- 1
  mirna_sde_optimized <- readRDS(paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_ms_1_alpha_coad.rds"))
  print(paste0("We are currently ",round(combo_progress/length(combo_used), digits = 4)*100,"% done with all combinations"))
  combo_progress <- combo_progress + 1
  for(gs in gene_sizes){
    my_cindices <- rep(0, 11)
    counter <- 1
    print(paste0("For this current combo we have searched ",round(gene_progress/length(gene_sizes), digits = 4)*100,"% of all gene sizes"))
    gene_progress <- gene_progress + 1
    for(ms in mirna_sde_optimized){
      cox_model <- cox_model_fitter(my.seed = 1,
                                    my.alpha = 1,
                                    cox.predictors = ms,
                                    cox.df = cox_df,
                                    gene.num = gs,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    my.filename = paste0("cc_singlecell_ms_coad_alpha_1_coefs_gene_size_",gs,"_combo_",c,".csv"))
      
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
    top_cindices <- c(top_cindices, top_cindex)
  }
  
  gene_size_df <- data.frame(gene_sizes, top_cindices)
  colnames(gene_size_df) <- c("gene_size", "c_index")
  write.csv(gene_size_df, file = "cc_singlecell_ms_optimal_gene_df.csv")
  
}







#CC Singlecell MM grid search
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
#                 "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
#                 "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
#                 "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
#                 "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
#                 "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")
# 
# alpha_value <- c(0.0,0.5,1.0)
# 
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# # top_cindices <- c()
# 
# for(a in alpha_value){
#   top_cindices <- c()
#   for(c in combo_used){
#     my_cindices <- rep(0, 11)
#     counter <- 1
#     load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
#     mirna.genes <- mirna.ranking
#     mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                                 second.metric = mad.genes,
#                                                 my.filename = paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_mm_",a,"_alpha_",data_set,".rds"))
# 
# 
# 
# 
# 
# 
# 
#     for(mm in mirna_sde_optimized[1:11]){
#       cox_model <- cox_model_fitter(my.seed = 1,
#                                     my.alpha = a,
#                                     my.dataset = "COAD",
#                                     cox.predictors = mm,
#                                     cox.df = cox_df,
#                                     gene.num = 1900,
#                                     tumor.stage = FALSE,
#                                     tumor.n = FALSE,
#                                     tumor.m = FALSE,
#                                     my.filename = paste0("cc_singlecell_mm_",data_set,"_alpha_",a,"_coefs_",gs,"_genes_",mm,"_index.csv"))
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
#     #print(my_cindices[top_index_used])
#     #print(top_index_used)
#     top_cindices <- c(top_cindices, top_cindex)
#   }
# 
#   write.csv(top_cindices, file = paste0("top_cindices_alpha_",a,"_cc_singlecell_mm_coad_used_combo_",c,"_index_",top_index_used,"_df.csv"))
# 
# 
# 
# }
# 
# 
# #CC Singlecell MMS
# combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
#                 "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
#                 "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
#                 "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
#                 "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
#                 "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")
# 
# 
# 
# 
# alpha_value <- c(0.0,0.5,1.0)
# 
# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "coad"
# # top_cindices <- c()
# 
# for(a in alpha_value){
#   top_cindices <- c()
#   for(c in combo_used){
#     my_cindices <- rep(0, 121)
#     counter <- 1
#     load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
#     mirna.genes <- mirna.ranking
#     mirna_mad_sde_optimized <- three_weight_optimizer(first.metric = mirna.genes,
#                                                 second.metric = mad.genes,
#                                                 third.metric = sde.genes,
#                                                 my.filename = paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_mms_",a,"_alpha_",data_set,".rds"))
# 
# 
# 
# 
# 
# 
# 
#     for(mms in mirna_mad_sde_optimized[1:121]){
#       cox_model <- cox_model_fitter(my.seed = 1,
#                                     my.alpha = a,
#                                     my.dataset = "COAD",
#                                     cox.predictors = mms,
#                                     cox.df = cox_df,
#                                     gene.num = 1900,
#                                     tumor.stage = FALSE,
#                                     tumor.n = FALSE,
#                                     tumor.m = FALSE,
#                                     my.filename = paste0("cc_singlecell_mms_",data_set,"_alpha_",a,"_coefs_",gs,"_genes_",mms,"_index.csv"))
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
#     #print(my_cindices[top_index_used])
#     #print(top_index_used)
#     top_cindices <- c(top_cindices, top_cindex)
#   }
# 
#   write.csv(top_cindices, file = paste0("top_cindices_alpha_",a,"_cc_singlecell_mms_coad_used_combo_index_df.csv"))
# 
# 
# 
# }
# 
