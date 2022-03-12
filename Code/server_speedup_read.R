#Name: server_speedup_read.R
#Purpose: To exploit the power of the lab servers for READ dataset

#mad_sde_read_optimized <- readRDS(file = "mad_sde_read_optimized.rds")
cox_df <- readRDS(file = "read_df_finished.rds")
#load(file = "800_1010_targets.RData", verbose = TRUE)
#mirna.genes <- mirna.ranking
#sde.genes <- readRDS(file = "sde_colon_and_rectal_cancer.rds")
mad.genes <- readRDS(file = "mad_colon_and_rectal_cancer.rds")


#Getting the needed functions
source("cox_model.R")
source("model_optimizer.R")
# #Elastic-net penalized cox model----
# c_index <- rep(0, 59)
# gene_sizes <- seq(100, 3000, 50)
# 
# for(gs in gene_sizes){
#   for(ms in mad_sde_read_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = ms,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("mad_sde_read_coefs_",gs,"_genes_",ms,"_index.csv"))
#     
#     #Getting the top concordance index from the cross validation and then rounding
#     #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#     #the c_index list with the result
#     top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#     c_index[which(gene_sizes==gs)] <- top_cindex
#     
#   }
#   
#   
# }
# 
# mad_sde_finished_df <- as.data.frame(cbind(gene_sizes, c_index))
# 
# write.csv(mad_sde_finished_df, "mad_sde_read_df.csv")

#CC Singlecell MS
# mirna_sde_read_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                                second.metric = sde.genes,
#                                                my.filename = "cc_singlecell_ms_read_optimized.rds")
# 
# 
# c_index <- rep(0, 59)
# gene_sizes <- seq(100, 3000, 50)
# 
# for(gs in gene_sizes){
#   for(ms in mirna_sde_read_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = ms,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("cc_singlecell_ms_read_coefs_",gs,"_genes_",ms,"_index.csv"))
# 
#     #Getting the top concordance index from the cross validation and then rounding
#     #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#     #the c_index list with the result
#     top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#     c_index[which(gene_sizes==gs)] <- top_cindex
# 
#   }
# 
# 
# }
# 
# cc_singlecell_ms_finished_df <- as.data.frame(cbind(gene_sizes, c_index))
# 
# write.csv(cc_singlecell_ms_finished_df, "cc_singlecell_ms_read_df.csv")
# 
# 


#CC Singlecell MM
# mirna_mad_read_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                                  second.metric = mad.genes,
#                                                  my.filename = "cc_singlecell_mm_read_optimized.rds")
# 
# 
# c_index <- rep(0, 59)
# gene_sizes <- seq(100, 3000, 50)
# 
# for(gs in gene_sizes){
#   for(mm in mirna_mad_read_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = mm,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("cc_singlecell_mm_read_coefs_",gs,"_genes_",mm,"_index.csv"))
#     
#     #Getting the top concordance index from the cross validation and then rounding
#     #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#     #the c_index list with the result
#     top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#     c_index[which(gene_sizes==gs)] <- top_cindex
#     
#   }
#   
#   
# }
# 
# cc_singlecell_mm_finished_df <- as.data.frame(cbind(gene_sizes, c_index))
# 
# write.csv(cc_singlecell_mm_finished_df, "cc_singlecell_mm_read_df.csv")


#Grid search TCGA-READ----
#CC Singlecell MS grid search
# combo_used <- c("1000_10", "1000_110", "1000_210", "1000_310", "1000_410", "1000_510", "1000_610", "1000_710", "1000_810", "1000_910", "1000_1010",
#                 "900_10", "900_110", "900_210", "900_310", "900_410", "900_510", "900_610", "900_710", "900_810", "900_910", "900_1010",
#                 "800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
#                 "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
#                 "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
#                 "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
#                 "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
#                 "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
#                 "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
#                 "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")

# mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
# data_set <- "read"
# top_cindices <- c()
# 
# for(c in combo_used){
#   my_cindices <- rep(0, 11)
#   counter <- 1
#   load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
#   mirna.genes <- mirna.ranking
#   mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                               second.metric = sde.genes,
#                                               my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_ms_",data_set,".rds"))
#   
#   
#   
#   
#   
#   
#   
#   for(ms in mirna_sde_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = ms,
#                                   cox.df = cox_df,
#                                   gene.num = 350,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("cc_singlecell_ms_",data_set,"_coefs_",gs,"_genes_",ms,"_index.csv"))
#     
#     #Getting the top concordance index from the cross validation and then rounding
#     #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
#     #the c_index list with the result
#     current_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
#     my_cindices[counter] <- current_cindex
#     counter <- counter + 1
#     
#   }
#   
#   top_cindex <-max(my_cindices)
#   top_index <- which(my_cindices==top_cindex)
#   print(top_index)
#   print(top_cindex)
#   top_index_used <- top_index[1]
#   print(my_cindices[top_index_used])
#   #print(top_index_used)
#   
#   top_cindices <- c(top_cindices, top_cindex)
#   
#   
#   
#   
# }
# 
# 
# write.csv(top_cindices, file = paste0("top_cindices_mirna_used_combo_",c,"_index_",top_index_used,"_df.csv"))



#CC Singlecell MM grid search
combo_used <- c("1000_10", "1000_110", "1000_210", "1000_310", "1000_410", "1000_510", "1000_610", "1000_710", "1000_810", "1000_910", "1000_1010",
                "900_10", "900_110", "900_210", "900_310", "900_410", "900_510", "900_610", "900_710", "900_810", "900_910", "900_1010",
                "800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010",
                "700_10", "700_110", "700_210", "700_310", "700_410", "700_510", "700_610", "700_710", "700_810", "700_910", "700_1010",
                "600_10", "600_110", "600_210", "600_310", "600_410", "600_510", "600_610", "600_710", "600_810", "600_910", "600_1010",
                "500_10", "500_110", "500_210", "500_310", "500_410", "500_510", "500_610", "500_710", "500_810", "500_910", "500_1010",
                "400_10", "400_110", "400_210", "400_310", "400_410", "400_510", "400_610", "400_710", "400_810", "400_910", "400_1010",
                "300_10", "300_110", "300_210", "300_310", "300_410", "300_510", "300_610", "300_710", "300_810", "300_910", "300_1010",
                "200_10", "200_110", "200_210", "200_310", "200_410", "200_510", "200_610", "200_710", "200_810", "200_910", "200_1010",
                "100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")




mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
data_set <- "read"
top_cindices <- c()

for(c in combo_used){
  my_cindices <- rep(0, 11)
  counter <- 1
  load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
  mirna.genes <- mirna.ranking
  mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
                                              second.metric = mad.genes,
                                              my.filename = paste0("Optimizations/Optimization_",c,"_targets_cc_singlecell_mm_0_alpha_",data_set,".rds"))
  
  
  
  
  
  
  
  for(mm in mirna_sde_optimized[1:11]){
    cox_model <- cox_model_fitter(my.seed = 1,
                                  cox.predictors = mm,
                                  cox.df = cox_df,
                                  gene.num = 350,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  my.alpha = 0,
                                  my.dataset = "READ",
                                  my.filename = paste0("cc_singlecell_mm_read_alpha_0",data_set,"_coefs_",gs,"_genes_",mm,"_index.csv"))
    
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
  print(my_cindices[top_index_used])
  
  top_cindices <- c(top_cindices, top_cindex)
  
  
  
  
}


write.csv(top_cindices, file = paste0("top_cindices_alpha_0_mirna_mad_read_used_combo_",c,"_index_",top_index_used,"_df.csv"))



