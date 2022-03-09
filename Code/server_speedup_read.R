#Name: server_speedup_read.R
#Purpose: To exploit the power of the lab servers for READ data

#mad_sde_coad_optimized <- readRDS(file = "mad_sde_coad_optimized.rds")
cox_df <- readRDS(file = "read_df_finished.rds")
#load(file = "800_1010_targets.RData", verbose = TRUE)
#mirna.genes <- mirna.ranking
sde.genes <- readRDS(file = "sde_colon_and_rectal_cancer.rds")
#mad.genes <- readRDS(file = "mad_colon_and_rectal_cancer.rds")

#Getting the needed functions
source("cox_model.R")
source("model_optimizer.R")
# #Elastic-net penalized cox model----
# c_index <- rep(0, 59)
# gene_sizes <- seq(100, 3000, 50)
# 
# for(gs in gene_sizes){
#   for(ms in mad_sde_coad_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = ms,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("mad_sde_coad_coefs_",gs,"_genes_",ms,"_index.csv"))
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
# write.csv(mad_sde_finished_df, "mad_sde_coad_df.csv")

#CC Singlecell MS
# mirna_sde_coad_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                                second.metric = sde.genes,
#                                                my.filename = "cc_singlecell_ms_coad_optimized.rds")
# 
# 
# c_index <- rep(0, 59)
# gene_sizes <- seq(100, 3000, 50)
# 
# for(gs in gene_sizes){
#   for(ms in mirna_sde_coad_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = ms,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("cc_singlecell_ms_coad_coefs_",gs,"_genes_",ms,"_index.csv"))
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
# write.csv(cc_singlecell_ms_finished_df, "cc_singlecell_ms_coad_df.csv")
# 
# 


#CC Singlecell MM
# mirna_mad_coad_optimized <- two_weight_optimizer(first.metric = mirna.genes,
#                                                  second.metric = mad.genes,
#                                                  my.filename = "cc_singlecell_mm_coad_optimized.rds")
# 
# 
# c_index <- rep(0, 59)
# gene_sizes <- seq(100, 3000, 50)
# 
# for(gs in gene_sizes){
#   for(mm in mirna_mad_coad_optimized[1:11]){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.predictors = mm,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   my.filename = paste0("cc_singlecell_mm_coad_coefs_",gs,"_genes_",mm,"_index.csv"))
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
# write.csv(cc_singlecell_mm_finished_df, "cc_singlecell_mm_coad_df.csv")


#Grid search TCGA-READ----
#CC Singlecell MS grid search
combo_used <- c("800_10", "800_110", "800_210", "800_310", "800_410", "800_510", "800_610", "800_710", "800_810", "800_910", "800_1010")
mirna_used <- gsub(combo_used[1], pattern = "_10", replacement = "")
counter <- 1
c_index <- rep(0, 121)
data_set <- "read"

for(c in combo_used){
  load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/Mirna/Inputs/",c,"_targets.RData"), verbose = TRUE)
  mirna.genes <- mirna.ranking
  mirna_sde_optimized <- two_weight_optimizer(first.metric = mirna.genes,
                                              second.metric = sde.genes,
                                              my.filename = paste0("Optimization_",c,"_targets_cc_singlecell_ms_",data_set,".rds"))
  
  
  
  
  
  
  
  for(ms in mirna_sde_optimized[1:11]){
    cox_model <- cox_model_fitter(my.seed = 1,
                                  cox.predictors = ms,
                                  cox.df = cox_df,
                                  gene.num = 350,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  my.filename = paste0("cc_singlecell_ms_",data_set,"_coefs_",gs,"_genes_",ms,"_index.csv"))
    
    #Getting the top concordance index from the cross validation and then rounding
    #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
    #the c_index list with the result
    top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
    c_index[counter] <- top_cindex
    counter <- counter + 1
  }
  
  
}

write.csv(c_index, paste0("cc_singlecell_ms_grid_search_",data_set,"_df_",mirna_used,".csv"))










