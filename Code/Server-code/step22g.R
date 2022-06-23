#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")


#Getting ideal gene number for SDE metric on TCGA-GBM
gene_sizes <- seq(100, 3000, 50)
sde_cindices <- rep(0, length(gene_sizes))
sde.genes <- readRDS("sde_gbm_cancer.rds")
load("gbm_df.RData", verbose = TRUE)

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = sde.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("sde_gbm_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mad_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  sde_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
sde_cindices_gbm_df <- as.data.frame(cbind(gene_sizes, sde_cindices))
sde_cindices_gbm_df$method <- rep("SDE", nrow(sde_cindices_gbm_df))
colnames(sde_cindices_gbm_df)[2] <- "c_index"

write.csv(sde_cindices_gbm_df, file = "sde_cindices_gbm.csv")

