#Loading needed packages----
library(doParallel)
library(glmnet)
library(survival)
library(tidyverse)
#Setting the directory for the lab server----
setwd("/home/awillems/Projects/CC_Singlecell")
#Registering multiple cores----
registerDoParallel(12)
#Loading the .RData files to simplify analysis----
load(file = "Data/Exported_data/R_objects/COA_data_se.RData")
load(file = "Data/Exported_data/R_objects/COA_data_df.RData")
load(file = "Data/Exported_data/R_objects/bulk_rna_df.RData")
load(file = "Data/Exported_data/R_objects/all_tumor_cells_fpkm_denoised_cleaned.RData")
load(file = "Data/Exported_data/R_objects/gene_expression_info.RData")
load(file = "Data/Exported_data/R_objects/integrated-gene-lists-for-two-metrics.RData") #This is the SDES + MiRNA metrics
load(file = "Data/Exported_data/R_objects/all_tumor_cells_fpkm_denoised_df_cleaned.RData")
load(file = "Data/Exported_data/R_objects/gene_lists_to_test_MiRNA_SDES.RData")
load(file = "Data/Exported_data/R_objects/all_intersections_cleaned-MiRNA-SDES.RData")
load(file = "Data/Exported_data/R_objects/df_for_train_test_split.RData")
load(file = "Data/Exported_data/R_objects/merged_df_replaced.RData")

#Cox models----
set.seed(1)
cox_models <- vector(mode = "list", length = length(all_intersections_cleaned))
f_objects <- vector(mode = "list", length = length(all_intersections_cleaned))
lambdas <- vector(mode = "list", length = length(all_intersections_cleaned))
c_indicies <- vector(mode = "list", length = length(all_intersections_cleaned))
all_coefs <- vector(mode = "list", length = length(all_intersections_cleaned))
all_active_indicies <- vector(mode = "list", length = length(all_intersections_cleaned))
all_active_coefs <- vector(mode = "list", length = length(all_intersections_cleaned))
all_formulas <- vector(mode = "list", length = length(all_intersections_cleaned))

#Making the y of the cox model----
my_y <- Surv(time = df_for_train_test_split$days.to.last.follow.up, event = df_for_train_test_split$vital.status)

#The cox model loop----
for (x in seq(1:length(all_intersections_cleaned))){
  current_formula_data <- all_intersections_cleaned[[x]]
  current_formula_data <- as.vector(current_formula_data)
  #For just MAD metric
  # current_formula_data <- rownames(t(mad.ranking.subset.df))
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))
  
  #current_formula_data <- current_formula_data[-107]
  
  my_formula <- paste("~", paste(current_formula_data[1:length(current_formula_data)], collapse = "+"))
  f <- as.formula(my_formula)
  
  
  my_x <- model.matrix(f, df_for_train_test_split)
  cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = 10, type.measure = "C", maxit=100000, family="cox", parallel = TRUE)
  current_lambdas <- cv_fit$lambda
  current_c_indicies <- cv_fit$cvm
  
  
  #Adding the 10-fold cox model and other relevant info to lists
  cox_models[[x]] <- cv_fit
  f_objects[[x]] <- f
  lambdas[[x]] <- current_lambdas
  c_indicies[[x]] <- current_c_indicies
  
  
  #Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family =  "cox", maxit = 100000)
  Coefficients <- coef(fit, s = cv_fit$lambda.min)
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  
  #Adding the relevant data bits to lists for further processing
  all_coefs[[x]] <- Coefficients
  all_active_indicies[[x]] <- Active.Index
  all_active_coefs[[x]] <- Active.Coefficients
  all_formulas[[x]] <- my_formula
}

save(all_formulas, file = "Data/Data-from-pipeline/all_formulas_MiRNA_SDES.RData")
save(cox_models, file = "Data/Data-from-pipeline/cox_models_MiRNA_SDES.RData")
save(f_objects, file = "Data/Data-from-pipeline/f_objects_MiRNA_SDES.RData")
save(lambdas, file = "Data/Data-from-pipeline/lambdas_MiRNA_SDES.RData")
save(c_indicies, file = "Data/Data-from-pipeline/c_indicies_MiRNA_SDES.RData")
save(all_coefs, file = "Data/Data-from-pipeline/all_coefs_MiRNA_SDES.RData")
save(all_active_coefs, file = "Data/Data-from-pipeline/all_active_coefs_MiRNA_SDES.RData")
save(all_active_indicies, file = "Data/Data-from-pipeline/all_active_indicies_MiRNA_SDES.RData")

for(x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  print(current_model)
}
