#Author: Andrew Willems <awillems@vols.utk.edu>.
#Purpose: To take scRNA-seq colon cancer data and begin analysis.
#R version: 4.0.2.

#Loading needed packages----
library(doParallel)
library(glmnet)
library(survival)
library(tidyverse)

#Setting the directory for the lab server----
setwd("/home/awillems/Projects/CC_Singlecell")
#Registering multiple cores----
registerDoParallel(12)
#Copying in the functions that are helpful from the Li et al. paper. ----
#Function used to rank my linear model with just two metrics
two_metric_geneRank <- function(ranking1 = NULL,
                               ranking2  = NULL,
                               a1        = 1,
                               a2        = 0){
  gns = c(names(ranking1),names(ranking2))
  gns = unique(gns)
  res = rep(0, length(gns))
  names(res) = gns
  for(i in names(res)) {
    res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2
  }
  res = res[order(res, decreasing = T)]
  res
}

#Function used to rank my linear model with three metrics
three_metric_geneRank <- function(ranking1 = NULL, 
                                  ranking2 = NULL, 
                                  ranking3 = NULL, 
                                  a1       = 1,
                                  a2       = 0,
                                  a3       = 0){
  gns = c(names(ranking1),names(ranking2),names(ranking3))
  gns = unique(gns)
  res = rep(0, length(gns))
  names(res) = gns
  for(i in names(res)) {
    res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2+getRank(ranking3, i)*a3
  }
  res = res[order(res, decreasing = T)]
  res
}

#Helper function for geneRank functions
getRank <- function(ranking = NULL, 
                    gn      = NULL){
  if (gn %in% names(ranking)) {
    return(ranking[gn])
  }
  else return(0.0)
}
#Loading the .RData files to simplify analysis----
load(file = "Data/Exported_data/R_objects/COA_data_se.RData")
load(file = "Data/Exported_data/R_objects/COA_data_df.RData")
load(file = "Data/Exported_data/R_objects/bulk_rna_df.RData")
load(file = "Data/Exported_data/R_objects/all_tumor_cells_fpkm_denoised_cleaned.RData")
load(file = "Data/Exported_data/R_objects/gene_expression_info.RData")
load(file = "Data/Exported_data/R_objects/mad.ranking.RData")
load(file = "Data/Exported_data/R_objects/vim.sdes.ranking.RData")
#load(file = "Data/Exported_data/R_objects/mirna.ranking.1364.mirnas.RData") #For the mirna that come from miRmap + dbDEMC (high) but a larger subset (1364 miRNAs)
#load(file = "Data/Exported_data/R_objects/mirna.ranking.all.three.dbs.demc.low.RData", verbose = TRUE) #For the mirna that comes from 3 mirn dbs (low version of dbDEMC) 
#load(file = "Data/Exported_data/R_objects/mirna.ranking.all.three.dbs.RData", verbose = TRUE) #For the mirna that comes from 3 mirn dbs (high version of dbDEMC)
#load(file = "Data/Exported_data/R_objects/mirna.ranking.RData")
#load(file = "Data/Exported_data/R_objects/integrated-gene-lists-for-two-metrics.RData") #This is the SDES + MiRNA metrics
load(file = "Data/Exported_data/R_objects/all_tumor_cells_fpkm_denoised_df_cleaned.RData")
#load(file = "Data/Exported_data/R_objects/gene_lists_to_test_MiRNA_SDES.RData")
#load(file = "Data/Exported_data/R_objects/all_intersections_cleaned-MiRNA-SDES.RData")
load(file = "Data/Exported_data/R_objects/df_for_train_test_split.RData")
load(file = "Data/Exported_data/R_objects/merged_df_replaced.RData")
#load(file = "Data/Exported_data/R_objects/mad.ranking.subset.df.RData")
#load(file = "Data/Lab-server-data/All-three-metrics/sorted_gene_lists_all_three_metrics.RData", verbose = TRUE)



#Optimization for just 2 metrics----
two_weight_optimizer <- function(my.start        =0,
                                 my.finish       =1,
                                 step.size       =0.1, 
                                 my.index        =1,
                                 my.list.length  =11,
                                 first.metric    =mad.ranking,
                                 second.metric   =vim.sdes.ranking,
                                 my.filename     ="Data/Exported-data/R-objects/two_metric_optimization.RData"){
  
  #Setting up weights, loop-indexer, and the list that will store the results
  weights <- seq(from = my.start, to=my.finish, by=step.size)
  df_index <- my.index
  integrated_gene_lists <- vector(mode = "list", length = my.list.length)
  
  #Doing the grid search
  for (x in weights) {
    current_ranking <- two_metric_geneRank(ranking1 = first.metric, ranking2 = second.metric,  a1=x, a2=1-x)
    current_ranking <- as.data.frame(current_ranking)
    integrated_gene_lists[[as.character(df_index)]] <- current_ranking
    df_index <- df_index + 1
  }
  
  #Saving the output of the grid search to a .RData
  save(integrated_gene_lists, file = my.filename)
  return(integrated_gene_lists)
}

#Optimization for all 3 metrics (MAD, SDE and miRNA)----
# weights <- seq(from = 0, to=1, by=0.1)
# a3_weights <- seq(from = 0, to=1, by=0.1)
# a3 <- 0
# df_index <- 1
# integrated_gene_lists <- vector(mode = "list", length = length(weights)*length(a3_weights))
# 
# 
# 
# for (x in a3_weights){
#   for (y in weights) {
#     current_ranking <- three_metric_geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3= y)
#     current_ranking <- as.data.frame(current_ranking)
#     integrated_gene_lists[[df_index]] <- current_ranking
#     df_index <- df_index + 1
#   }
# }
# 
# save(integrated_gene_lists, file = "Data/Exported_data/R_objects/integrated_gene_lists_for_all_three_metricsall_intersections_all_three_metrics_1800_gene_subset.RData")

#Subsetting the integrated lists to just the top N number of genes for each list----
gene_lists_to_test <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(integrated_gene_lists))){
  current_gene_list <- integrated_gene_lists[[x]]
  current_gene_list$GeneName <- rownames(current_gene_list)
  gene_names_to_test <- head(current_gene_list, n = 1000)
  gene_lists_to_test[[x]] <- gene_names_to_test
}

save(gene_lists_to_test, file = "Data/Exported_data/R_objects/gene_lists_to_test_all_intersections_MAD_SDES.RData")
#save(gene_lists_to_test, file = "Data/Exported_data/R_objects/gene_lists_to_test_mad_sdes.RData")

#Shrinking the scRNA-seq dataframe down to just the common genes in each of my integrated
#lists and it
genes_of_interest <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(integrated_gene_lists))){
  current_list <- gene_lists_to_test[[x]]
  common_genes <- intersect(colnames(all_tumor_cells_fpkm_denoised_df_cleaned), current_list$GeneName)
  current_genes <- subset(all_tumor_cells_fpkm_denoised_df_cleaned, select = common_genes) 
  genes_of_interest[[x]] <- current_genes
}

save(genes_of_interest, file = "Data/Exported_data/R_objects/genes_of_interest_MAD_SDES.RData")
#save(genes_of_interest, file = "Data/Exported_data/R_objects/genes_of_interest_mad_sdes.RData")

#Now getting the intersection of my gene signature list and the bulk dataframe
#from TCGA's list
all_intersections <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(genes_of_interest))){
  current_set <- genes_of_interest[[x]]
  current_intersect <- intersect(colnames(gene_expression_info), colnames(current_set))
  all_intersections[[x]] <- current_intersect
}

save(all_intersections, file = "Data/Exported_data/R_objects/all_intersections_MAD_SDES.RData")
#save(all_intersections, file = "Data/Exported_data/R_objects/all_intersections_mad_sdes.RData")

#Ensuring that any 'bad' characters in my gene list are removed
all_intersections_cleaned <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(all_intersections))){
  genes_in_bulk_RNA <- all_intersections[[x]]
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="-",replacement=".")
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="_", replacement=".")
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="/", replacement=".")
  all_intersections_cleaned[[x]] <- genes_in_bulk_RNA
}

save(all_intersections_cleaned, file = "Data/Exported_data/R_objects/all_intersections_cleaned_MAD_SDES.RData")
#save(all_intersections_cleaned, file = "Data/Exported_data/R_objects/all_intersections_cleaned_mad_sdes.RData")




#Cox models----
set.seed(1)
cox_models <- vector(mode = "list", length = length(all_intersections_cleaned))
f_objects <- vector(mode = "list", length = length(all_intersections_cleaned))
lambdas <- vector(mode = "list", length = length(all_intersections_cleaned))
c_indicies <- vector(mode = "list", length = length(all_intersections_cleaned))
all_coefs <- vector(mode = "list", length = length(all_intersections_cleaned))
all_active_coefs <- vector(mode = "list", length = length(all_intersections_cleaned))
all_formulas <- vector(mode = "list", length = length(all_intersections_cleaned))

#Making the y for the cox model----
my_y <- Surv(time = df_for_train_test_split$days.to.last.follow.up, event = df_for_train_test_split$vital.status)

#The cox model loop----
for (x in seq(1:length(all_intersections_cleaned))){
  current_formula_data <- all_intersections_cleaned[[x]]
  current_formula_data <- as.vector(current_formula_data)
  #For just MAD metric
  # current_formula_data <- rownames(t(mad.ranking.subset.df))
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))

  #For just SDES metric
  # current_formula_data <- rownames(t(vim.sdes.ranking.subset.df))
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))

  #For just MiRNA metric
  current_formula_data <- rownames(t(mirna.ranking.subset.df))
  current_formula_data <- intersect(current_formula_data, colnames(merged_df))

  
  #For just scDD method
  # scdd_subset <- head(res$gene, n=150)
  # current_formula_data <- scdd_subset
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))
  # 
  # #For just scDD method cell lines
  # scdd_subset_cell <- head(scdd_res_cell_line$gene, n=150)
  # current_formula_data <- scdd_subset_cell
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))
  # 
  # #For just DEsingle CC patient data
  # des_subset <- head(rownames(des_results), n=150)
  # current_formula_data <- des_subset
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))
  # 
  # #For just DEsingle cell-line data
  # des_subset_cell <- head(rownames(results.sig.h1.gm), n=150)
  # current_formula_data <- des_subset_cell
  # current_formula_data <- intersect(current_formula_data, colnames(merged_df))
  # 
  
  
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
  all_active_coefs[[x]] <- Active.Coefficients
  all_formulas[[x]] <- my_formula
}

save(all_formulas, file = "Data/Data-from-pipeline/all_formulas_MAD_SDES.RData")
save(cox_models, file = "Data/Data-from-pipeline/cox_models_MAD_SDES.RData")
save(f_objects, file = "Data/Data-from-pipeline/f_objects_MAD_SDES.RData")
save(lambdas, file = "Data/Data-from-pipeline/lambdas_MAD_SDES.RData")
save(c_indicies, file = "Data/Data-from-pipeline/c_indicies_MAD_SDES.RData")
save(all_coefs, file = "Data/Data-from-pipeline/all_coefs_MAD_SDES.RData")
save(all_active_coefs, file = "Data/Data-from-pipeline/all_active_coefs_MAD_SDES.RData")

for(x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  print(current_model)
}
