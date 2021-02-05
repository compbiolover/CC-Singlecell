#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Contains just code needed to optimize linear model on
#lab server

#Setting the working directory to where our project files are----
setwd("/home/awillems/Projects/CC_Singlecell")
#Loading the needed packages----
library(glmnet)
library(survival)
library(tidyverse)
library(utils)

#Loading needed functions----
# Ranking genes by three measurements
# geneRank <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, a1 = 1,
#                      a2 = 0, a3 = 0){
#   gns = c(names(ranking1),names(ranking2),names(ranking3))
#   gns = unique(gns)
#   res = rep(0, length(gns))
#   names(res) = gns
#   for(i in names(res)) {
#     res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2+getRank(ranking3, i)*a3
#   }
#   #res=res/sum(res)
#   #res = res[order(res, decreasing = T)]
#   res
# }
# 
# getRank <- function(ranking = NULL, gn = NULL){
#   if (gn %in% names(ranking)) {
#     return(ranking[gn])
#   }
#   else return(0.0)
# }
# 


geneRank <- function(ranking1 = NULL,
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
   #res=res/sum(res) 
   res = res[order(res, decreasing = T)]
   res
 }
 
getRank <- function(ranking = NULL,
                    gn      = NULL){
   if (gn %in% names(ranking)) {
     return(ranking[gn])
   }
   else return(0.0)
 }


#Loading the needed files for the optimization----
load(file = "Data/all_intersections_cleaned.RData")
load(file = "Data/all_intersections.RData")
load(file = "Data/all_tumor_cells_fpkm_denoised.RData")
load(file = "Data/bulk_rna_df.RData")
load(file = "Data/COA_data_df.RData")
load(file = "Data/COA_data_se.RData")
load(file = "Data/COA_survival_data.RData")
load(file = "Data/df_for_train_test_split.RData")
load(file = "Data/finished_sets.RData")
load(file = "Data/gene_expression_info.RData")
load(file = "Data/gene_lists_to_test.RData")
load(file = "Data/genes_of_interest.RData")
load(file = "Data/mad.ranking.RData")
load(file = "Data/mirna.ranking.RData")
load(file = "Data/vim.sdes.ranking.RData")

#Code for the grid search linear model optimization----
# weights <- seq(from = 0, to=1, by=0.1)
# df_index <- 1
# integrated_gene_lists <- list()
# a3_weights <- seq(from = 0, to=1, by=0.1)
# a3 <- 0
# 
# for (x in a3_weights){
#   print(x)
#   for (y in weights) {
#     print(y)
#     current_ranking <- geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3=y)
#     current_ranking <- as.data.frame(current_ranking)
#     integrated_gene_lists[[df_index]] <- current_ranking
#     df_index <- df_index + 1
#   }
# }

#Code for just two metrics----
weights <- seq(from = 0, to=1, by=0.1)
df_index <- 1
integrated_gene_lists <- list()

for (x in weights) {
  print(x)
  current_ranking <- geneRank(ranking1 = vim.sdes.ranking, ranking2 = mirna.ranking,  a1=x, a2=1-x)
  current_ranking <- as.data.frame(current_ranking)
  integrated_gene_lists[[df_index]] <- current_ranking
  df_index <- df_index + 1
}


# 
# gene_lists_to_test <- list()
# all_tumor_cells_fpkm_denoised_df <- as.data.frame(all_tumor_cells_fpkm_denoised)
# all_tumor_cells_fpkm_denoised_df <- t(all_tumor_cells_fpkm_denoised_df)
# for (x in seq(1:length(integrated_gene_lists))){
#   print(x)
#   current_gene_list <- as.data.frame(integrated_gene_lists[[x]])
#   current_gene_list$GeneName <- rownames(current_gene_list)
#   gene_names_to_test <- head(current_gene_list, n = 900)
#   gene_lists_to_test[[x]] <- gene_names_to_test
# }

# 
# #Changing the colnames of the signle cell dataframe to the simple gene name so
# #that subsetting works
# current_colname_split <- strsplit(rownames(all_tumor_cells_fpkm_denoised_df), "_")
# finished_gene_list <- c()
# current_list <- current_colname_split
# for (y in seq(1:length(current_list))){
#   #print(current_list[[y]][2])
#   finished_gene_list <- c(finished_gene_list, current_list[[y]][2])
# }
# 
# all_tumor_cells_fpkm_denoised_df <- t(all_tumor_cells_fpkm_denoised_df)
# colnames(all_tumor_cells_fpkm_denoised_df) <- finished_gene_list
# 
# genes_of_interest <- list()
# for (x in seq(1:length(integrated_gene_lists))){
#   current_list <- gene_lists_to_test[[x]]
#   current_list <- as.data.frame(current_list)
#   common_genes <- intersect(colnames(all_tumor_cells_fpkm_denoised_df), current_list$GeneName)
#   current_genes <- subset(all_tumor_cells_fpkm_denoised_df, select = common_genes) 
#   genes_of_interest[[x]] <- current_genes
# }
# 
# all_split_colnames <- list()
# for (x in seq(1:length(integrated_gene_lists))){
#   for (y in seq(1:length(colnames(genes_of_interest[[x]])))){
#     current_df <- genes_of_interest[[x]]
#     current_df <- as.data.frame(current_df)
#     current_colname_split <- strsplit(colnames(current_df), "_")
#     all_split_colnames[[x]] <- current_colname_split
#   }
#   
# }
# 
# finished_sets <- list()
# for (x in seq(1:length(integrated_gene_lists))) {
#   finished_gene_list <- c()
#   current_list <- all_split_colnames[[x]]
#   for (y in seq(1:length(current_list))){
#     finished_gene_list <- c(finished_gene_list, current_list[[y]])
#   }
#   finished_sets[[x]] <-finished_gene_list
# }
# 
# #Making the dataframe that contains just the genes by samples that we need from bulk RNA-seq data
# gene_expression_info <- COA_data_se@assays@data@listData[["HTSeq - Counts"]]
# rownames(gene_expression_info) <- COA_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
# colnames(gene_expression_info) <- COA_data_se@colData@rownames
# gene_expression_info <- t(gene_expression_info)
# 
# 
# all_intersections <- list()
# for (x in seq(1:length(finished_sets))){
#   current_set <- finished_sets[[x]]
#   current_intersect <- intersect(colnames(gene_expression_info), current_set)
#   all_intersections[[x]] <- current_intersect
# }
# 


#Merging the two dataframes together into a larger dataframe that we can use for the Cox PH
bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
COA_data_df_unique <- subset(COA_data_df, select = unique(colnames(COA_data_df)))
merged_df <- merge(bulk_rna_df_unique, COA_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

#Making the dataframe with the survival info
survival_df <- subset(merged_df, select=vital_status)
time <- as.numeric(merged_df$days_to_last_follow_up)
survival_df$time <- cbind(time)
survival_df <- na.omit(survival_df)

# all_intersections_cleaned <- list()
# for (x in seq(1:length(all_intersections))){
#   genes_in_bulk_RNA <- all_intersections[[x]]
#   genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="-",replacement=".")
#   genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="_", replacement=".")
#   genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="/", replacement=".")
#   all_intersections_cleaned[[x]] <- genes_in_bulk_RNA
# }

#save(all_intersections_cleaned, file = "Data/Exported-data/R-objects/all_intersections_cleaned.RData")
rows_to_remove <- setdiff(rownames(merged_df), rownames(survival_df))
merged_df <- merged_df[!(row.names(merged_df) %in% rows_to_remove), ]


colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern="-", replacement=".")
colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern="_", replacement=".")
colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern="/", replacement=".")

#Cox models----
cox_models <- list()
f_objects <- list()
lambdas <- list()
c_indicies <- list()
prediction_res <- list()
all_coefs <- list()
all_active_coefs <- list()
all_formulas <- list()

for (x in seq(1:length(all_intersections_cleaned))){
  current_formula_data <- all_intersections_cleaned[[x]]
  current_formula_data <- as.vector(current_formula_data)
  current_formula_data <- current_formula_data[-107]
  
  my_formula <- paste("~", paste(current_formula_data[1:length(current_formula_data)], collapse = "+"))
  f <- as.formula(my_formula)
  
  
  
  my_x <- model.matrix(f, df_for_train_test_split)
  my_y <- Surv(time = df_for_train_test_split$days.to.last.follow.up, event = df_for_train_test_split$vital.status)
  cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = 10, type.measure = "C", maxit=100000, family="cox")
  current_lambdas <- cv_fit$lambda
  current_c_indicies <- cv_fit$cvm
  
  
  
  cox_models[[x]] <- cv_fit
  f_objects[[x]] <- f
  lambdas[[x]] <- current_lambdas
  c_indicies[[x]] <- current_c_indicies
  
  
  #Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family =  "cox", maxit = 100000)
  Coefficients <- coef(fit, s = cv_fit$lambda.min)
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  all_coefs[[x]] <- Coefficients
  all_active_coefs[[x]] <- Active.Coefficients
  
}
save(all_formulas, file = "Data/Data-from-pipeline/all_formulas.RData")
save(cox_models, file = "Data/Data-from-pipeline/cox_models.RData")
save(f_objects, file = "Data/Data-from-pipeline/f_objects.RData")
save(lambdas, file = "Data/Data-from-pipeline/lambdas.RData")
save(c_indicies, file = "Data/Data-from-pipeline/c_indicies.RData")

for(x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  print(current_model)
}

all_c_indexes_finished <- list()
all_lambdas <-list()
min_lambda <- list()
for (x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  all_c_indexes_finished[[x]] <- current_model$cvm
  all_lambdas[[x]] <- current_model$lambda
  min_lambda[[x]] <- current_model$lambda.min
}
