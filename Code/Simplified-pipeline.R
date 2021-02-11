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
#load(file = "Data/Exported_data/R_objects/vim.sdes.ranking.RData")
load(file = "Data/Exported_data/R_objects/mirna.ranking.RData")
#load(file = "Data/Exported_data/R_objects/integrated-gene-lists-for-two-metrics.RData") #This is the SDES + MiRNA metrics
load(file = "Data/Exported_data/R_objects/all_tumor_cells_fpkm_denoised_df_cleaned.RData")
#load(file = "Data/Exported_data/R_objects/gene_lists_to_test_MiRNA_SDES.RData")
#load(file = "Data/Exported_data/R_objects/all_intersections_cleaned-MiRNA-SDES.RData")
load(file = "Data/Exported_data/R_objects/df_for_train_test_split.RData")
load(file = "Data/Exported_data/R_objects/merged_df_replaced.RData")
#load(file = "Data/Exported_data/R_objects/mad.ranking.subset.df.RData")
#load(file = "Data/Lab-server-data/All-three-metrics/sorted_gene_lists_all_three_metrics.RData")



#Optimization for just 2 metrics----
weights <- seq(from = 0, to=1, by=0.1)
df_index <- 1
integrated_gene_lists <- vector(mode = "list", length = length(weights))

for (x in weights) {
  current_ranking <- two_metric_geneRank(ranking1 = mad.ranking, ranking2 = mirna.ranking,  a1=x, a2=1-x)
  current_ranking <- as.data.frame(current_ranking)
  integrated_gene_lists[[df_index]] <- current_ranking
  df_index <- df_index + 1
}

save(integrated_gene_lists, file = "Data/Exported_data/R_objects/integrated_gene_lists_for_mad_mirna.RData")

#Optimization for all 3 metrics (MAD and SDE and miRNA)----
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
# save(integrated_gene_lists, file = "Data/Exported_data/R_objects/integrated_gene_lists_for_all_three_metrics.RData")

#Subsetting the integrated lists to just the top N number of genes for each list----
gene_lists_to_test <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(integrated_gene_lists))){
  current_gene_list <- integrated_gene_lists[[x]]
  current_gene_list$GeneName <- rownames(current_gene_list)
  gene_names_to_test <- head(current_gene_list, n = 900)
  gene_lists_to_test[[x]] <- gene_names_to_test
}

#save(gene_lists_to_test, file = "Data/Exported_data/R_objects/gene_lists_to_test_all_three_metrics.RData")
save(gene_lists_to_test, file = "Data/Exported_data/R_objects/gene_lists_to_test_mad_mirna.RData")

#Shrinking the scRNA-seq dataframe down to just the common genes in each of my integrated
#lists and it
genes_of_interest <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(integrated_gene_lists))){
  current_list <- gene_lists_to_test[[x]]
  common_genes <- intersect(colnames(all_tumor_cells_fpkm_denoised_df_cleaned), current_list$GeneName)
  current_genes <- subset(all_tumor_cells_fpkm_denoised_df_cleaned, select = common_genes) 
  genes_of_interest[[x]] <- current_genes
}

#save(genes_of_interest, file = "Data/Exported_data/R_objects/genes_of_interest_all_three_metrics.RData")
save(genes_of_interest, file = "Data/Exported_data/R_objects/genes_of_interest_mad_mirna.RData")

#Now getting the intersection of my gene signature list and the bulk dataframe
#from TCGA's list
all_intersections <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(genes_of_interest))){
  current_set <- genes_of_interest[[x]]
  current_intersect <- intersect(colnames(gene_expression_info), colnames(current_set))
  all_intersections[[x]] <- current_intersect
}

#save(all_intersections, file = "Data/Exported_data/R_objects/all_intersections_all_three_metrics.RData")
save(all_intersections, file = "Data/Exported_data/R_objects/all_intersections_mad_mirna.RData")

#Ensuring that any 'bad' characters in my gene list are removed
all_intersections_cleaned <- vector(mode = "list", length = length(integrated_gene_lists))
for (x in seq(1:length(all_intersections))){
  genes_in_bulk_RNA <- all_intersections[[x]]
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="-",replacement=".")
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="_", replacement=".")
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="/", replacement=".")
  all_intersections_cleaned[[x]] <- genes_in_bulk_RNA
}

#save(all_intersections_cleaned, file = "Data/Exported_data/R_objects/all_intersections_cleaned_all_three_metrics.RData")
save(all_intersections_cleaned, file = "Data/Exported_data/R_objects/all_intersections_cleaned_mad_mirna.RData")


#Cox models----
set.seed(1)
cox_models <- vector(mode = "list", length = length(all_intersections_cleaned))
f_objects <- vector(mode = "list", length = length(all_intersections_cleaned))
lambdas <- vector(mode = "list", length = length(all_intersections_cleaned))
c_indicies <- vector(mode = "list", length = length(all_intersections_cleaned))
#prediction_res <- vector(mode = "list", length = length(all_intersections_cleaned))
all_coefs <- vector(mode = "list", length = length(all_intersections_cleaned))
all_active_coefs <- vector(mode = "list", length = length(all_intersections_cleaned))
all_formulas <- vector(mode = "list", length = length(all_intersections_cleaned))

#Making the y of the cox model----
my_y <- Surv(time = df_for_train_test_split$days.to.last.follow.up, event = df_for_train_test_split$vital.status)

#Mad list subset (still need to integrate)
#current_formula_data <- intersect(colnames(mad.ranking.subset.df), colnames(df_for_train_test_split))

#The cox model loop----
for (x in seq(1:length(all_intersections_cleaned))){
  current_formula_data <- all_intersections_cleaned[[x]]
  current_formula_data <- as.vector(current_formula_data)
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

for(x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  print(current_model)
}

#KM Curves----
# #Saving just the active gene names
# active_genes <-rownames(Coefficients)[Active.Index]
# surv_gene_df=merged_df[,active_genes]
# #metric.coef <- as.list(mad.ranking.subset)
# metric.coef <- active_genes
# integrated_gene_list_info <- integrated_gene_lists[[11]]
# integrated_gene_list_info <- t(integrated_gene_list_info)
# common_genes <- intersect(metric.coef, colnames(integrated_gene_list_info))
# integrated_gene_list_info <- subset(integrated_gene_list_info, select=common_genes)
# #metric.coef <- as.data.frame(metric.coef)
# metric.coef <- t(metric.coef)
# colnames(metric.coef)[1] <- "Score"
# gene_weight <- as.data.frame(integrated_gene_list_info)
# #gene_weight <- as.data.frame(metric.coef)
# gene_weight <- t(gene_weight)
# #gene_weight <- subset(gene_weight, select=colnames(surv_gene_df))
# #gene_weight <- t(gene_weight)
# gene_weight <- sort(gene_weight, decreasing = TRUE)
# weight_function=function(x){crossprod(as.numeric(x),gene_weight)}
# surv_gene_df <- surv_gene_df[,common_genes]
# trainScore=apply(surv_gene_df,1,weight_function)
# risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
# surv_gene_df <- cbind(risk, surv_gene_df)
# vital.status <- merged_df$vital.status
# days.to.last.follow.up <- merged_df$days.to.last.follow.up
# surv_gene_df <- cbind(vital.status, surv_gene_df)
# surv_gene_df <- cbind(days.to.last.follow.up, surv_gene_df)
# km_fit <- survfit(Surv(days.to.last.follow.up, vital.status) ~ risk, data = surv_gene_df)
# 
# #P-value calculation for KM curves
# diff=survdiff(Surv(days.to.last.follow.up, vital.status) ~risk,data = surv_gene_df)
# pValue=1-pchisq(diff$chisq,df=1)
# pValue=signif(pValue,4)
# pValue=format(pValue, scientific = TRUE)
# 
# #KM Curves plotting code
# surPlot<-ggsurvplot(km_fit,
#                     data=surv_gene_df,
#                     pval=paste0("p=",pValue),
#                     pval.size=4,
#                     legend.labs=c("High risk", "Low risk"),
#                     legend.title="Risk",
#                     xlab="Time(months)",
#                     palette=c("red", "blue"))
# 
# 
# 
# 
# 
# 
