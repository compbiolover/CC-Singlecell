#Author: Andrew Willems <awillems@vols.utk.edu>.
#Purpose: To take scRNA-seq colon cancer data and begin analysis.
#R version: 4.0.2.

#Loading needed packages----
library(glmnet);packageVersion("glmnet")
library(hoardeR);packageVersion("hoardeR")
library(monocle3);packageVersion("monocle3")
library(Rmagic);packageVersion("Rmagic")
library(SummarizedExperiment);packageVersion("SummarizedExperiment")
library(survival);packageVersion("survival")
library(survminer);packageVersion("survminer")
library(switchde);packageVersion("switchde")
library(TCGAbiolinks);packageVersion("TCGAbiolinks")
library(tidyverse);packageVersion("tidyverse")

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

#Loading the single-cell data files. ----
all_tumor_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_tumor_all_cells_FPKM.csv")


#Loading the .RData files to simplify analysis----
load(file = "Data/Exported-data/R-objects/COA_data_se.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/COA_data_df.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/bulk_rna_df.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/all_tumor_cells_fpkm_denoised_cleaned.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/mad.ranking.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/vim.sdes.ranking.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/mirna.ranking.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/integrated-gene-lists-for-two-metrics.RData", verbose = TRUE) #This is the SDES + MiRNA metrics
load(file = "Data/Exported-data/R-objects/all_tumor_cells_fpkm_denoised_df_cleaned.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/gene_lists_to_test_MiRNA_SDES.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/all_intersections_cleaned-MiRNA-SDES.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/df_for_train_test_split.RData", verbose = TRUE)
load(file = "Data/Exported-data/R-objects/mad.ranking.subset.df.RData", verbose = TRUE)


#Cox models----
set.seed(1)
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
  #current_formula_data <- intersect(colnames(mad.ranking.subset.df), colnames(df_for_train_test_split))
  current_formula_data <- as.vector(current_formula_data)
  #current_formula_data <- current_formula_data[-107]
  
  my_formula <- paste("~", paste(current_formula_data[1:length(current_formula_data)], collapse = "+"))
  all_formulas[[x]] <- my_formula
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


#KM Curves----
#Saving just the active gene names
active_genes <-rownames(Coefficients)[Active.Index]
surv_gene_df=merged_df[,active_genes]
#metric.coef <- as.list(mad.ranking.subset)
metric.coef <- active_genes
integrated_gene_list_info <- integrated_gene_lists[[11]]
integrated_gene_list_info <- t(integrated_gene_list_info)
common_genes <- intersect(metric.coef, colnames(integrated_gene_list_info))
integrated_gene_list_info <- subset(integrated_gene_list_info, select=common_genes)
#metric.coef <- as.data.frame(metric.coef)
metric.coef <- t(metric.coef)
colnames(metric.coef)[1] <- "Score"
gene_weight <- as.data.frame(integrated_gene_list_info)
#gene_weight <- as.data.frame(metric.coef)
gene_weight <- t(gene_weight)
#gene_weight <- subset(gene_weight, select=colnames(surv_gene_df))
#gene_weight <- t(gene_weight)
gene_weight <- sort(gene_weight, decreasing = TRUE)
weight_function=function(x){crossprod(as.numeric(x),gene_weight)}
surv_gene_df <- surv_gene_df[,common_genes]
trainScore=apply(surv_gene_df,1,weight_function)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
surv_gene_df <- cbind(risk, surv_gene_df)
vital.status <- merged_df$vital.status
days.to.last.follow.up <- merged_df$days.to.last.follow.up
surv_gene_df <- cbind(vital.status, surv_gene_df)
surv_gene_df <- cbind(days.to.last.follow.up, surv_gene_df)
km_fit <- survfit(Surv(days.to.last.follow.up, vital.status) ~ risk, data = surv_gene_df)

#P-value calculation for KM curves
diff=survdiff(Surv(days.to.last.follow.up, vital.status) ~risk,data = surv_gene_df)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

#KM Curves plotting code
surPlot<-ggsurvplot(km_fit,
                    data=surv_gene_df,
                    pval=paste0("p=",pValue),
                    pval.size=4,
                    legend.labs=c("High risk", "Low risk"),
                    legend.title="Risk",
                    xlab="Time(months)",
                    palette=c("red", "blue"))






