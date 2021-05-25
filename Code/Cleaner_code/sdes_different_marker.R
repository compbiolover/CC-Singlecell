#sdes_different_marker.R

#Load the different marker files----
load("cc_tumor_fpkm_sde_snai1.RData", verbose = TRUE)
snai1_sde <- sde.genes
load("cc_tumor_fpkm_sde_twist1.RData", verbose = TRUE)
twist1_sde <- sde.genes
load("cc_tumor_fpkm_sde_twist1-redo.RData", verbose = TRUE)
twist1_redo_sde <- sde.genes
load("cc_tumor_fpkm_sde_zeb1.RData", verbose = TRUE)
zeb1_sde <- sde.genes
load("cc_tumor_fpkm_sde.RData", verbose = TRUE)
vim_sde <- sde.genes
load("cox_df.RData", verbose = TRUE)


#Loading needed packages----
library(survival);packageVersion("survival")
library(survminer);packageVersion("survminer")

#Other required packages----
require(BiocGenerics)
require(doParallel)
require(glmnet)
require(parallel)
require(survival)

#Setting the number of processors on the----
#machine to speed up the fitting
num_of_cores <- parallel::detectCores()
registerDoParallel(cores = num_of_cores)


my_gene_list <- list(snai1=snai1_sde, twist1=twist1_sde, twist1_redo=twist1_redo_sde, zeb1=zeb1_sde, vim=vim_sde)


#Getting access to cox function from other script----
if(!exists("cox_model_fitter", mode="function")) source("cox_model.R")

#Setting the gene size to a constant value----
cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1

for (y in my_gene_list) {
  print(counter)
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = y , tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  
}

print("Saving Different SDES Marker Results")
save(cox_models, file = "different_sdes_markers_rc_patients.RData")

my_df <- data.frame(concordance_index=my_cindicies)
write.csv(my_df, file = "different_sdes_markers_rc_patients.csv")




