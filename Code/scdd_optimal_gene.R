#Name: scdd_optimal_gene.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: To find the ideal number of genes for the scDD method for datasets


#Loading needed packages----
library(scDD)
library(tidyverse)

#Loading needed data----
#load("coad_df.RData", verbose = TRUE)

#Looping over each unique input----
scdd_gene_nums <- seq(1, 570, 50)
cindices <- list()
gene_sizes <- list()
for(gene_num in scdd_gene_nums){
  current_num <- read.csv(paste0("Data/Data-from-Cleaner-code/scdd_top_",gene_num,"_genes_cc_patients.csv"))
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.df = cox_df,
                                gene.num = gene_num,
                                cox.predictors = current_num$gene,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                regular.cox = FALSE,
                                save.regular.cox.genes = FALSE,
                                remove.stage = c("tumor.stage1","tumor.stage2","tumor.stage3", "tumor.stage4"),
                                remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2", "ajcc.n3"),
                                my.filename = paste0("Data/TCGA-COAD/Data/Other-methods/scdd/Outputs/scdd_",gene_num, "genes.csv"))
  
  
  
  c_finder <-cox_model$CV$index[1]
  current_c <- cox_model$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  cindices <- c(cindices, current_c)
  gene_sizes <- c(gene_sizes, gene_num)
  
}



#Binding all of the c-indices results together with all of the gene sizes
cindices_df <- cbind(gene_sizes, cindices)

#Writing out the finished c-indicies data frame
write.csv(cindices_df, file = "Data/Other-methods/scdd/Outputs/cindices_df_scdd_rc.csv")
