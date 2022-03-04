#Name: gene_size_new.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: To find the optimal number of genes for each of our methods

#Loading needed functions----
source("Code/cox_model.R")
source("Code/model_optimizer.R")


#Loading needed files----
cox_df <- readRDS("read_df_finished_v2.rds")

#Looping through all of the optimized alpha model files----
files <- list.files("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Weighted-Optimizations")

for(f in files){
  current_file <- paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Weighted-Optimizations/", f)
  load(file = current_file, verbose = TRUE)
  mad_mirna_sdes_optimized <- integrated_gene_lists
  
  gene_size <- seq(100, 3000, by = 50)
  cindices <- c()
  gene_sizes <- c()
  active_coefs <- c()
  
  for (gs in gene_size){
    for(iw in mad_mirna_sdes_optimized){
      cox_model <- cox_model_fitter(my.seed = 1,
                                    cox.df = cox_df,
                                    gene.num = gs,
                                    cox.predictors = iw,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    regular.cox = FALSE,
                                    save.regular.cox.genes = FALSE,
                                    remove.stage = c("tumor.stage1",
                                                     "tumor.stage2",
                                                     "tumor.stage3",
                                                     "tumor.stage4"),
                                    remove.n.stage = c("ajcc.n0", "ajcc.n1",
                                                       "ajcc.n2", "ajcc.n3"),
                                    my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Gene-size/READ_",gene_size,"_model_weight_active_gene_",iw,".csv"))
      
      
      c_finder <-cox_model$CV$index[1]
      current_c <- cox_model$CV$cvm[c_finder]
      current_c <- round(current_c, digits = 4)
      cindices <- c(cindices, current_c)
      gene_sizes <- c(gene_sizes, gs)
      current_coefs <- length(cox_model$`Active Coefficients`)
      active_coefs <- c(active_coefs, current_coefs)
      
    }
    
  }
  
  cindices_df <- cbind(gene_sizes, cindices)
  
  write.csv(cindices_df, "/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Gene-size/cindices_df.csv")
  write.csv(active_coefs, "/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Gene-size/active_coefs_df.csv")
  
  
  
  
  
}