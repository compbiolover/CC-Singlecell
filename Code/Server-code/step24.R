#Loading the needed functions from their respective files
source("cox_model.R")
source("model_optimizer.R")

load("gbm_df.RData", verbose = TRUE)



#Getting ideal gene number for miRNA metric on TCGA-GBM
#For this metric we don't know which combination of miRNA and miRNA targets
#will yield the best result in advance so we are trying 3 different miRNA-
#miRNA target combinations that encompass specific areas of our grid search.
#They are low miRNA number and low miRNA target number, medium miRNA number and
#medium miRNA target number, and high miRNA number and high miRNA target number.
#Once the ideal miRNA-miRNA target pair is known from the grid search we will
#also include a gene size search for it here. The top performer is 0.7512
gene_sizes <- seq(100, 3000, 50)
mirna_high_cindices <- rep(0, length(gene_sizes))


#Loading the high miRNA-miRNA target file
# load(file = "300_1010_targets.RData",
#      verbose = TRUE)


# load(file = "200_1010_targets.RData",
#      verbose = TRUE)


# load(file = "100_110_targets.RData",
#      verbose = TRUE)

load(file = "100_1010_targets.RData",
     verbose = TRUE)

# load(file = "100_910_targets.RData",
#      verbose = TRUE)

high.mirna.genes <- mirna.ranking

for(gs in gene_sizes){
  cox_model <- cox_model_fitter(my.seed = 1,
                                cox.predictors = high.mirna.genes,
                                cox.df = cox_df,
                                gene.num = gs,
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                my.filename = paste0("high_mirna_gbm_coefs_",gs,"_genes.csv"))
  
  #Getting the top concordance index from the cross validation and then rounding
  #it to 4 digits to follow cv.glmnet reporting convention. Finally, we update
  #the mirna_high_cindices list with the result
  top_cindex <- round(cox_model$CV$cvm[cox_model$CV$index[1]], digits = 4)
  mirna_high_cindices[which(gene_sizes==gs)] <- top_cindex
  
}

#Binding the gene size and c-index vectors together to get the finished data
#frame
mirna_high_cindices_read_df <- as.data.frame(cbind(gene_sizes,
                                                   mirna_high_cindices))

colnames(mirna_high_cindices_read_df) [2] <- "c_index"
mirna_high_cindices_read_df$mirna_type <- rep("high",
                                              nrow(mirna_high_cindices_read_df))

write.csv(mirna_high_cindices_read_df, file = "mirna_high_cindices_gbm_across_gene_size_100_910_mirna.csv")
