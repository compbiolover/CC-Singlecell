#Name:gene_size.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Look at Concordance index across input gene size

#Needed functions----
cox_model_fitter <- function(my.seed       =1,
                             cox.df        =NULL,
                             gene.num      =,
                             cox.predictors=NULL,
                             tumor.stage   =FALSE,
                             tumor.n       =FALSE,
                             tumor.m       =FALSE,
                             regular.cox   =TRUE,
                             save.regular.cox.genes =TRUE,
                             remove.stage   = c("tumor.stage1","tumor.stage2","tumor.stage3", "tumor.stage4"),
                             remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2", "ajcc.n3"),
                             my.filename   ="my_saved_genes.csv"){
  
  #Doing input sanity checks----
  if(missing(my.seed)){
    stop("You must specify a seed for reproducible analysis.")
  }
  
  if(missing(cox.df)){
    stop("You must include a data.frame to carry out cox analysis.")
  }
  
  if(missing(cox.predictors)){
    stop("You must include predictors for your cox model.")
  }
  
  if (class(my.seed)!="numeric"){
    stop("You must specify a number for your seed.")
  }
  
  if(class(cox.df)!="data.frame"){
    stop("Your cox.df parameter must be of type 'data.frame'.")
  }
  
  #If the packages are installed, they----
  #will be loaded. If they are not,
  #the packages will be installed
  #from CRAN and then loaded.
  require(BiocGenerics)
  require(doParallel)
  require(glmnet)
  require(lmtest)
  require(parallel)
  require(survival)
  require(survminer)
  
  #Setting the seed for reproducible----
  #output
  set.seed(my.seed)
  
  #Setting the number of processors on the----
  #machine to speed up the fitting
  num_of_cores <- parallel::detectCores()
  registerDoParallel(cores = num_of_cores)
  
  #Making the list that will store the----
  #data we will return
  cox_data <- vector(mode = "list", length = 9)
  
  #The predictors for the cox model----
  if(is.numeric(cox.predictors)==TRUE){
    my_predictors <- names(cox.predictors)
  }else{
    my_predictors <- rownames(cox.predictors)
  }
  
  if(is.character(cox.predictors)==TRUE){
    my_predictors <- cox.predictors
  }
  
  
  my_predictors <- head(my_predictors, n=gene.num)
  # my_predictors <- sapply(my_predictors, gsub, pattern="-",replacement=".")
  # my_predictors <- unlist(my_predictors)
  # colname_changes <- sapply(colnames(cox.df), gsub, pattern="-",replacement=".")
  # colname_changes <- sapply(colnames(cox.df), gsub, pattern="_",replacement=".")
  # colname_changes <- sapply(colnames(cox.df), gsub, pattern="/",replacement=".")
  # colname_changes <- unlist(colname_changes)
  # colnames(cox.df) <- colname_changes
  my_predictors <- intersect(my_predictors, colnames(cox.df))
  my_predictors <- paste("~", paste(my_predictors[1:length(my_predictors)], collapse = "+"))
  
  
  if(tumor.stage==TRUE & tumor.n==FALSE & tumor.m==FALSE){
    my_predictors <- paste(my_predictors, "tumor.stage", sep = "+")
    my_predictors <- as.formula(my_predictors)
  }else if (tumor.stage==TRUE & tumor.n==TRUE & tumor.m==FALSE){
    my_predictors <- paste(my_predictors,"ajcc.n", sep = "+")
    my_predictors <- as.formula(my_predictors)
  }else if (tumor.stage==TRUE & tumor.n==TRUE & tumor.m==TRUE){
    my_predictors <- paste(my_predictors, "ajcc.m", sep = "+")
    my_predictors <- as.formula(my_predictors)
  }else{
    my_predictors <- as.formula(my_predictors)
    #print("This is the genes only predictor....")
  }
  
  my_x <- model.matrix(my_predictors, cox.df)
  #The response object for the cox model----
  my_y <- Surv(time = cox.df$days.to.last.follow.up, event = cox.df$vital.status)
  #The 10-fold cross-validation fit----
  cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = 10, type.measure = "C", maxit=100000, family="cox", parallel = TRUE)
  
  
  #Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family =  "cox", maxit = 100000)
  Coefficients <- coef(fit, s = cv_fit$lambda.min)
  
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  active_genes <-rownames(Coefficients)[Active.Index]
  
  if(tumor.stage==TRUE & tumor.n==FALSE & tumor.m==FALSE){
    print("This is the genes + tumor stage predictor...")
    active_genes <- active_genes[!active_genes %in% remove.stage]
    active_genes <- c(active_genes, "tumor.stage")
  }
  
  if(tumor.stage==TRUE & tumor.n==TRUE & tumor.m==FALSE){
    print("This is the genes + tumor stage + n stage predictor...")
    active_genes <- active_genes[!active_genes %in% remove.n.stage]
    active_genes <- c(active_genes, "tumor.stage", "ajcc.n")
  }
  
  #Getting survival
  my_surv <- survival::survfit(cv_fit, s = "lambda.min", x = my_x, y = my_y)
  
  #Regular Cox
  if(regular.cox==TRUE){
    active_predictors <-paste(active_genes, collapse = "+")
    regular_cox_df <- cox.df[,active_genes]
    regular_cox_df$days.to.last.follow.up <- cox.df$days.to.last.follow.up
    regular_cox_df$vital.status <- cox.df$vital.status
    my_formula <- paste("~", paste(active_genes[1:length(active_genes)], collapse = "+"))
    regular_cox <- coxph(Surv(time = regular_cox_df$days.to.last.follow.up, event = regular_cox_df$vital.status)~., data = regular_cox_df)
    cox_data[["Regular Cox"]] <- regular_cox
    cox_data[["Summary of Regular Cox"]] <- summary(regular_cox)
    if(save.regular.cox.genes==TRUE){
      summ_regular_cox <- summary(regular_cox)
      write.csv(summ_regular_cox$coefficients, file = my.filename)
    }
    
  }
  
  #Adding the relevant data bits to list to return
  cox_data[["CV"]] <- cv_fit
  cox_data[["Coefficients"]] <- Coefficients
  cox_data[["Active Coefficients"]] <- Active.Coefficients
  cox_data[["Active Index"]] <- Active.Index
  cox_data[["Active Genes"]] <- active_genes
  cox_data[["Predictors"]] <- my_predictors
  cox_data[["Predicted Survival"]] <- my_surv
  
  
  #Returning our finished output----
  return(cox_data)
  
}


#Function used to rank my linear model with just two metrics----
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

#Function used to rank my linear model with three metrics----
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
  integrated_gene_lists <- list()
  
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
three_weight_optimizer <- function(my.start        =0,
                                   my.finish       =1,
                                   step.size       =0.1,
                                   a3.start        =0,
                                   my.index        =1,
                                   my.list.length  =121,
                                   first.metric    =mad.ranking,
                                   second.metric   =vim.sdes.ranking,
                                   third.metric    =mirna.ranking,
                                   my.filename     ="Data/Exported-data/R-objects/three_metric_optimization.RData"){
  
  #Setting up weights, the starting value for the a3 metric, loop-indexer, and the list that will store the results
  weights <- seq(from = my.start, to=my.finish, by=step.size)
  a3_weights <- seq(from = my.start, to=my.finish, by=step.size)
  a3 <- a3.start
  df_index <- my.index
  integrated_gene_lists <- list()
  
  #Doing the grid search
  for (x in a3_weights){
    for (y in weights) {
      current_ranking <- three_metric_geneRank(ranking1 = first.metric, ranking2 = second.metric, ranking3 = third.metric,  a1=x, a2=1-(x+a3), a3= y)
      current_ranking <- as.data.frame(current_ranking)
      integrated_gene_lists[[as.character(df_index)]] <- current_ranking
      df_index <- df_index + 1
    }
  }
  
  #Saving the output of the grid search to a .RData
  save(integrated_gene_lists, file = my.filename)
  return(integrated_gene_lists)
}


#Needed files----
load("gbm_df.RData", verbose = TRUE)
load("mad.RData", verbose = TRUE)
load("sde.RData", verbose = TRUE)
#load("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Inputs/RData/800_510_targets.RData", verbose = TRUE)
#load("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Inputs/RData/800_510_targets.RData", verbose = TRUE)
#load("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Weighted-Optimizations/Optimization_500_810_targets.RData", verbose = TRUE)
#mad_mirna_sdes_optimized <- integrated_gene_lists
#Weighting----
files <- list.files("/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Inputs/RData/")
for(f in files[18:33]){
  current_file <- paste0("/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Inputs/RData/", f)
  mirna.ranking <- load(file = current_file)
  mad_mirna_sdes_optimized <- three_weight_optimizer(first.metric = mad.genes,
                                                     second.metric = sde.genes,
                                                     third.metric = mirna.ranking,
                                                     my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Weighted-Optimizations/Optimization_",f))

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
                                    remove.stage = c("tumor.stage1","tumor.stage2","tumor.stage3", "tumor.stage4"),
                                    remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2", "ajcc.n3"),
                                    my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Outputs/Gene-size/GBM_",gene_size,"_model_weight_active_gene_",iw,".csv"))


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

  write.csv(cindices_df, "/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Outputs/Gene-size/cindices_df.csv")
  write.csv(active_coefs, "/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Outputs/Gene-size/active_coefs_df.csv")



}
# mad_mirna_sdes_optimized <- three_weight_optimizer(first.metric = mad.genes,
#                                            second.metric = sde.genes,
#                                            third.metric = mirna.ranking,
#                                            my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Weighted-Optimizations/Optimization_1000_10_targets.RData"))
# 
# 

#Doing 2 weight optimization----
# mad_mirna_sdes_optimized <- two_weight_optimizer(first.metric = mad.genes,
#                                                    second.metric = mirna.ranking,
#                                                    my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Weighted-Optimizations/Optimization_CC_singlecell_MM.RData"))
# 


#Gene size vs. Concordance index----
# gene_size <- seq(100, 3000, by = 50)
# cindices <- c()
# gene_sizes <- c()
# active_coefs <- c()
# 
# for (gs in gene_size){
#   for(iw in mad_mirna_sdes_optimized){
#     cox_model <- cox_model_fitter(my.seed = 1,
#                                   cox.df = cox_df,
#                                   gene.num = gs,
#                                   cox.predictors = iw,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   regular.cox = FALSE,
#                                   save.regular.cox.genes = FALSE,
#                                   remove.stage = c("tumor.stage1","tumor.stage2","tumor.stage3", "tumor.stage4"),
#                                   remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2", "ajcc.n3"),
#                                   my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Outputs/Gene-size/READ_",gene_size,"_model_weight_active_gene_",iw,".csv"))
# 
# 
#     c_finder <-cox_model$CV$index[1]
#     current_c <- cox_model$CV$cvm[c_finder]
#     current_c <- round(current_c, digits = 4)
#     cindices <- c(cindices, current_c)
#     gene_sizes <- c(gene_sizes, gs)
#     current_coefs <- length(cox_model$`Active Coefficients`)
#     active_coefs <- c(active_coefs, current_coefs)
# 
#   }
# 
# }
# 

#Just gene size for individual metrics-----
#load("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Inputs/RData/500_810_targets.RData", verbose = TRUE)
# for (gs in gene_size){
#   cox_model <- cox_model_fitter(my.seed = 1,
#                                 cox.df = cox_df,
#                                 gene.num = gs,
#                                 cox.predictors = mirna.ranking,
#                                 tumor.stage = FALSE,
#                                 tumor.n = FALSE,
#                                 tumor.m = FALSE,
#                                 regular.cox = FALSE,
#                                 save.regular.cox.genes = FALSE,
#                                 remove.stage = c("tumor.stage1","tumor.stage2","tumor.stage3", "tumor.stage4"),
#                                 remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2", "ajcc.n3"),
#                                 my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Gene-size/mirna_alone",gene_size,"_model_weight_active_gene_.csv"))
#   
#   
#   
#   c_finder <-cox_model$CV$index[1]
#   current_c <- cox_model$CV$cvm[c_finder]
#   current_c <- round(current_c, digits = 4)
#   cindices <- c(cindices, current_c)
#   gene_sizes <- c(gene_sizes, gs)
#   
# }




cindices_df <- cbind(gene_sizes, cindices)

write.csv(cindices_df, "/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Gene-size/cindices_df_mirna_mad_only.csv")
write.csv(active_coefs, "/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Gene-size/active_coefs_df_mirna_mad.csv")

df_to_plot <- data.frame(data=aggregate(x = read_mad_only$cindices,
                                        by = list(read_mad_only$gene_sizes),
                                        FUN = mean))

write.csv(df_to_plot, "Data/Other-methods/desingle/Outputs/cindices_df_to_plot_desingle.csv")


