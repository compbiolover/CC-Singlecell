#Loading needed packages----
library(ggplot2)
library(lmtest)
library(tidyverse)
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
  counter <- 1
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
                                   my.filename     ="Data/Exported-data/R-objects/three_metric_optimization.RData",
                                   show.perc       =TRUE){
  
  #Setting up weights, the starting value for the a3 metric, loop-indexer, and the list that will store the results
  weights <- seq(from = my.start, to=my.finish, by=step.size)
  a3_weights <- seq(from = my.start, to=my.finish, by=step.size)
  a3 <- a3.start
  df_index <- my.index
  integrated_gene_lists <- list()
  counter <- 1
  total_length <- 121
  
  #Doing the grid search
  for (x in a3_weights){
    for (y in weights) {
      current_ranking <- three_metric_geneRank(ranking1 = first.metric, ranking2 = second.metric, ranking3 = third.metric,  a1=x, a2=1-(x+a3), a3= y)
      current_ranking <- as.data.frame(current_ranking)
      integrated_gene_lists[[as.character(df_index)]] <- current_ranking
      if(show.perc==TRUE){
        total_percent_done <- counter/total_length*100
        total_percent_done <- round(total_percent_done, digits = 2)
        print(paste0(total_percent_done, "% of the total grid search is done."))
        counter <- counter + 1
      }
      df_index <- df_index + 1
    }
  }
  
  #Saving the output of the grid search to a .RData
  save(integrated_gene_lists, file = my.filename)
  return(integrated_gene_lists)
}

#cox_model_fitter----
cox_model_fitter <- function(my.seed       =1,
                             cox.df        =NULL,
                             gene.num      =900,
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
  cox_data <- list()
  
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
  my_predictors <- sapply(my_predictors, gsub, pattern="-",replacement=".")
  my_predictors <- unlist(my_predictors)
  colname_changes <- sapply(colnames(cox.df), gsub, pattern="-",replacement=".")
  colname_changes <- sapply(colnames(cox.df), gsub, pattern="_",replacement=".")
  colname_changes <- sapply(colnames(cox.df), gsub, pattern="/",replacement=".")
  colname_changes <- unlist(colname_changes)
  colnames(cox.df) <- colname_changes
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



#Loading the metric files----
#load("coad_df.RData", verbose = TRUE)
# load(file = "mad.RData", verbose = TRUE)
# load(file = "sde.RData", verbose = TRUE)
#load(file = "mirna_genes_20_100_targets.RData", verbose = TRUE)

#MiRNAs----
#For COAD data set
mirnas <- seq(1, 1560, by = 50)
#For READ dataa set
#mirnas <- seq(1, 1074, by = 50)
load("coad_df.RData", verbose = TRUE)
load(file = "mad.RData", verbose = TRUE)
load(file = "sde.RData", verbose = TRUE)
#load(file = "mirna_genes_40_100_targets_colorect.RData", verbose = TRUE)
#mirna.ranking <- my_gene_df_finished

#For all miRNAs where we keep # of miRNAs constant but vary miRNA targets
# for(y in mirnas[1:22]){
#   for(x in mirnas[1:22]){
#     mirna.ranking <- load(file = paste0("MiRNA/mirna_genes_global_search_mirna_fill_in_", y, "_",x,"_targets.RData"), verbose = TRUE)
#     mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
#                                                        second.metric = mirna.ranking,
#                                                        third.metric = sde.genes,
#                                                        my.filename = paste0("MiRNA/Global-heatmap/Weighted-Optimizations/three_weight_optimized_global_search_fill_in_",y,"_mirnas_",x,"_targets.RData"))
# 
# 
#     #load(file = "/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Weighted-Optimizations/three_weight_optimized_global_search_fill_in_1_mirnas_1_targets.RData", verbose = TRUE)
#     
#     #mad_sdes_mirna_optimized <- integrated_gene_lists
#     
#     #Cox Model
#     cox_models <- list()
#     my_cindicies <- c()
#     counter <- 1
#     
#     for (z in mad_sdes_mirna_optimized[1:121]) {
#       current_weight <- z
#       current_cox <- cox_model_fitter(my.seed = 1,
#                                       cox.df = cox_df,
#                                       gene.num = 1800,
#                                       cox.predictors = current_weight,
#                                       tumor.stage = FALSE,
#                                       tumor.n = FALSE,
#                                       tumor.m = FALSE,
#                                       regular.cox = FALSE,
#                                       save.regular.cox.genes = FALSE,
#                                       my.filename = paste0("three_weight_optimized_global_search_fill_in_",y,"_mirna_",x,"_top.csv"))
#       
#       cox_models[[as.character(counter)]] <- current_cox
#       perc_done <- (counter/121)*100
#       print(perc_done)
#       counter <- counter + 1
#       
#       #Storing all of the c-index values in a vector that we can use later to build the plot
#       c_finder <-current_cox$CV$index[1]
#       current_c <- current_cox$CV$cvm[c_finder]
#       current_c <- round(current_c, digits = 4)
#       my_cindicies <- c(my_cindicies, current_c)
#     }
#     
#     
#     #write.csv(my_cindicies, file = paste0("three_weight_cindices_",x,"_mirna_",x,"_top.csv"))
#     top_cindex <-max(my_cindicies)
#     top_index <- which(my_cindicies==top_cindex)
#     print(top_index)
#     top_index_used <- top_index[1]
#     
#     #Now doing it for just the top c-index combo
#     cox_models <- list()
#     my_cindicies <- c()
#     counter <- 1
#     
#     for (z in mad_sdes_mirna_optimized[top_index_used]) {
#       current_weight <- z
#       current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
#                                       cox.predictors = current_weight,
#                                       tumor.stage = FALSE,
#                                       tumor.n = FALSE,
#                                       tumor.m = FALSE,
#                                       regular.cox = TRUE,
#                                       save.regular.cox.genes = TRUE,
#                                       my.filename = paste0("MiRNA/Global-heatmap/Top-cox-models/cc_singlecell_mms_global_",y,"_mirnas_",x,"_mirna_targets_top_performer_index_",top_index_used,".csv"))
#       
#       cox_models[[as.character(counter)]] <- current_cox
#       
#       
#       #Storing all of the c-index values in a vector that we can use later to build the plot
#       c_finder <-current_cox$CV$index[1]
#       current_c <- current_cox$CV$cvm[c_finder]
#       current_c <- round(current_c, digits = 4)
#       my_cindicies <- c(my_cindicies, current_c)
#       cox_models$`1`$CV
#       c_index_df <- data.frame(c_index=my_cindicies)
#       write.csv(c_index_df, file = paste0("MiRNA/Global-heatmap/Finished-outputs/cindices_global_heatmap_",y,"_mirna_",x,"_mirna_targets_index_",top_index_used,".csv"))
#       
#       counter <- counter + 1
#     }
#     
#   }
#   
#   
#   
# }

#For all miRNAs at a changing size
for(x in mirnas[1:32]){
  mirna.ranking <- load(file = paste0("MiRNA/Global-heatmap/Initial-inputs/mirna_genes_global_search_mirna_", x, "_",x,"_targets.RData"), verbose = TRUE)
  mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
                                                     second.metric = mirna.ranking,
                                                     third.metric = sde.genes,
                                                     my.filename = paste0("MiRNA/Global-heatmap/Weighted-Optimizations/three_weight_optimized_global_search_",x,"_mirnas_",x,"_targets.RData"))

#Cox Model----
  cox_models <- list()
  my_cindicies <- c()
  counter <- 1
  
  load(file =paste0("MiRNA/Global-heatmap/Weighted-Optimizations/three_weight_optimized_global_search_",x,"_mirnas_",x,"_targets.RData"), verbose = TRUE)
  
  for (y in mad_sdes_mirna_optimized[1:121]) {
    current_weight <- y
    current_cox <- cox_model_fitter(my.seed = 1,
                                    cox.df = cox_df,
                                    gene.num = 1800,
                                    cox.predictors = current_weight,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    regular.cox = FALSE,
                                    save.regular.cox.genes = FALSE,
                                    my.filename = paste0("three_weight_optimized_global_search_",x,"_mirna_",x,"_top.RData"))

    cox_models[[as.character(counter)]] <- current_cox
    perc_done <- (counter/121)*100
    print(paste0(round(perc_done, digits = 2), "% done with Concordance Index outputs"))
    counter <- counter + 1

    #Storing all of the c-index values in a vector that we can use later to build the plot
    c_finder <-current_cox$CV$index[1]
    current_c <- current_cox$CV$cvm[c_finder]
    current_c <- round(current_c, digits = 4)
    my_cindicies <- c(my_cindicies, current_c)
  }


  #write.csv(my_cindicies, file = paste0("three_weight_cindices_",x,"_mirna_",x,"_top.csv"))
  top_cindex <-max(my_cindicies)
  top_index <- which(my_cindicies==top_cindex)
  print(top_index)
  top_index_used <- top_index[1]


  #Now doing it for just the top c-index combo
  cox_models <- list()
  my_cindicies <- c()
  counter <- 1

  for (y in mad_sdes_mirna_optimized[top_index_used]) {
    current_weight <- y
    current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
                                    cox.predictors = current_weight,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    regular.cox = TRUE,
                                    save.regular.cox.genes = TRUE,
                                    my.filename = paste0("MiRNA/Global-heatmap/Top-cox-models/three_weight_optimized_global_search_mirnas_",x,"_",x,"_targets_top_index_",top_index_used,".RData"))

    cox_models[[as.character(counter)]] <- current_cox
    counter <- counter + 1

    #Storing all of the c-index values in a vector that we can use later to build the plot
    c_finder <-current_cox$CV$index[1]
    current_c <- current_cox$CV$cvm[c_finder]
    current_c <- round(current_c, digits = 4)
    my_cindicies <- c(my_cindicies, current_c)
    cox_models$`1`$CV
    c_index_df <- data.frame(c_index=my_cindicies)
    write.csv(c_index_df, file = paste0("MiRNA/Global-heatmap/Finished-outputs/cindices_global_heatmap_",x,"_","_mirna_",x,"_mirna_targets_top_index_",top_index_used,".csv"))
  }






}


#Optimizing the three metric model----
# mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
#                                                    second.metric = mirna.ranking,
#                                                    third.metric = sde.genes,
#                                                    my.filename = "three_weight_optimized_40_mirna_100_colorect.RData")

#Loading the merged data frame----
#For TCGA-COAD
#load("coad_df.RData", verbose = TRUE)
# 
# calculated_days <- merged_df$days.to.death - merged_df$days.to.last.follow.up
# calculated_days[calculated_days==0]=1
# merged_df$days.to.last.follow.up <- ifelse(is.na(calculated_days), merged_df$days.to.last.follow.up, calculated_days)
# merged_df$days.to.last.follow.up <- ifelse(merged_df$days.to.last.follow.up==0, 1, merged_df$days.to.last.follow.up)
# cox_time <- merged_df$days.to.last.follow.up
# cox_event <- merged_df$vital.status
# cox_tumor <- merged_df$ajcc_pathologic_stage
# cox_tumor_m <- merged_df$ajcc.pathologic.m
# cox_tumor_n <- merged_df$ajcc.pathologic.n
# cox_gender <- merged_df$gender
# cox_eth <- merged_df$ethnicity
# cox_race <- merged_df$race
# cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
# cox_df$days.to.last.follow.up <- cox_time
# cox_df$vital.status <- cox_event
# cox_df$tumor.stage <- cox_tumor
# cox_df$ajcc.m <- cox_tumor_m
# cox_df$ajcc.n <- cox_tumor_n
# cox_df$race <- cox_race
# cox_df$ethnicity <- cox_eth
# cox_df$gender <- cox_gender
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iv", replacement = 4)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iii", replacement = 3)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage ii", replacement = 2)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage i", replacement = 1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
# cox_df <- filter(cox_df, !tumor.stage=="not reported")
# cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]


#For TCGA-READ
#load("read_df.RData", verbose = TRUE)

#For TCGA-LUSC
#load("lusc_df.RData", verbose = TRUE)

#For TCGA-COAD + TCGA-READ
#load("coad_and_read_df.RData", verbose = TRUE)

#For TCGA-GBM
#load("gbm_df.RData", verbose = TRUE)

#Cox Model----
# load("three_weight_optimized_20_mirna_100.RData", verbose = TRUE)
# mad_sdes_mirna_optimized <- integrated_gene_lists
# 
# cox_models <- list()
# my_cindicies <- c()
# counter <- 1
# 
# for (x in mad_sdes_mirna_optimized[1:121]) {
#   current_weight <- x
#   current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
#                                   cox.predictors = current_weight,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   regular.cox = FALSE,
#                                   save.regular.cox.genes = FALSE,
#                                   my.filename = "three_weight_optimized_output_400_mirna_100_top.csv")
# 
#   cox_models[[as.character(counter)]] <- current_cox
#   perc_done <- (counter/121)*100
#   print(perc_done)
#   counter <- counter + 1
# 
#   #Storing all of the c-index values in a vector that we can use later to build the plot
#   c_finder <-current_cox$CV$index[1]
#   current_c <- current_cox$CV$cvm[c_finder]
#   current_c <- round(current_c, digits = 4)
#   my_cindicies <- c(my_cindicies, current_c)
# }
# 
# top_cindex <-max(my_cindicies)
# top_index <- which(my_cindicies==top_cindex)
# print(top_index)
# 
# 
# #Now doing it for just the top Concordance index result----
# cox_models <- list()
# my_cindicies <- c()
# counter <- 1
# 
# for (x in mad_sdes_mirna_optimized[top_index]) {
#   current_weight <- x
#   current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
#                                   cox.predictors = current_weight,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   regular.cox = FALSE,
#                                   save.regular.cox.genes = FALSE ,
#                                   my.filename = "three_weight_optimized_output_40_mirna_100_top_colorect.csv")
# 
#   cox_models[[as.character(counter)]] <- current_cox
#   counter <- counter + 1
# 
#   #Storing all of the c-index values in a vector that we can use later to build the plot
#   c_finder <-current_cox$CV$index[1]
#   current_c <- current_cox$CV$cvm[c_finder]
#   current_c <- round(current_c, digits = 4)
#   my_cindicies <- c(my_cindicies, current_c)
# }
# 
# cox_models$`1`$CV
# c_index_df <- data.frame(c_index=my_cindicies, active_coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
# write.csv(c_index_df, file = "CC_Singlecell_MMS_Mirna_800_10_Targets_colorect.csv")

#For plotting the coefficients----
# coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
# coef_df_sub <- filter(coef_df, abs(coefs)>0.1)
# 
# 
# 
# colnames(coef_df_sub)[2] <- "Gene"
# coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
#   geom_col()+
#   theme_bw()+
#   ggtitle("CC Singlecell MMS COAD Cox Coefficients")+
#   ylab("Coefficients")+
#   xlab("Gene")+
#   theme(legend.position = "none",
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
#   scale_y_continuous(expand = expansion(mult = c(0, .1)))+
#   coord_flip()
# 
# 
# #coef_plot
# ggsave("CC_Singlecell_MMS_Mirna_800_10_Targets.png", plot = coef_plot)

