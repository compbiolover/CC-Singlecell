#Name: cox_model.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Fit an arbitrarily large number of predictors (genes) to a
#cox model to predict survival time


cox_model_fitter <- function(my.seed       = 1,
                             my.alpha      = 1,
                             my.dataset    = "COAD",
                             cox.df        =NULL,
                             gene.num      =1800,
                             cox.predictors=NULL,
                             tumor.stage   =FALSE,
                             tumor.n       =FALSE,
                             tumor.m       =FALSE,
                             n.folds       =10,
                             remove.stage   = c("tumor.stage1","tumor.stage2",
                                                "tumor.stage3", "tumor.stage4"),
                             remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2",
                                                "ajcc.n3"),
                             my.filename   ="~/Desktop/my_models_active_coefs.csv",
                             verbose = TRUE){
  
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
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(doParallel))
  suppressMessages(library(glmnet))
  suppressMessages(library(lmtest))
  suppressMessages(library(parallel))
  suppressMessages(library(survival))
  suppressMessages(library(survminer))
  
  #Setting the seed for reproducible output----
  set.seed(my.seed)
  
  
  #Setting the number of processors on the----
  #machine to speed up the fitting
  num_of_cores <- parallel::detectCores()
  registerDoParallel(cores = num_of_cores)
  
  #Making the list that will store the----
  #data we will return
  cox_data <- vector(mode = "list", length = 12)
  
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
  my_predictors <- paste("~", paste(my_predictors[1:length(my_predictors)], 
                                    collapse = "+"))
  
  
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
    if(verbose==TRUE){
      print("This is the genes only predictor....")
    }
  }
  
  my_x <- model.matrix(my_predictors, cox.df)
  
  
  
  #The response object for the cox model----
  my_y <- Surv(time = cox.df$days.to.last.follow.up, event = cox.df$vital.status)
  
  my_foldid<-sample(1:10,size=length(my_y),replace=TRUE)
  
  # if(my.dataset=="COAD"){
  #   my_foldid <- readRDS(file = "coad_df_fold_id.rds")
  # }else if (my.dataset=="READ"){
  #   my_foldid <- readRDS(file = "read_df_fold_id.rds")
  # }else if(my.dataset=="GBM"){
  #   my_foldid <- readRDS(file = "gbm_df_fold_id.rds")
  # }
  
  #Saving the x and y matrices for use to do inference and confidence interval
  #construction later
  # if(save.x==TRUE){
  #   write.csv(my_x, file = x.filename)
  # }
  # 
  # if(save.y==TRUE){
  #   write.csv(my_y, file = y.filename)
  # }

  
  #The 10-fold cross-validation fit----
  cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = n.folds, type.measure = "C",
                      maxit=100000, family="cox", parallel = TRUE,
                      alpha = 1)
  
  
  #Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family =  "cox", maxit = 100000)
  
  Coefficients <- coef(fit, s = cv_fit$lambda.min)
  
  
  
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  active_genes <-rownames(Coefficients)[Active.Index]
 
  
  #Saving the coefficients of the model
  active_coefs_df <- cbind(active_genes, Active.Coefficients)
  write.csv(active_coefs_df, file = my.filename)
  
  #Assessing the performance of the 10-fold cross-validation
  #on the entire data set----
  # model_perf<- assess.glmnet(cv_fit, newx = my_x, newy = my_y)
  
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
  

  
  #Adding the relevant data bits to list to return
  cox_data[["CV"]] <- cv_fit
  #cox_data[["Cox Performance"]] <- model_perf
  cox_data[["Coefficients"]] <- Coefficients
  cox_data[["Active Coefficients"]] <- Active.Coefficients
  cox_data[["Active Index"]] <- Active.Index
  cox_data[["Active Genes"]] <- active_genes
  cox_data[["Predictors"]] <- my_predictors
  cox_data[["Finished Coefficients"]] <- active_coefs_df
  
  
  #Returning our finished output----
  return(cox_data)
  
  }
  
 




