#gene_size_code_for_server.R
#Purpose: Exploit the speed of the lab server to fit a lot of cox models
#at different gene sizes

#Loading in needed dataframe----
load("cox_df.RData", verbose = TRUE)
dim(cox_df)

#These packages are required----
library(doParallel);packageVersion("doParallel")
library(glmnet);packageVersion("glmnet")
library(parallel);packageVersion("parallel")
library(survival);packageVersion("survival")
library(survminer);packageVersion("survminer")


#Loading in more needed data files----
load("mirna_mad_optimized.RData", verbose = TRUE)
mirna_mad_optimized <- integrated_gene_lists

load("sde_mirna_optimized.RData", verbose = TRUE)
mirna_sde_optimized <- integrated_gene_lists

#cox_model_fitter----
cox_model_fitter <- function(my.seed       =1,
                             cox.df        =NULL,
                             gene.num      =900,
                             cox.predictors=NULL,
                             tumor.stage   =FALSE,
                             tumor.m       =FALSE,
                             tumor.n       =FALSE){
  
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
  require(parallel)
  require(survival)
  
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
  if(tumor.stage==TRUE){
    my_predictors <- paste(my_predictors, "tumor.stage", sep = "+")
    if(tumor.m==TRUE){
      my_predictors <- paste(my_predictors,"ajcc.n", sep = "+")
      if(tumor.n==TRUE){
        my_predictors <- paste(my_predictors, "ajcc.n", sep = "+")
        my_predictors <- as.formula(my_predictors)
      }else{
        my_predictors <- as.formula(my_predictors)
      }
    }else{
      my_predictors <- as.formula(my_predictors)
    }
  }else{
    if(tumor.m==TRUE){
      my_predictors <- paste(my_predictors, "ajcc.m", sep = "+")
      my_predictors <- as.formula(my_predictors)
    }else{
      my_predictors <- as.formula(my_predictors)
    }
  }
  #my_predictors <- as.formula(my_predictors)
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
  
  
  
  #Adding the relevant data bits to list to return
  cox_data[["CV"]] <- cv_fit
  cox_data[["Coefficients"]] <- Coefficients
  cox_data[["Active Coefficients"]] <- Active.Coefficients
  cox_data[["Active Index"]] <- Active.Index
  cox_data[["Predictors"]] <- my_predictors
  
  #Returning our finished output----
  return(cox_data)
}

print("Doing Mirna + MAD Cox model fitting")

#Cox model for mirna + mad-----
cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1
gene_sizes <- seq(from=100, to=1800, by=50)
for (y in gene_sizes){
  print(y)
  for (x in mirna_mad_optimized) {
    current_weight <- x
    current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
    cox_models[[as.character(counter)]] <- current_cox
    counter <- counter + 1
    
    #Storing all of the c-index values in a vector that we can use later to build the plot
    c_finder <-current_cox$CV$index[1]
    current_c <- current_cox$CV$cvm[c_finder]
    current_c <- round(current_c, digits = 4)
    my_cindicies <- c(my_cindicies, current_c)
    my_gene_sizes <- c(my_gene_sizes, y)
    
  }
}


my_df <- data.frame(gene_num=my_gene_sizes, concordance_index=my_cindicies)
my_df$cindex_se <- sqrt(my_df$concordance_index/sqrt(length(my_df$concordance_index)))
write.csv(my_df, file = "mirna_mad_cc_patients_gene_size_data.csv")

df_to_plot <- data.frame(data=aggregate(x = my_df$concordance_index,              
                                        by = list(my_df$gene_num),              
                                        FUN = mean))    


colnames(df_to_plot) <- c("gene_num", "concordance_index")

print("Getting ready to save Mirna + MAD Plotting dataframe")

save(df_to_plot, file = "mirna_mad_cc_patients_gene_size_df_to_plot.RData")


print("Doing Mirna + SDES Cox model fitting")


#Optimal cox model code----
cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1
gene_sizes <- seq(from=100, to=1800, by=50)
current_weight <- mirna_sde_optimized$`3`

for (y in gene_sizes) {
  print(y)
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = y, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  my_gene_sizes <- c(my_gene_sizes, y)
  
}


#Cox model for mirna + sde-----
cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1
gene_sizes <- seq(from=100, to=1800, by=50)
for (y in gene_sizes){
  print(y)
  for (x in mirna_sde_optimized) {
    current_weight <- x
    current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = y, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
    cox_models[[as.character(counter)]] <- current_cox
    counter <- counter + 1
    
    #Storing all of the c-index values in a vector that we can use later to build the plot
    c_finder <-current_cox$CV$index[1]
    current_c <- current_cox$CV$cvm[c_finder]
    current_c <- round(current_c, digits = 4)
    my_cindicies <- c(my_cindicies, current_c)
    my_gene_sizes <- c(my_gene_sizes, y)
    
  }
}


my_df <- data.frame(gene_num=my_gene_sizes, concordance_index=my_cindicies)
my_df$cindex_se <- sqrt(my_df$concordance_index/sqrt(length(my_df$concordance_index)))
write.csv(my_df, file = "mirna_sde_cc_patients_gene_size_data.csv")

df_to_plot <- data.frame(data=aggregate(x = my_df$concordance_index,              
                                        by = list(my_df$gene_num),              
                                        FUN = mean))    


colnames(df_to_plot) <- c("gene_num", "concordance_index")

print("Getting ready to save Mirna + SDES Plotting dataframe")

save(df_to_plot, file = "mirna_sde_cc_patients_gene_size_df_to_plot.RData")


