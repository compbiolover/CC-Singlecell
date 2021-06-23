#Name: Cox_model_function.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Fit an arbitrarily large number of predictors (genes) to a
#cox model to predict survival time

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
                             remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2"),
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
    print("This is the genes only predictor....")
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
    active_genes <- active_genes[!active_genes %in% remove.stage]
    active_genes <- c(active_genes, "tumor.stage")
  }
  
  if(tumor.stage==TRUE & tumor.n==TRUE & tumor.m==FALSE){
    active_genes <- active_genes[!active_genes %in% remove.n.stage]
    active_genes <- c(active_genes, "ajcc.n")
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



#risk_score_calculator-----
risk_score_calculator <- function(my.file="Data/Data-from-Cleaner-code/Regular_cox_model_outputs/coad_and_read_regular_cox_genes_1800.csv", tumor.data=FALSE, n.data=FALSE, my.title="Finished KM Plot", cox.df=cox_df){
  #Required packages----
  require(survival)
  
  #Required functions----
  risk_converter <-function(my.name=my_genes[1], my.data=risk_df, my.med.exp=gene_info){
    if(my.med.exp[counter]>0){
      my.data[,my.name]<- ifelse(my.data[,my.name]> my.med.exp[,my.name], 1, 0)
    }else{
      my.data[,my.name]<- ifelse(my.data[,my.name]< my.med.exp[,my.name], -1, 0)
    }
    return(my.data[,my.name])
  }
  km_plotter <- function(km.fit        =km_fit,
                         data.source   =surv_gene_df,
                         p.value       =TRUE,
                         pval.digits   =4,
                         confidence.int=FALSE,
                         legend.labs   =c("High risk", "Low risk"),
                         legend.title  ="Risk",
                         x.lab         ="Time (days)",
                         plot.title    ="MiRNA",
                         color.pal     =c("red", "blue")){
    
    #Loading the needed package. If not installed it is automatically done.----
    require(survminer)
    
    #KM Curves plotting code----
    sur_Plot<-ggsurvplot(km.fit,
                         data=data.source,
                         pval=TRUE,
                         pval.method = FALSE,
                         pval.size=pval.digits,
                         conf.int = confidence.int,
                         legend.labs=legend.labs,
                         legend.title=legend.title,
                         xlab=x.lab,
                         title=plot.title,
                         palette=color.pal)
    
    #Returning our finished KM plot----
    return(sur_Plot)
  }
  
  
  
  #Read in the data----
  my_file <- read.csv(my.file)
  colnames(my_file)[c(1,6)] <- c("gene","p.value")
  if(tumor.data==FALSE & n.data==FALSE){
    risk_df <- cox_df[,my_file$gene]
  }
  
  #Only for if there is tumor data included----
  if(tumor.data==TRUE & n.data==FALSE){
    tumor_info <- filter(my_file, gene=="tumor.stage1" | gene=="tumor.stage2" | gene=="tumor.stage3" | gene=="tumor.stage4")
    tumor_names <- tumor_info$gene
    print(tumor_names)
    tumor_info <- length(tumor_info$gene)
    just_genes <- length(my_file$gene) - tumor_info
    risk_df <- cox.df[,my_file$gene[1:just_genes]]
    for(x in tumor_names){
      if(x=="tumor.stage1"){
        risk_df$tumorstage1 <- ifelse(cox_df$tumor.stage==1, 5,0)
      }else if(x=="tumor.stage2"){
        risk_df$tumorstage2 <- ifelse(cox_df$tumor.stage==2, 10,0)
      }else if(x=="tumor.stage3"){
        risk_df$tumorstage3 <- ifelse(cox_df$tumor.stage==3, 15,0)
      }else if(x=="tumor.stage4"){
        risk_df$tumorstage4 <- ifelse(cox_df$tumor.stage==4, 30,0)
      }
    }
    #risk_df$tumorstage1 <- ifelse(cox_df$tumor.stage==1, 5,0)
    # risk_df$tumorstage2 <- ifelse(cox_df$tumor.stage==2, 10,0)
    # risk_df$tumorstage3 <- ifelse(cox_df$tumor.stage==3, 15,0)
    # risk_df$tumorstage4 <- ifelse(cox_df$tumor.stage==4, 30,0)
    risk_df <- as.matrix(risk_df)
    gene_sign <- ifelse(my_file$coef>0, 1, -1)
    print(length(gene_sign))
    print(dim(risk_df))
    View(risk_df)
    risk_df <- risk_df%*%diag(gene_sign)
    risk_df <- as.data.frame(risk_df)
    colnames(risk_df) <- my_file$gene
    gene_info <- data.frame(med_expression=apply(risk_df, 2, median))
    gene_info <- t(gene_info)
    gene_info <- as.data.frame(gene_info)
    my_genes <- colnames(risk_df[1:length(colnames(risk_df))])
    risk_df$vital.status <- cox.df$vital.status
    risk_df$time <- cox.df$days.to.last.follow.up
    counter <- 1
    my_converted_scores <- list()
    total_risk_length <- length(colnames(risk_df))
    risk_length <- total_risk_length-2
    for(x in my_genes){
      current_risk <- risk_converter(my.name = my_genes[counter], my.data = risk_df[1:risk_length], my.med.exp = gene_info[1:risk_length])
      my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
      counter <- counter + 1
    }
    #Making a dataframe of the converted scores
    converted_df <- data.frame(my_df=1:dim(cox.df)[1])
    counter <- 1
    for (x in my_converted_scores) {
      converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
      counter <- counter + 1
    }
    
    converted_df[,"my_df"] <- NULL
    rownames(converted_df) <- rownames(cox.df)
    converted_df$vital.status <- cox.df$vital.status
    converted_df$time <- cox.df$days.to.last.follow.up
    converted_df <- apply(converted_df, 2, as.numeric)
    patient_risks <- rowSums(x=converted_df[,1:risk_length])
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- patient_risks
    converted_df <- apply(converted_df, 2, as.numeric)
    median_risk <- median(abs(converted_df[,1:risk_length]))
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- ifelse(converted_df$risk>median_risk, "high", "low")
    km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)
    
    km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    
    
  #For just the genes only----
  }else{
    risk_df <- cox.df[,my_file$gene]
    risk_df <- as.matrix(risk_df)
    gene_sign <- ifelse(my_file$coef>0, 1,-1)
    risk_df <- risk_df%*%diag(gene_sign)
    risk_df <- as.data.frame(risk_df)
    colnames(risk_df) <- my_file$gene
    gene_info <- data.frame(med_expression=apply(risk_df, 2, median))
    gene_info <- t(gene_info)
    gene_info <- as.data.frame(gene_info)
    my_genes <- colnames(risk_df[1:length(colnames(risk_df))])
    risk_df$vital.status <- cox.df$vital.status
    risk_df$time <- cox.df$days.to.last.follow.up
    counter <- 1
    my_converted_scores <- list()
    total_risk_length <- length(colnames(risk_df))
    risk_length <- total_risk_length-2
    for(x in my_genes){
      current_risk <- risk_converter(my.name = my_genes[counter], my.data = risk_df[1:risk_length], my.med.exp = gene_info[1:risk_length])
      my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
      counter <- counter + 1
    }
    #Making a dataframe of the converted scores
    converted_df <- data.frame(my_df=1:dim(cox.df)[1])
    counter <- 1
    for (x in my_converted_scores) {
      converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
      counter <- counter + 1
    }
    
    converted_df[,"my_df"] <- NULL
    rownames(converted_df) <- rownames(cox.df)
    converted_df$vital.status <- cox.df$vital.status
    converted_df$time <- cox.df$days.to.last.follow.up
    converted_df <- apply(converted_df, 2, as.numeric)
    patient_risks <- rowSums(x=converted_df[,1:risk_length])
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- patient_risks
    converted_df <- apply(converted_df, 2, as.numeric)
    median_risk <- median(abs(converted_df[,1:risk_length]))
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- ifelse(converted_df$risk>median_risk, "high", "low")
    km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)
    
    km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    
  }
  
  
}




# require(survival)
# my_file <- read.csv("Data/Data-from-Cleaner-code/Regular_cox_model_outputs/coad_and_read_regular_cox_genes_tumor_1800.csv")
# colnames(my_file)[c(1,6)] <- c("gene","p.value")
# 
# risk_df <- cox_df[,my_file$gene[1:240]]
# risk_df$tumorstage2 <- ifelse(cox_df$tumor.stage==2, 10,0)
# risk_df$tumorstage3 <- ifelse(cox_df$tumor.stage==3, 15,0)
# risk_df$tumorstage4 <- ifelse(cox_df$tumor.stage==4, 30,0)
# risk_df <- as.matrix(risk_df)
# gene_sign <- ifelse(my_file$coef>0, 1, -1)
# risk_df <- risk_df%*%diag(gene_sign)
# risk_df <- as.data.frame(risk_df)
# colnames(risk_df) <- my_file$gene
# gene_info <- data.frame(med_expression=apply(risk_df, 2, median))
# gene_info <- t(gene_info)
# gene_info <- as.data.frame(gene_info)
# risk_df$vital.status <- cox_df$vital.status
# risk_df$time <- cox_df$days.to.last.follow.up
# 
# 
# 
# my_genes <- colnames(risk_df[1:243])
# 
# counter <- 1
# my_converted_scores <- list()
# for(x in my_genes){
#   current_risk <- risk_converter(my.name = my_genes[counter], my.data = risk_df[1:243], my.med.exp = gene_info[1:243])
#   my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
#   counter <- counter + 1
# }
# 
# #Making a dataframe of the converted scores
# converted_df <- data.frame(my_df=1:673)
# counter <- 1
# for (x in my_converted_scores) {
#   converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
#   counter <- counter + 1
# }
# 
# converted_df[,"my_df"] <- NULL
# rownames(converted_df) <- rownames(cox_df)
# converted_df$vital.status <- cox_df$vital.status
# converted_df$time <- cox_df$days.to.last.follow.up
# converted_df <- apply(converted_df, 2, as.numeric)
# patient_risks <- rowSums(x=converted_df[,1:243])
# converted_df <- as.data.frame(converted_df)
# converted_df$risk <- patient_risks
# converted_df <- apply(converted_df, 2, as.numeric)
# median_risk <- median(abs(converted_df[,1:243]))
# converted_df <- as.data.frame(converted_df)
# converted_df$risk <- ifelse(converted_df$risk>median_risk, "high", "low")
# km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)
# 
# km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = "CC Singlecell MS + Tumor Stage COAD + READ")
# 

#Code for just testing tumor stage and n pathological state----
my_predictors <- paste("~", paste("ajcc.n", sep = "+"), collapse = "+")
my_predictors <- as.formula(my_predictors)
my_x <- model.matrix(my_predictors, cox_df)
my_y <- Surv(time = cox_df$days.to.last.follow.up, event = cox_df$vital.status)

cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = 10, type.measure = "C", maxit=100000, family="cox", parallel = TRUE)


