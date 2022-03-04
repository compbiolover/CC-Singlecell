#Name: auc_calculator.R
#Purpose: Develop a score based on our gene signature and then testing it with
#AUC-ROC

#Load needed functions----
cox_model_fitter <- function(my.seed       =1,
                             cox.df        =NULL,
                             gene.num      =1800,
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
  require(survivalROC)
  require(survminer)
  require(tidyverse)
  
  #Setting the seed for reproducible----
  #output
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
  Coefficients2 <- coef(cv_fit, s = "lambda.min")
  
  
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Index2 <- which(as.logical(Coefficients2) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  Active.Coefficients2 <- Coefficients2[Active.Index2]
  active_genes <-rownames(Coefficients)[Active.Index]
  active_genes2 <- rownames(Coefficients2)[Active.Index2]
  
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
  
  #Calculating risk score from coefficients and expression levels----
  risk_genes <- cox.df[,active_genes]
  betas <- Active.Coefficients
  genes <- active_genes
  risk_df <- as.data.frame(rbind(genes, betas))
  colnames(risk_df) <- risk_df["genes",]
  risk_df <- risk_df[-c(1), ]
  risk_df_finished <- risk_df[rep(seq_len(nrow(risk_df)), each = nrow(cox.df)), ]
  #risk_df_finished <- risk_df[rep(seq_len(nrow(risk_df)), each = 501), ]
  rownames(risk_df_finished) <- rownames(risk_genes)
  risk_df_finished %<>% 
    mutate_each(funs(if(is.character(.)) as.numeric(.) else .))
  finished_auc_df <- risk_genes * risk_df_finished
  finished_auc_df$score <- rowSums(finished_auc_df)
  #finished_auc_df$median_score <- median(finished_auc_df$score)
  #finished_auc_df$score_group <- ifelse(finished_auc_df$score>finished_auc_df$median_score, "high", "low")
  
  #AUC-ROC curve
  auc_times <- seq(1, 1825, 50)
  auc_return <- list()
  counter <- 1
  for(t in auc_times){
    auc_roc_curve <- survivalROC(Stime = cox.df$days.to.last.follow.up, status = cox.df$vital.status, marker = finished_auc_df$score, predict.time = t, method = "KM")
    auc_return[[counter]] <- auc_roc_curve
    counter <- counter + 1
    
  }
  #auc_roc_curve <- survivalROC(Stime = cox.df$days.to.last.follow.up, status = cox.df$vital.status, marker = finished_auc_df$score, predict.time = 1825, method = "KM")
  
  #Adding the relevant data bits to list to return
  cox_data[["CV"]] <- cv_fit
  cox_data[["Coefficients"]] <- Coefficients
  cox_data[["Coefficients 2"]] <- Coefficients2
  cox_data[["Active Coefficients"]] <- Active.Coefficients
  cox_data[["Active Coefficients 2"]] <- Active.Coefficients2
  cox_data[["Active Index"]] <- Active.Index
  cox_data[["Active Index 2"]] <- Active.Index2
  cox_data[["Active Genes"]] <- active_genes
  cox_data[["Active Genes 2"]] <- active_genes2
  cox_data[["Predictors"]] <- my_predictors
  cox_data[["Predicted Survival"]] <- my_surv
  cox_data[["AUC-ROC Data"]] <- finished_auc_df
  cox_data[["AUC-ROC Curve"]] <- auc_return 
  
  
  #Returning our finished output----
  return(cox_data)
  
}

#risk_score_calculator-----
risk_score_calculator <- function(my.file="Data/Data-from-Cleaner-code/Regular_cox_model_outputs/coad_and_read_regular_cox_genes_1800.csv", 
                                  tumor.data=FALSE, n.data=FALSE, my.title="Finished KM Plot", cox.df=cox_df, show.pval=TRUE, show.pval.method=FALSE){
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
                         pval=show.pval,
                         pval.method = show.pval.method,
                         pval.size=pval.digits,
                         conf.int = confidence.int,
                         legend.labs=legend.labs,
                         legend.title=legend.title,
                         xlab=x.lab,
                         title=plot.title,
                         palette=color.pal,
                         font.main= c(40, "bold"),
                         font.x= c(40, "bold"),
                         font.y=c(40, "bold"),
                         font.tickslab=c(40, "plain"),
                         font.legend=c(40, "plain"))
    
    
    #Returning our finished KM plot----
    return(sur_Plot)
  }
  
  
  #List to store the return objects in
  survival_return <- vector(mode = "list", length = 2)
  #Read in the data----
  my_file <- read.csv(my.file)
  colnames(my_file)[c(1,6)] <- c("gene","p.value")
  
  if(tumor.data==FALSE & n.data==FALSE){
    risk_df <- cox.df[,my_file$gene]
    risk_df <- as.matrix(risk_df)
    print(dim(risk_df))
    gene_sign <- ifelse(my_file$coef>0, 1,-1)
    print(length(gene_sign))
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
    #Making a data frame of the converted scores
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
    finished_plot <- km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    finished_plot
    
    survival_return[["KM Plot"]] <- finished_plot
  }
  
  #Only for if there is tumor data included----
  if(tumor.data==TRUE & n.data==FALSE){
    tumor_info <- filter(my_file, gene=="tumor.stage1" | gene=="tumor.stage2" | gene=="tumor.stage3" | gene=="tumor.stage4")
    tumor_names <- tumor_info$gene
    print(tumor_names)
    tumor_info <- length(tumor_info$gene)
    just_genes <- length(my_file$gene) - tumor_info
    if(just_genes==1){
      risk_df <- cox.df[,my_file$gene[1]]
      risk_df <- as.data.frame(risk_df)
    }else{
      risk_df <- cox.df[,my_file$gene[1:just_genes]]
    }
    
    
    for(x in tumor_names){
      if(x=="tumor.stage1"){
        risk_df$tumorstage1 <- ifelse(cox.df$tumor.stage==1, 5,0)
      }else if(x=="tumor.stage2"){
        risk_df$tumorstage2 <- ifelse(cox.df$tumor.stage==2, 10,0)
      }else if(x=="tumor.stage3"){
        risk_df$tumorstage3 <- ifelse(cox.df$tumor.stage==3, 15,0)
      }else if(x=="tumor.stage4"){
        risk_df$tumorstage4 <- ifelse(cox.df$tumor.stage==4, 30,0)
      }
    }
    risk_df <- as.matrix(risk_df)
    gene_sign <- ifelse(my_file$coef>0, 1, -1)
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
    #Making a data frame of the converted scores
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
    finished_plot <-km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    finished_plot
    survival_return[["KM Plot"]] <- finished_plot
    
  }else if(tumor.data==TRUE & n.data==TRUE){
    tumor_info <- filter(my_file, gene=="tumor.stage1" | gene=="tumor.stage2" | gene=="tumor.stage3" | gene=="tumor.stage4")
    n_info <- filter(my_file, gene=="ajcc.n1" | gene=="ajcc.n2" | gene=="ajcc.n3")
    tumor_names <- tumor_info$gene
    n_names <- n_info$gene
    print(tumor_names)
    print(n_names)
    tumor_info <- length(tumor_info$gene)
    n_info <- length(n_info$gene)
    just_genes <- length(my_file$gene) - tumor_info - n_info
    if(just_genes==1){
      risk_df <- cox.df[,my_file$gene[1]]
      risk_df <- as.data.frame(risk_df)
    }else{
      risk_df <- cox.df[,my_file$gene[1:just_genes]]
    }
    
    for(x in tumor_names){
      if(x=="tumor.stage1"){
        risk_df$tumorstage1 <- ifelse(cox.df$tumor.stage==1, 5,0)
      }else if(x=="tumor.stage2"){
        risk_df$tumorstage2 <- ifelse(cox.df$tumor.stage==2, 10,0)
      }else if(x=="tumor.stage3"){
        risk_df$tumorstage3 <- ifelse(cox.df$tumor.stage==3, 15,0)
      }else if(x=="tumor.stage4"){
        risk_df$tumorstage4 <- ifelse(cox.df$tumor.stage==4, 30,0)
      }
    }
    
    for(x in n_names){
      if(x=="ajcc.n0"){
        risk_df$ajcc.n0 <- ifelse(cox.df$ajcc.n==0, 5,0)
      }else if(x=="ajcc.n1"){
        risk_df$ajcc.n1 <- ifelse(cox.df$ajcc.n==1, 10,0)
      }else if(x=="ajcc.n2"){
        risk_df$ajcc.n2 <- ifelse(cox.df$ajcc.n==2, 15,0)
      }else if(x=="ajcc.n3"){
        risk_df$ajcc.n3 <- ifelse(cox.df$ajcc.n==3, 20,0)
      }
    }
    
    
    risk_df <- as.matrix(risk_df)
    #gene_sign <- ifelse(my_file$coef>0, 1, -1)
    
    print(length(gene_sign))
    print(dim(risk_df))
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
    #Making a data frame of the converted scores
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
    finished_plot <-km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    survival_return[["KM Plot"]] <- finished_plot
  }
  
  survival_return[["Survival DF"]] <- converted_df
  return(survival_return)
}







#Loading needed libraries
library(survivalROC)

#AUC-ROC plot
plot(cox_models$`1`$`AUC-ROC Curve`$FP, cox_models$`1`$`AUC-ROC Curve`$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(cox_models$`1`$`AUC-ROC Curve`$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
