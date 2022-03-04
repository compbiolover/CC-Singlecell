#Loading files----
# load("mad.RData", verbose = TRUE)
# load("sde.RData", verbose = TRUE)
#mirna.ranking <- readRDS(file = "/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Inputs/RData/mirna_test2.rds")
#load("~/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Inputs/RData/800-310-targets.RData", verbose = TRUE)


#Needed functions----
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
    gene_sign <- ifelse(my_file$coef>0, 1, -1)
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
    
  }
  return(finished_plot)
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





#Doing the actual optimization----
# combo_used <- c("100_10", "100_110", "100_210", "100_310", "100_410", "100_510", "100_610", "100_710", "100_810", "100_910", "100_1010")

# for(c in combo_used){
#   mirna.ranking <- readRDS(file = paste0("~/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Inputs/RData/",c,"_targets_mirna_ranking.rds"))
#   mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
#                                                      second.metric = mirna.ranking,
#                                                      third.metric = sde.genes,
#                                                      my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-GBM/Data/MiRNA/Global-heatmap/Weighted-Optimizations/Optimization_",c,"_targets.RData"))

#Actual code----
cox_df <- readRDS("read_df_finished_v2.rds")
files <- list.files("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Weighted-Optimizations")

for(f in files){
  cox_models <- list()
  my_cindicies <- c()
  top_cindices <- c()
  counter <- 1
  load(file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Weighted-Optimizations/",f), verbose = TRUE)
  
  mad_sdes_mirna_optimized <- integrated_gene_lists
  
  for (y in mad_sdes_mirna_optimized[1:121]) {
    current_weight <- y
    current_cox <- cox_model_fitter(my.seed = 1,
                                    cox.df = cox_df,
                                    gene.num = 1100,
                                    cox.predictors = current_weight,
                                    tumor.stage = FALSE,
                                    tumor.n = FALSE,
                                    tumor.m = FALSE,
                                    regular.cox = FALSE,
                                    save.regular.cox.genes = FALSE,
                                    my.filename = paste0("three_weight_optimized_global_search_",x,"_mirna_",x,"_top.csv"))
    
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
  
  
  top_cindex <-max(my_cindicies)
  top_index <- which(my_cindicies==top_cindex)
  print(top_index)
  print(my_cindicies)
  top_index_used <- top_index[1]
  print(my_cindicies[top_index_used])
  #print(top_index_used)
  
  top_cindices <- c(top_cindices, my_cindicies[top_index_used])
  
  # c_index_df <- data.frame(c_index=my_cindicies[top_index_used])
  # write.csv(c_index_df, file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-cindices/top_cindices_mirna_used_combo_",c,"_index_",top_index_used,"_changed_df.csv"))
  
  #Now doing it for just the top c-index combo
  # cox_models <- list()
  # my_cindicies <- c()
  # counter <- 1
  
  # for (y in mad_sdes_mirna_optimized[top_index_used]) {
  #   current_weight <- y
  #   current_cox <- cox_model_fitter(my.seed = 1,
  #                                   cox.df = cox_df,
  #                                   gene.num = 1800,
  #                                   cox.predictors = current_weight,
  #                                   tumor.stage = FALSE,
  #                                   tumor.n = FALSE,
  #                                   tumor.m = FALSE,
  #                                   regular.cox = FALSE,
  #                                   save.regular.cox.genes = FALSE,
  #                                   my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-coefficients/top_coefs_global_search_mirnas_",c,"_targets_mirna_top_index_",top_index_used,"_changed_df.csv"))
  #   
  #   cox_models[[as.character(counter)]] <- current_cox
  #   counter <- counter + 1
  #   
  #   #Storing all of the c-index values in a vector that we can use later to build the plot
  #   c_finder <-current_cox$CV$index[1]
  #   current_c <- current_cox$CV$cvm[c_finder]
  #   current_c <- round(current_c, digits = 4)
  #   my_cindicies <- c(my_cindicies, current_c)
  #   cox_models$`1`$CV
  #   #print(my_cindicies)
  #   c_index_df <- data.frame(c_index=my_cindicies)
  #   c_index_df
  #   write.csv(c_index_df, file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-cindices/top_cindices_mirna_used_combo_",c,"_index_",top_index_used,"_changed_df.csv"))
  # }
  # 
  
  
  
  
}

write.csv(top_cindices, file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-READ/Data/MiRNA/Global-heatmap/Outputs/Top-cindices/top_cindices_new_mirna.csv"))

  # cox_models <- list()
  # my_cindicies <- c()
  # counter <- 1
  # load(file = files[f], verbose = TRUE)
  # 
  # 
  # for (y in mad_sdes_mirna_optimized[1:121]) {
  #   current_weight <- y
  #   current_cox <- cox_model_fitter(my.seed = 1,
  #                                   cox.df = cox_df,
  #                                   gene.num = 1800,
  #                                   cox.predictors = current_weight,
  #                                   tumor.stage = FALSE,
  #                                   tumor.n = FALSE,
  #                                   tumor.m = FALSE,
  #                                   regular.cox = FALSE,
  #                                   save.regular.cox.genes = FALSE,
  #                                   my.filename = paste0("three_weight_optimized_global_search_",x,"_mirna_",x,"_top.csv"))
  # 
  #   cox_models[[as.character(counter)]] <- current_cox
  #   perc_done <- (counter/121)*100
  #   print(paste0(round(perc_done, digits = 2), "% done with Concordance Index outputs"))
  #   counter <- counter + 1
  # 
  #   #Storing all of the c-index values in a vector that we can use later to build the plot
  #   c_finder <-current_cox$CV$index[1]
  #   current_c <- current_cox$CV$cvm[c_finder]
  #   current_c <- round(current_c, digits = 4)
  #   my_cindicies <- c(my_cindicies, current_c)
  # }
  # 
  # 
  # top_cindex <-max(my_cindicies)
  # top_index <- which(my_cindicies==top_cindex)
  # print(top_index)
  # print(my_cindicies)
  # top_index_used <- top_index[1]
  # print(my_cindicies[top_index_used])
  # #print(top_index_used)
  # 
  # #Now doing it for just the top c-index combo
  # cox_models <- list()
  # my_cindicies <- c()
  # counter <- 1
  # 
  # for (y in mad_sdes_mirna_optimized[top_index_used]) {
  #   current_weight <- y
  #   current_cox <- cox_model_fitter(my.seed = 1,
  #                                   cox.df = cox_df,
  #                                   gene.num = 1800,
  #                                   cox.predictors = current_weight,
  #                                   tumor.stage = FALSE,
  #                                   tumor.n = FALSE,
  #                                   tumor.m = FALSE,
  #                                   regular.cox = FALSE,
  #                                   save.regular.cox.genes = FALSE,
  #                                   my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-coefficients/top_coefs_global_search_mirnas_",c,"_targets_dbdemc_mirna_top_index_",top_index_used,".csv"))
  # 
  #   cox_models[[as.character(counter)]] <- current_cox
  #   counter <- counter + 1
  # 
  #   #Storing all of the c-index values in a vector that we can use later to build the plot
  #   c_finder <-current_cox$CV$index[1]
  #   current_c <- current_cox$CV$cvm[c_finder]
  #   current_c <- round(current_c, digits = 4)
  #   my_cindicies <- c(my_cindicies, current_c)
  #   cox_models$`1`$CV
  #   c_index_df <- data.frame(c_index=my_cindicies)
  #   write.csv(c_index_df, file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-cindices/top_cindices_dbdemc_mirna_used_combo_",c,"_index_",top_index_used,".csv"))
  # }
  # 
  # 






#Just individual steps with no for loop----
# mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
#                                                    second.metric = mirna.ranking,
#                                                    third.metric = sde.genes,
#                                                    my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Weighted-Optimizations/Optimization-",combo_used,"-targets.RData"))
# 
# 
# 
#Actual code
# cox_models <- list()
# my_cindicies <- c()
# counter <- 1
# load(file = "coad_df.RData", verbose = TRUE)
# 
# 
# for (y in mad_sdes_mirna_optimized[1:121]) {
#   current_weight <- y
#   current_cox <- cox_model_fitter(my.seed = 1,
#                                   cox.df = cox_df,
#                                   gene.num = 1800,
#                                   cox.predictors = current_weight,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   regular.cox = FALSE,
#                                   save.regular.cox.genes = FALSE,
#                                   my.filename = paste0("three_weight_optimized_global_search_",x,"_mirna_",x,"_top.csv"))
#   
#   cox_models[[as.character(counter)]] <- current_cox
#   perc_done <- (counter/121)*100
#   print(paste0(round(perc_done, digits = 2), "% done with Concordance Index outputs"))
#   counter <- counter + 1
#   
#   #Storing all of the c-index values in a vector that we can use later to build the plot
#   c_finder <-current_cox$CV$index[1]
#   current_c <- current_cox$CV$cvm[c_finder]
#   current_c <- round(current_c, digits = 4)
#   my_cindicies <- c(my_cindicies, current_c)
# }
# 
# 
# top_cindex <-max(my_cindicies)
# top_index <- which(my_cindicies==top_cindex)
# print(top_index)
# print(my_cindicies)
# top_index_used <- top_index[1]
# print(my_cindicies[top_index_used])
# #print(top_index_used)
# 
# #Now doing it for just the top c-index combo
# cox_models <- list()
# my_cindicies <- c()
# counter <- 1
# 
# for (y in mad_sdes_mirna_optimized[top_index_used]) {
#   current_weight <- y
#   current_cox <- cox_model_fitter(my.seed = 1,
#                                   cox.df = cox_df,
#                                   gene.num = 1800,
#                                   cox.predictors = current_weight,
#                                   tumor.stage = FALSE,
#                                   tumor.n = FALSE,
#                                   tumor.m = FALSE,
#                                   regular.cox = TRUE,
#                                   save.regular.cox.genes = TRUE,
#                                   my.filename = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-coeffecients/top_coefs_global_search_mirnas_",combo_used,"_targets_top_index_",top_index_used,".csv"))
#   
#   cox_models[[as.character(counter)]] <- current_cox
#   counter <- counter + 1
#   
#   #Storing all of the c-index values in a vector that we can use later to build the plot
#   c_finder <-current_cox$CV$index[1]
#   current_c <- current_cox$CV$cvm[c_finder]
#   current_c <- round(current_c, digits = 4)
#   my_cindicies <- c(my_cindicies, current_c)
#   cox_models$`1`$CV
#   c_index_df <- data.frame(c_index=my_cindicies)
#   write.csv(c_index_df, file = paste0("/home/awillems/Projects/CC_Singlecell/TCGA-COAD/Data/MiRNA/Global-heatmap/Outputs/Top-cindices/top_cindices_mirna_combo_",combo_used,"_index_",top_index_used,".csv"))
# }

#Individual metrics

#MAD/SDE
# load("read_df.RData", verbose = TRUE)
# 
# mad_res <- cox_model_fitter(my.seed = 1,
#                             cox.df = cox_df,
#                             gene.num = 1800,
#                             cox.predictors = integrated_gene_lists[[27]],
#                             tumor.stage = FALSE,
#                             tumor.n = FALSE,
#                             tumor.m = FALSE,
#                             regular.cox = TRUE,
#                             save.regular.cox.genes = TRUE,
#                             my.filename = "~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-COAD/top_peforming_coefs_0.704.csv")
# 
# 
# 
# c_finder <-mad_res$CV$index[1]
# current_c <- mad_res$CV$cvm[c_finder]
# current_c <- round(current_c, digits = 4)
# c_index_df <- data.frame(c_index=current_c)
# write.csv(c_index_df, file = "Data/TCGA-READ/MiRNA/Global_Search/Outputs/top_cindices_mirna.csv")
# 


