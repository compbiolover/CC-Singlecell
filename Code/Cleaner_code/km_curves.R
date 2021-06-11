#Name: km_curves.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently take in my cox model
#output and plot the Kaplan-Meir curves.

#Doing the hazard ratio calculation----
hr_calculator <- function(model.coefs            =Coefficients, 
                          data                   =merged_df, 
                          my.remove              =c("tumor.stagestge i","tumor.stagestge ii", "tumor.stagestge iii", "tumor.stagestge iv", "ajcc.nN0","ajcc.nN1", "ajcc.nN2", "ajcc.nN3", "ajcc.nNX"), 
                          include.cat.data       =TRUE, 
                          tumor.stage            =TRUE, 
                          n.stage                =TRUE,
                          early.tumor.stage.multi=10,
                          early.n.stage.multi    =10,
                          late.tumor.stage.multi =10000, 
                          late.n.stage.multi     =10000){
  require(tidyverse)
  require(survival)
  require(survminer)
  
  hr_return_list <- list()

  
  
  
  if (include.cat.data==TRUE){
    print("Categorical data included in HR calculation.")
    if(tumor.stage==TRUE & n.stage==FALSE){
      print("Tumor stage data included in HR calculation.")
      data$tumor.stage <- gsub(data$tumor.stage, pattern="stge iv", replacement=4)
      data$tumor.stage <- gsub(data$tumor.stage, pattern="stge iii", replacement=3)
      data$tumor.stage <- gsub(data$tumor.stage, pattern="stge ii", replacement=2)
      data$tumor.stage <- gsub(data$tumor.stage, pattern="stge i", replacement=1)
      data$tumor.stage <- gsub(data$tumor.stage, pattern=4, replacement=4*late.tumor.stage.multi)
      data$tumor.stage <- gsub(data$tumor.stage, pattern=3, replacement=3*late.tumor.stage.multi)
      data$tumor.stage <- gsub(data$tumor.stage, pattern=2, replacement=2*early.tumor.stage.multi)
      data$tumor.stage <- gsub(data$tumor.stage, pattern=1, replacement=1*early.tumor.stage.multi)
      
      Active.Index <- which(as.logical(model.coefs) != 0)
      Active.Coefficients  <- model.coefs[Active.Index]
      active_genes <-rownames(model.coefs)[Active.Index]
      print("Active length genes before....")
      print(length(active_genes))
      print("Active genes before removal....")
      print(active_genes)
      active_genes <- active_genes[!active_genes %in% my.remove]
      print("Active length genes after....")
      print(length(active_genes))
      print("Active genes after removal....")
      print(active_genes)
      surv_gene_df <- data[,active_genes]
      gene_expr <- mean(as.matrix(surv_gene_df))
      if(is.null(dim(surv_gene_df))==TRUE){
        surv_gene_df <- as.data.frame(surv_gene_df)
      }
      patient_gene_expr<- surv_gene_df[1:length(rownames(surv_gene_df)),]
      rm(surv_gene_df)
      active_genes <- c(active_genes, "tumor.stage")
      print(active_genes)
      #####Check this also!!!
      
      if(tumor.stage==TRUE & n.stage==TRUE){
        print("N stage data included in HR calculation.")
        data$ajcc.n <- gsub(data$ajcc.n, pattern="N0", replacement=1)
        data$ajcc.n <- gsub(data$ajcc.n, pattern="N1", replacement=2)
        data$ajcc.n <- gsub(data$ajcc.n, pattern="N2", replacement=3)
        data$ajcc.n <- gsub(data$ajcc.n, pattern="N3", replacement=4)
        data$ajcc.n <- gsub(data$ajcc.n, pattern=4, replacement=4*late.n.stage.multi)
        data$ajcc.n <- gsub(data$ajcc.n, pattern=3, replacement=3*late.n.stage.multi)
        data$ajcc.n <- gsub(data$ajcc.n, pattern=2, replacement=2*early.n.stage.multi)
        data$ajcc.n <- gsub(data$ajcc.n, pattern=1, replacement=1*early.n.stage.multi)
        
        Active.Index <- which(as.logical(model.coefs) != 0)
        Active.Coefficients  <- model.coefs[Active.Index]
        active_genes <-rownames(model.coefs)[Active.Index]
        active_genes <- active_genes[!active_genes %in% my.remove]
        surv_gene_df <- data[,active_genes]
        gene_expr <- mean(as.matrix(surv_gene_df))
        if(is.null(dim(surv_gene_df))==TRUE){
          surv_gene_df <- as.data.frame(surv_gene_df)
          #View(surv_gene_df)
        }
        patient_gene_expr<- surv_gene_df[1:length(rownames(surv_gene_df)),]
        rm(surv_gene_df)
        active_genes <- c(active_genes, "ajcc.n", "tumor.stage")
        print(active_genes)
        
      }else{
        print("N stage data is NOT included in HR calculation.")
      }
      
      
    }else{
      print("Tumor stage data is NOT included in HR calculation.")
    }
    
    
  }else{
    Active.Index <- which(as.logical(model.coefs) != 0)
    Active.Coefficients  <- model.coefs[Active.Index]
    active_genes <-rownames(model.coefs)[Active.Index]
    active_genes <- active_genes[!active_genes %in% my.remove]
  }
  
  surv_gene_df=data[,active_genes]
  if(is.null(dim(surv_gene_df))==TRUE){
    surv_gene_df <- as.data.frame(surv_gene_df)
  }
  surv_gene_df <- apply(surv_gene_df, c(1,2), as.numeric)
  #View(surv_gene_df)
  beta <-Active.Coefficients
  gene_expr <- mean(as.matrix(surv_gene_df))
  #print(class(surv_gene_df))
  print(dim(surv_gene_df))
  patient_gene_expr<- surv_gene_df[1:length(rownames(surv_gene_df)),]
  subtracted_value <- patient_gene_expr - gene_expr
  #print(dim(subtracted_value))
  print(length(beta))
  print(beta)
  hr_calc <- beta*(subtracted_value)
  #print(hr_calc)
  if(is.null(dim(hr_calc))){
    risk <- hr_calc
  }else{
    risk <- as.vector(apply(hr_calc,1,sum))
  }
  
  med_hr_value <- median(as.matrix(hr_calc))
  risk <- ifelse(risk>med_hr_value, "high", "low")
  surv_gene_df <- cbind(risk, surv_gene_df)
  vital.status <- data$vital.status
  days.to.last.follow.up <- data$days.to.last.follow.up
  surv_gene_df <- cbind(vital.status, surv_gene_df)
  surv_gene_df <- cbind(days.to.last.follow.up, surv_gene_df)
  surv_gene_df <- as.data.frame(surv_gene_df)
  #View(surv_gene_df)
  surv_gene_df$days.to.last.follow.up <- as.numeric(surv_gene_df$days.to.last.follow.up)
  surv_gene_df$vital.status <- as.numeric(surv_gene_df$vital.status)
  km_fit <- survfit(Surv(days.to.last.follow.up, vital.status) ~ risk, data = surv_gene_df)
  
  hr_return_list[["DF"]] <- surv_gene_df
  hr_return_list[["KM"]] <- km_fit
  hr_return_list[["beta"]] <- beta
  hr_return_list[["hr_calc"]] <- hr_calc
  hr_return_list[["risk"]] <- risk
  hr_return_list[["patient_gene_exp"]] <- patient_gene_expr
  hr_return_list[["subtracted_value"]] <- subtracted_value
  return(hr_return_list)
}

#KM p-value calculator----
km_pvalue_calculator <- function(surv.time         =time,
                                 surv.status       =status,
                                 surv.predictors   =risk,
                                 surv.df           =dataframe,
                                 num.sig.figs      =4,
                                 scientific.p.value=TRUE){
  
  #Loading required packages for the function.---- 
  #If they aren't already installed
  #they are installed automatically.
  require(survival)
  
  #First checking to make sure we have all the inputs----
  if(missing(surv.time)){
    stop("You must specify a time variable for the analysis.")
  }
  
  if(missing(surv.status)){
    stop("You must specify a status variable for the analysis.")
  }
  
  if(missing(surv.predictors)){
    stop("You must specify predictors for the analysis.")
  }
  
  if(missing(surv.df)){
    stop("You must specify a data.frame that contains the variables to carry out the analysis.")
  }
  
  if(missing(num.sig.figs)){
    stop("You need to specifiy a number of digits for the p-value.")
  }
  
  if(missing(scientific.p.value)){
    stop("You need to specifiy if you want the p-value in scientific notation.")
  }
  
  
  #P-value calculation for KM curves----
  colname_changes <- sapply(colnames(surv.df), gsub, pattern="-",replacement=".")
  colname_changes <- sapply(colnames(surv.df), gsub, pattern="_",replacement=".")
  colname_changes <- sapply(colnames(surv.df), gsub, pattern="/",replacement=".")
  colname_changes <- unlist(colname_changes)
  colnames(surv.df) <- colname_changes
  diff=survdiff(Surv(surv.time, surv.status) ~ surv.predictors, data = surv.df)
  p_Value=1-pchisq(diff$chisq,df=1)
  p_Value=signif(p_Value,num.sig.figs)
  p_Value=format(p_Value, scientific = scientific.p.value)
  
  #Returning the p_Value output----
  return(p_Value)
}

#KM plotter----
km_plotter <- function(km.fit        =km_fit,
                       data.source   =surv_gene_df,
                       p.value       =p_Value,
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
                       pval=paste0("p=",p.value),
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
