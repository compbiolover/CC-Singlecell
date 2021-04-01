#Name: km_curves.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently take in my cox model
#output and plot the Kaplan-Meir curves.

#Loading needed packages----
library(survival)
library(survminer)

#Doing the hazard ratio calculation----
hr_calculator <- function(model.coefs=Coefficients, data=merged_df, my.remove=NULL){
  require(survival)
  require(survminer)
  
  hr_return_list <- list()

  Active.Index <- which(as.logical(model.coefs) != 0)
  Active.Coefficients  <- model.coefs[Active.Index]
  active_genes <-rownames(model.coefs)[Active.Index]
  print(active_genes)
  active_genes <- active_genes[!active_genes %in% my.remove]
  surv_gene_df=data[,active_genes]
  beta <-Active.Coefficients
  gene_expr <- mean(as.matrix(surv_gene_df))
  # print(class(surv_gene_df))
  print(dim(surv_gene_df))
  patient_gene_expr<- surv_gene_df[1:length(rownames(surv_gene_df)),]
  subtracted_value <- patient_gene_expr - gene_expr
  hr_calc <- beta*(subtracted_value)
  risk <- as.vector(apply(hr_calc,1,sum))
  med_hr_value <- median(as.matrix(hr_calc))
  risk <- ifelse(risk>med_hr_value, "high", "low")
  surv_gene_df <- cbind(risk, surv_gene_df)
  vital.status <- data$vital.status
  days.to.last.follow.up <- data$days.to.last.follow.up
  surv_gene_df <- cbind(vital.status, surv_gene_df)
  surv_gene_df <- cbind(days.to.last.follow.up, surv_gene_df)
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
