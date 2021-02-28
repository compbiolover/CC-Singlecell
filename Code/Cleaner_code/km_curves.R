#Name: km_curves.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently take in my cox model
#output and plot the Kaplan-Meir curves.

#Loading needed packages----
library(survival)
library(survminer)

#Doing the hazard ratio calculation----
hr_calculator <- function(model.coefs=Coefficients, data=merged_df){
  hr_return_list <- list()
  
  Active.Index <- which(as.logical(model.coefs) != 0)
  Active.Coefficients  <- model.coefs[Active.Index]
  active_genes <-rownames(model.coefs)[Active.Index]
  surv_gene_df=merged_df[,active_genes]
  beta <-Active.Coefficients
  gene_expr <- mean(as.matrix(surv_gene_df))
  patient_gene_expr<- surv_gene_df[1:length(rownames(surv_gene_df)),]
  subtracted_value <- patient_gene_expr - gene_expr
  hr_calc <- beta*(subtracted_value)
  risk <- as.vector(apply(hr_calc,1,sum))
  med_hr_value <- median(as.matrix(hr_calc))
  risk <- ifelse(risk>med_hr_value,"high","low")
  surv_gene_df <- cbind(risk, surv_gene_df)
  vital.status <- merged_df$vital.status
  days.to.last.follow.up <- merged_df$days.to.last.follow.up
  surv_gene_df <- cbind(vital.status, surv_gene_df)
  surv_gene_df <- cbind(days.to.last.follow.up, surv_gene_df)
  km_fit <- survfit(Surv(days.to.last.follow.up, vital.status) ~ risk, data = surv_gene_df)
  
  hr_return_list[["DF"]] <- surv_gene_df
  hr_return_list[["KM"]] <- km_fit
  return(hr_return_list)
}

#Only for simulation data----
alive_patients <- filter(merged_df, vital.status == 0)
dead_patients <- filter(merged_df, vital.status == 1)


surv_gene_df_alive <- filter(surv_gene_df, vital.status==0)
surv_gene_df_dead <- filter(surv_gene_df, vital.status==1)
surv_gene_df_alive$days.to.last.follow.up <- NULL
surv_gene_df_alive$vital.status <- NULL
positive_multi <- runif(42, min = 1000, max = 10000)
negative_multi <- runif(42, min = 0, max = 0.1)
surv_gene_df_alive <- surv_gene_df_alive * positive_multi
surv_gene_df_dead$days.to.last.follow.up <- NULL
surv_gene_df_dead$vital.status <- NULL
surv_gene_df_dead <- surv_gene_df_dead * negative_multi
surv_gene_df_remade <- surv_gene_df_alive %>% full_join(surv_gene_df_dead, by = c(colnames(surv_gene_df)[3:length(colnames(surv_gene_df))]))
surv_gene_df <- surv_gene_df_remade


# #Setting the value of the days.to.last.follow.up column based on the value of the vital.status column
# merged_df$days.to.last.follow.up[merged_df$vital.status==0]=runif(length(alive_patients$vital.status), min = 5000, max = 10000)
# merged_df$days.to.last.follow.up[merged_df$vital.status==1]=runif(length(dead_patients$vital.status), min = 0, max = 0.1)



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
