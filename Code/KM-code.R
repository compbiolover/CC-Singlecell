#KM Curves----

#Loading needed packages----
library(survival)
library(survminer)

#These are examples I found of how to calculate the birinarizing of the expression data
#trainScore=apply(surv_gene_df,1,weight_function)
#risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))


#Saving just the active gene names
#For MAD metric currently
active_genes <-rownames(Coefficients)[Active.Index]
surv_gene_df=merged_df[,active_genes]
beta <-Active.Coefficients
gene_expr <- mean(as.matrix(surv_gene_df))
#patient_gene_expr_1 <- as.vector(apply(surv_gene_df,1,sum)) #Might be a problem here. Don't sum up. This is a confirmed problem
patient_gene_expr<- surv_gene_df[1:length(rownames(surv_gene_df)),]
subtracted_value <- patient_gene_expr - gene_expr
hr_calc <- beta*(subtracted_value)
#hr_calc_2 <- beta*(patient_gene_expr_1-gene_expr) #Testing original method. Original method confirmed to be faulty (-1195.41 calculated value for first patient. Should be a row with a score for each gene in this patient)
risk <- as.vector(apply(hr_calc,1,sum))
med_hr_value <- median(as.matrix(hr_calc))
risk <- ifelse(risk>med_hr_value,"high","low")
#risk_test <- ifelse(hr_calc>median(as.matrix(hr_calc)),"high", "low")
#risk <- as.vector(ifelse(hr_calc>median(as.matrix(hr_calc)),"high", "low"))
surv_gene_df <- cbind(risk, surv_gene_df)
vital.status <- merged_df$vital.status
days.to.last.follow.up <- merged_df$days.to.last.follow.up
surv_gene_df <- cbind(vital.status, surv_gene_df)
surv_gene_df <- cbind(days.to.last.follow.up, surv_gene_df)
km_fit <- survfit(Surv(days.to.last.follow.up, vital.status) ~ risk, data = surv_gene_df)

#P-value calculation for KM curves
diff=survdiff(Surv(days.to.last.follow.up, vital.status) ~risk,data = surv_gene_df)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue

#KM Curves plotting code
surPlot<-ggsurvplot(km_fit,
                    data=surv_gene_df,
                    pval=paste0("p=",pValue),
                    pval.size=4,
                    conf.int = FALSE,
                    legend.labs=c("High risk", "Low risk"),
                    legend.title="Risk",
                    xlab="Time (days)",
                    title="All 3 Metrics",
                    palette=c("red", "blue"))
