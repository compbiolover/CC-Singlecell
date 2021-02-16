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
patient_gene_expr <- as.vector(apply(surv_gene_df,1,sum))
hr_calc <- beta*(patient_gene_expr-gene_expr)
risk <- as.vector(ifelse(hr_calc>median(hr_calc),"high", "low"))
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
                    palette=c("red", "blue"))
