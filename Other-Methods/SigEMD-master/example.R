## SigEMD needs log2(TPM+1) as input data.


setwd("/Users/peng/Google Drive/singleCellRNA/SigEMD") # set your working path
# SigEMD is based on R package "aod","arm","fdrtool","lars"
library(aod)
library(arm)
library(fdrtool)
library(lars)
source("FunImpute.R")
source("SigEMDHur.R")
source("SigEMDnonHur.R")
source("plot_sig.R")

# Load exampleData: SigEMD needs log2(TPM+1) as input data, which rows represents genes and columns represents samples. 
#"condition" is also need. 

load("exampleData.RData")
data <- dataclean(data)
databinary<- databin(data)
names(condition) <- colnames(data)

############################# Imputation #############################
## Imputation can be chosen by the user. It is not a necessary step in SigEMD, but we suggest to impute the missing zero values.
## A lasso regression model is fitted to impute the missing values.
## First, identify a set of genes need to be imputed (Hur_gene) and a set of genes (genes_use) used to fit the lasso model.
Hur_gene<- idfyImpgene(data,databinary,condition)
genes_use<- idfyUsegene(data,databinary,condition,ratio = 0.5) # ratio is default as 0.8, but it can be set by users to obtain appropriate number of genes_use.

datimp <- FunImpute(object = data, genes_use = (genes_use), genes_fit = (Hur_gene),dcorgene = NULL) 
data<-datimp$alldat 


############################# Idntify diferentially expressed genes #############################

# Call calculate_single. 
#Here we only use 5 permutations as an example, but in actual experiments using at least 100 permutations is advised. 
#If previours "Imputation" step is applied, then variable "Hur_gene" is can be direcly used. 
#If the user did not impute the data, then "Hur_gene" can be set to NULL. 
results<- calculate_single(data =  data,condition =  condition,Hur_gene = Hur_gene, binSize=0.2,nperm=5)
# Set Hur_gene to NULL if imputation isn't applied.
# results<- calculate_single(data =  data,condition =  condition,Hur_gene = NULL, binSize=0.2,nperm=5)

# The function returns a matrix with EMD score and the adjust p-value of each gene.
emd<- results$emdall
head(emd)

# The function plot_emd_density_sig will display the density distributions of each of the groups for a given gene.
plot_emd_density_sig(results,"PNKD")

# We can plot the gene with the largest EMD score:
plot_emd_density_sig(results,rownames(emd[order(-emd[,"emd"]),])[1])
