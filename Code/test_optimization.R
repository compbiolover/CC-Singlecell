#test_optimization.R
#Testing out ways to speed up optimization of the linear model

#Loading input lists-----
load("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-COAD/MAD/mad.RData")
load("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-COAD/SDE/sde.RData")

linear_weights <- seq(0,1, 0.1)

gene_names <- unique(c(names(sde.genes)), names(mad.genes))

lin_list <- rep(0, length(gene_names))

names(lin_list) <- gene_names
