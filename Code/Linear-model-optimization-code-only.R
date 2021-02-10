#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Contains just code needed to optimize linear model on
#lab server

#Loading needed packages----
library(tidyverse)
#Loading in the functions that I need----
#Ranking genes by three measurements
geneRank <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, a1 = 1,
                     a2 = 0, a3 = 0){
  gns = c(names(ranking1),names(ranking2),names(ranking3))
  gns = unique(gns)
  res = rep(0, length(gns))
  names(res) = gns
  for(i in names(res)) {
    res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2+getRank(ranking3, i)*a3
  }
  res
}

getRank <- function(ranking = NULL, gn = NULL){
  if (gn %in% names(ranking)) {
    return(ranking[gn])
  }
  else return(0.0)
}


#Setting the working directory----
setwd("/home/awillems/Projects/CC_Singlecell")
#Loading in the needed data files----
load(file = "Data/mad.ranking.RData")
load(file = "Data/vim.sdes.ranking.RData")
load(file = "Data/mirna.ranking.RData")

#Doing the optimization for the model----
set.seed(1)
weights <- seq(from = 0, to=1, by=0.1)
a3_weights <- seq(from = 0, to=1, by=0.1)
a3 <- 0
df_index <- 1
integrated_gene_lists_all_three_metrics <- vector(mode = "list", length =length(weights) * length(a3_weights))


for (x in a3_weights){
  for (y in weights) {
    current_ranking <- geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3=y)
    current_ranking <- as.data.frame(current_ranking)
    integrated_gene_lists_all_three_metrics[[df_index]] <- current_ranking
    df_index <- df_index + 1
  }
}

#Saving the output of the initial optimization
save(integrated_gene_lists_all_three_metrics, file = "Data/integrated_gene_lists_all_three_metrics.RData")

#Resetting these two variables for this new for loop
df_index <- 1

#Sorting all of the lists from greatest gene score to least gene score
for(x in integrated_gene_lists_all_three_metrics){
  current_list <- x
  current_list <- arrange(current_list, desc(current_ranking))
  integrated_gene_lists_all_three_metrics[[df_index]] <- current_list
  df_index <- df_index + 1
}

#Saving the sorted lists
save(integrated_gene_lists_all_three_metrics, file = "Data/sorted_gene_lists_all_three_metrics.RData")