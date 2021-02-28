#Name: function_tester.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: A simple script to test all of the
#cleaner functions that I have made in other
#files

#Setting the working directory----
#on the server
setwd("/home/awillems/Projects/CC_Singlecell")

#Copying in the functions that----
#are helpful from Li et al.
#Function used to rank my linear model with just two metrics
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

#Function used to rank my linear model with three metrics
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
#three_weight_optimizer----
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
  integrated_gene_lists <- vector(mode = "list", length = my.list.length)
  
  #Doing the grid search
  for (x in a3_weights){
    for (y in weights) {
      current_ranking <- three_metric_geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3= y)
      current_ranking <- as.data.frame(current_ranking)
      integrated_gene_lists[[as.character(df_index)]] <- current_ranking
      df_index <- df_index + 1
    }
  }
  
  #Saving the output of the grid search to a .RData
  save(integrated_gene_lists, file = my.filename)
  return(integrated_gene_lists)
}
# weights <- seq(from = 0, to=1, by=0.1)
# a3_weights <- seq(from = 0, to=1, by=0.1)
# a3 <- 0
# df_index <- 1
# integrated_gene_lists <- vector(mode = "list", length = length(weights)*length(a3_weights))
# 
# for (x in a3_weights){
#   for (y in weights) {
#     current_ranking <- three_metric_geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3= y)
#     current_ranking <- as.data.frame(current_ranking)
#     integrated_gene_lists[[df_index]] <- current_ranking
#     df_index <- df_index + 1
#   }
# }
# 
# save(integrated_gene_lists, file = "Data/Exported_data/R_objects/integrated_gene_lists_for_all_three_metricsall_intersections_all_three_metrics_1800_gene_subset.RData")



#two_weight_optimizer----
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
  integrated_gene_lists <- vector(mode = "list", length = my.list.length)
  
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

#Loading the data files needed----
#for three_weight_optimizer test
load("Data/Exported_data/R_objects/mad.ranking.RData")
load("Data/Exported_data/R_objects/vim.sdes.ranking.RData")
load("Data/Exported_data/R_objects/mirna.ranking.RData")

#Test of three_weight_optimizer----
# test <- three_weight_optimizer(my.start       = 0,
#                                my.finish      = 1,
#                                step.size      = 0.1,
#                                a3.start       = 0,
#                                my.index       = 1,
#                                my.list.length = 121,
#                                first.metric   = mad.ranking,
#                                second.metric  = vim.sdes.ranking,
#                                third.metric   = mirna.ranking,
#                                my.filename    = "Data/Exported_data/R_objects/Function_tests/three_metric_optimization.RData")

#Test of two_weight_optimizer----
test <- two_weight_optimizer(my.start         = 0,
                               my.finish      = 1,
                               step.size      = 0.1,
                               my.index       = 1,
                               my.list.length = 11,
                               first.metric   = mad.ranking,
                               second.metric  = vim.sdes.ranking,
                               my.filename    = "Data/Exported_data/R_objects/Function_tests/two_metric_optimization.RData")
