#Name: model_optimizer.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Optimize the lists of genes from our
#integrated lists by tuning the weights of our
#linear model via grid search.

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



#Attempt at optimization-----
getIntersection <- function(gene.list1 = NULL,
                            gene.list2 = NULL,
                            gene.list3 = NULL){
  
  return_list <- vector(mode = "list", length = 3)
  
  common_genes <- intersect(c(names(gene.list1), names(gene.list2)), names(gene.list3))
  mad.genes.sub <- mad.genes[common_genes]
  sde.genes.sub <- sde.genes[common_genes]
  mirna.genes.sub <- mirna.ranking[common_genes]
  
  return_list[[1]] <- mad.genes.sub
  return_list[[2]] <- sde.genes.sub
  return_list[[3]] <- mirna.genes.sub
  return(return_list)
}




#Optimization for just 2 metrics----
two_weight_optimizer <- function(my.start        =0,
                                 my.finish       =1,
                                 step.size       =0.1,
                                 first.metric    =mad.genes,
                                 second.metric   =sde.genes,
                                 my.filename     ="Data/Reproducible-results/Data/two_metric_optimization.rds"){
  
  #Setting up weights, loop-indexer, and the list that will store the results
  weights <- seq(from = my.start, to=my.finish, by=step.size)
  integrated_gene_lists <- vector(mode = "list", length = 11)
  list_num <- seq(from = 1, to=length(weights), by=1)
  
  #Doing the grid search
  for (w in weights) {
    current_ranking <- two_metric_geneRank(ranking1 = first.metric,
                                           ranking2 = second.metric,
                                           a1=w, a2=1-w)
    current_ranking <- as.data.frame(current_ranking)
    integrated_gene_lists[[list_num[which(weights==w)]]] <- current_ranking
  }
  
  #Saving the output of the grid search to a .rds
  saveRDS(integrated_gene_lists, file = my.filename)
  return(integrated_gene_lists)
}


#Old 2 weight optimizer----
two_weight_optimizer_old <- function(my.start        =0,
                                 my.finish       =1,
                                 step.size       =0.1,
                                 first.metric    =mad.genes,
                                 second.metric   =sde.genes,
                                 my.filename     ="Data/Reproducible-results/Data/two_metric_optimization.RData"){
  
  #Setting up weights, loop-indexer, and the list that will store the results
  weights <- seq(from = my.start, to=my.finish, by=step.size)
  integrated_gene_lists <- vector(mode = "list", length = 11)
  counter <- 1
  
  #Doing the grid search
  for (w in weights) {
    current_ranking <- two_metric_geneRank(ranking1 = first.metric,
                                           ranking2 = second.metric,
                                           a1=w, a2=1-w)
    current_ranking <- as.data.frame(current_ranking)
    integrated_gene_lists[[counter]] <- current_ranking
    counter <- counter + 1
  }
  
  #Saving the output of the grid search to a .rds
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
                                   my.filename     ="Data/Exported-data/R-objects/three_metric_optimization.rds"){
  
  #Setting up weights, the starting value for the a3 metric, loop-indexer, and the list that will store the results
  weights <- seq(from = my.start, to=my.finish, by=step.size)
  a3_weights <- seq(from = my.start, to=my.finish, by=step.size)
  a3 <- a3.start
  df_index <- my.index
  integrated_gene_lists <- vector(mode = "list", length = 121)
  
  #Doing the grid search
  for (x in a3_weights){
    for (y in weights) {
      current_ranking <- three_metric_geneRank(ranking1 = first.metric, ranking2 = second.metric, ranking3 = third.metric,  a1=x, a2=1-(x+a3), a3= y)
      current_ranking <- as.data.frame(current_ranking)
      integrated_gene_lists[[df_index]] <- current_ranking
      df_index <- df_index + 1
    }
  }
  
  #Saving the output of the grid search to a .rds file
  saveRDS(integrated_gene_lists, file = my.filename)
  return(integrated_gene_lists)
}

