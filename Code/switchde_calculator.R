#Name: switchde_calculator.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently build Inference of 
#switch-like differential expression along
#single-cell RNA seq trajectories (switchde)
#outputs

#The switchde (SDE) metric ----
switchde_calculator <- function(denoised.sc=sc, 
                                pseudo.time=pt, zero.inflated=FALSE){
  
  #Loading/installing needed package----
  require(switchde)
  require(tidyverse)
  
  #Actual metric----
  sdes <- switchde(denoised.sc, as.numeric(pseudo.time$Pseudotime), verbose = TRUE, zero_inflated = zero.inflated)
  sde.filtered <- filter(sdes, qval < 0.05)
  index <- order(abs(sde.filtered$k), decreasing = T)
  vim.sdes.rank <- sde.filtered[index,]
  vim.sdes.ranking <- vim.sdes.rank$k
  names(vim.sdes.ranking) <- vim.sdes.rank$gene
  vim.sdes.ranking<-abs(vim.sdes.ranking)/sum(abs(vim.sdes.ranking))
  
  #Return object----
  return(vim.sdes.ranking)
}
