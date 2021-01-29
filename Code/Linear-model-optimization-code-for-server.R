#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Containsn just code needed to optimize linear model on
#lab server

#Loading the needed packages----

#Loading needed functions----
# Ranking genes by three measurements
geneRank <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, a1 = 1,
                     a2 = 0, a3 = 0){
  gns = c(names(ranking1),names(ranking2),names(ranking3))
  gns = unique(gns)
  res = rep(0, length(gns))
  names(res) = gns
  for(i in names(res)) {
    res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2+getRank(ranking3, i)*a3
  }
  #res=res/sum(res)
  #res = res[order(res, decreasing = T, method = "radix")]
  res
}

getRank <- function(ranking = NULL, gn = NULL){
  if (gn %in% names(ranking)) {
    return(ranking[gn])
  }
  else return(0.0)
}