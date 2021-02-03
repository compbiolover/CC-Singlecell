weights <- seq(from = 0, to=1, by=0.1)
df_index <- 1
integrated_gene_lists <- list()
a3_weights <- seq(from = 0, to=1, by=0.1)


for (x in a3_weights){
  print(x)
  for (y in weights) {
    print(y)
    current_ranking <- geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3= y)
    current_ranking <- as.data.frame(current_ranking)
    integrated_gene_lists[[df_index]] <- current_ranking
    df_index <- df_index + 1
  }
}
