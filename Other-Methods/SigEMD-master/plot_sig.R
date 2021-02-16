library(ggplot2)
plot_emd_density_sig <- function(emdobj, gene_name) {
  
  if (gene_name %in% emdobj$Hur_gene){
    gene.data <- emdobj$Hur$data[gene_name,]
    outcomes <- emdobj$Hur$outcomes
  } else {
    idx<- which(emdobj$nonHur_gene==gene_name)
    gene.data <- emdobj$nonHur$data[[idx]]
    outcomes <- emdobj$nonHur$outcomes[[idx]]
  }
  
  
  classes <- unique(outcomes)
  
  emd_score <- emdobj$emdall[gene_name, "emd"]
  q_value<- emdobj$emdall[gene_name, "padjust"]
  # to appease CRAN
  group <- NULL
  exp <- NULL
  
  df<-data.frame(row.names=colnames(gene.data), group=outcomes, exp=gene.data)
  
  title <- paste(gene_name, "\n", "(emd score = ",
                 round(emd_score, 2), ")","\n","(adjust pvalue = ",
                 round(q_value, 3), ")",sep="")
  
  ggplot(df, aes(exp, fill=group)) + geom_density(alpha=0.5) +
    xlab("data")  + ggtitle(title) +
    theme(axis.text=element_text(size=24),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size=24))
} 
