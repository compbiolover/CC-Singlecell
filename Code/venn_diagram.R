#venn_diagram.R
library(ggvenn)

all_data <- list(`CC Singlecell MM` = head(rownames(mirna_mad_optimized$`6`), n=500),
                 `Active Cox Genes` = cox_6_read_cc_mm$active_genes)

read_venn_cc_singlecell_mm <- ggvenn(all_data, c("CC Singlecell MM", "Active Cox Genes"), show_percentage = FALSE, set_name_size = 3)+
  ggtitle("TCGA-READ")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

big_venn <- ggarrange(coad_venn_cc_singlecell_mm, read_venn_cc_singlecell_mm, coad_venn_cc_singlecell_ms, read_venn_cc_singlecell_ms, nrow = 2, ncol = 2, labels = c("A.", "B.", "C.", "D."))
