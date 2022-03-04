#committee_meeting1.R
#Code for figures in committee meeting report 1

#Loading needed packages----
library(ggpubr)
library(grid)
library(gridExtra)
library(pheatmap)
library(png)
library(reshape2)
library(tidyverse)


#random matrix of unfiltered data----
ran_nums <- runif(100, min = 0, max = 500)
ran_mat <- matrix(ran_nums, nrow = 10, ncol = 10)
gene_heatmap <- pheatmap(ran_mat, cluster_cols = FALSE, cluster_rows = FALSE, labels_row = NULL, labels_col = NULL, angle_col = 0, legend = FALSE)


#random filtered heatmap----
ran_nums <- runif(60, min = 0, max = 250)
ran_mat <- matrix(ran_nums, nrow = 6, ncol = 10)
gene_heatmap2 <- pheatmap(ran_mat, cluster_rows = FALSE, cluster_cols = FALSE, labels_row = NULL, labels_col = NULL, legend = FALSE)


#Psuedotime----
cds_output <- cell_dataset_builder(vim.genes = c("VIM"), cell.data = cc_tumor_fpkm$denoised_sc_dataframe, cell.meta = cc_tumor_fpkm$cds_gene_names, graph.root = 2.5, show.traj = TRUE, point.size = 2.5)
test <- cds_output$`PT Graph`
test <- test + theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none",
  axis.title.x = element_text(family = "sans", size = 40),
  axis.title.y = element_text(family = "sans", size = 40))


#Mirna interaction map
ran_nums <- runif(100, min = 0, max = 8)
ran_mat <- matrix(ran_nums, nrow = 10, ncol = 10)
mirna_interactions <- pheatmap(ran_mat, cluster_rows = FALSE, cluster_cols = FALSE, labels_row = NULL, labels_col = NULL, legend = FALSE)





#Bar plot of different genes in our mirna list----
mirna_list <- read.csv("~/Documents/PhD Program/Hong Lab/Committee-meetings/Meeting_1/Figures/bigger_mirnas_for_heatmap.csv")
colnames(mirna_list) <- c("gene", "score")
mirna_list <- mirna_list[1:10,]


mirna_barplot <-ggplot(mirna_list, mapping = aes(x=gene, y=score))+
  geom_col(aes(fill=gene))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size = 14), panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("MiRNA Gene Scores")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

lm <-rbind(c(1,2))
big_mirna_graph <- grid.arrange(grobs = list(mirna_heatmap3[[4]], mirna_barplot),
                                layout_matrix = lm, top = textGrob("MiRNA Metric",gp=gpar(fontsize=20,font=2)), left = "Genes")

big_mirna_graph



#Merging the KM curves and set sizes graphs together into a single figure----
big_graph <- ggarrange(ggarrange(coad_risk$plot, read_risk$plot, coad_coef_plot, read_coef_plot, nrow = 2, ncol = 2, labels = c("A.", "B.", "C.", "D.")), set_bar, nrow = 2, heights = c(5:1))
big_graph

