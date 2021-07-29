#c-index-graphic.R
#Simply reading in the file and plotting it.
#Loading needed packages----
library(ggplot2)
library(ggpubr)

#Setting the working directory----
setwd("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/")

#Reading in the file----
my_file <- read.csv("Completed_KM_Curves/Final-c-index-and-km-curve-document.csv")
my_file <- subset(my_file, select=c(1:4, 19:23))
my_file <- filter(my_file, scRNA.seq.dataset != "GSE81861_Cell_line_hct116")
my_file$pvalue_trans <- -log10(my_file$pvalue_finished_corrected_risk)
res_aov_1 <- aov(C.index_at_optimal_genes_corrected_risk~Method, data = my_file)
aov_sum_1 <- summary(res_aov_1)
pvalue_to_plot_1 <- round(aov_sum_1[[1]][["Pr(>F)"]][1], digits = 5)
tukey_aov_1 <- TukeyHSD(res_aov_1)


#genes only
my_file <- filter(my_file, Method=="MiRNA + SDES" | Method== "MiRNA + MAD" | Method== "scDD (genes only)" | Method== "DESeq2 (genes only)" | Method== "DESeq2 (genes only)" | Method=="edgeR (genes only)" | Method=="DEsingle (genes only)" | Method=="MAD" | Method=="MiRNA" | Method=="SDES")
#genes + tumor stage
my_file <- filter(my_file, Method=="MiRNA + SDES + Tumor stage" | Method== "MiRNA + MAD + Tumor stage" | Method== "scDD (genes + Tumor stage)" | Method== "DESeq2 (genes only + Tumor stage)" | Method=="edgeR (genes only + Tumor stage)" | Method=="DEsingle (genes + Tumor stage)" | Method=="MAD + Tumor stage" | Method=="MiRNA + Tumor stage" | Method=="SDES + Tumor stage")
#genes + tumor and N stage
my_file <- filter(my_file, Method=="MiRNA + SDES + Tumor stage + N stage" | Method== "MiRNA + MAD + Tumor stage + N stage" | Method== "scDD (genes + Tumor stage + N stage)" | Method== "DESeq2 (genes only + Tumor stage + N stage)" | Method=="edgeR (genes only + Tumor stage + N stage)" | Method=="DEsingle (genes + Tumor stage + N stage)" | Method=="MAD + Tumor stage + N stage" | Method=="MiRNA + Tumor stage + N stage" | Method=="SDES + Tumor stage + N stage")
#Making the names cleaner
my_file <- my_file %>% group_by(Method)
my_file[c(1,3,6), "Method"] <- "CC Singlecell MS"
my_file[c(2,4,5), "Method"] <- "CC Singlecell MM"
my_file[7:9,"Method"] <- "scDD"
my_file[10:12,"Method"] <- "DEsingle"
my_file[13:15, "Method"] <- "DESeq2"
my_file[16:18, "Method"] <- "edgeR"
my_file[19:21, "Method"] <- "MAD"
my_file[22:24, "Method"] <- "SDES"
my_file[25:27, "Method"] <- "MiRNA"


#Making nicer graphs----
#my_file <- filter(my_file, Bulk.dataset== "TCGA-COAD" | Bulk.dataset=="TCGA-READ")
#my_comparisons <- list(c("DEsingle, MiRNA + MAD"), c("DEsingle, MiRNA + SDES"), c("scDD, MiRNA + MAD"), c("scDD, MiRNA + SDES"))
#my_file$Method <- factor(my_file$Method, levels = c("MAD","SDES","MiRNA","MiRNA + SDES","MiRNA + MAD","DEsingle","scDD","DESeq2","edgeR"))
my_file$Method <- factor(my_file$Method, levels = c("MAD", "SDES", "MiRNA", "CC Singlecell MS", "CC Singlecell MM", "DEsingle", "scDD", "DESeq2", "edgeR"))
my_file$Bulk.dataset <- factor(my_file$Bulk.dataset, levels = c("TCGA-COAD", "TCGA-READ", "TCGA-COAD + TCGA-READ"))

#C-index graph----
cindex_bar3 <- ggplot(my_file, aes(x=Method, y = C.index_at_optimal_genes_corrected_risk, fill=Bulk.dataset))+
  geom_bar(stat= "identity",position = position_dodge())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("Gene Signature & All Clinical Data C-indices")+
  xlab("Method")+
  ylab("Concordance Index")+
  scale_fill_discrete("Bulk Dataset")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  geom_hline(yintercept = 0.70)



cindex_bar3 <- cindex_bar3 + coord_cartesian(ylim = c(0.5, 0.7))
cindex_bar3


#P-value graph----
pvalue_bar3 <- ggplot(my_file, aes(x=Method, y =pvalue_trans, fill=Bulk.dataset))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("Gene Signature & All Clinical Data P-values")+
  xlab("Method")+
  ylab("-log(P-value)")+
  scale_fill_discrete("Bulk Dataset")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  geom_hline(yintercept = -log10(0.05))


pvalue_bar3

#All graphs together----
big_graph <- ggarrange(cindex_bar, pvalue_bar, cindex_bar2, pvalue_bar2, cindex_bar3, pvalue_bar3, labels = c("A.", "B.", "C.", "D.", "E.", "F."), ncol = 1, nrow = 6, common.legend = TRUE)



#Set size of active genes----
set_bar <- ggplot(my_file, aes(x=Method, y = glmnet_active_genes_optimal_corrected_risk, fill=Bulk.dataset))+
  geom_bar(stat= "identity",position = position_dodge())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("Genes")+
  xlab("Method")+
  ylab("Set Size")+
  scale_fill_discrete("Bulk Dataset")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))


set_bar


big_set_graph <- ggarrange(set_bar, set_bar2, set_bar3, labels = c("A.", "B.", "C."), ncol = 3, nrow = 1, common.legend = TRUE)
big_set_graph


#Intersect plot----
#TODO


#Coefficients graph----



#Boxplot with dotplot----
# p<- ggplot(my_file, aes(x=Method, y=cindex_finished, fill=Method))+
#   geom_boxplot(data=my_file, aes(x=Method, y=cindex_finished))+
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#         legend.position = "none", panel.background = element_blank(),
#         axis.line = element_line(colour = "grey"))+
#   ggtitle("Gene Signatures Only")+
#   xlab("Method")+
#   ylab("Concordance Index")
# 
# p
# 
# 
# finished_plot <- p+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
#   stat_compare_means(method = "anova", label.x = 1, label.y = 0.56)+
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "MiRNA + MAD", vjust = -0.8)
#   #scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
#   #scale_x_discrete(guide = guide_axis(n.dodge = 3))
#   
#Power analysis----
summary_df <- my_file %>% dplyr::group_by(Method) %>% dplyr::summarise(m=median(C.index_at_1800_genes))
groupmedians <- summary_df$m

power_analysis <- power.anova.test(groups = length(groupmedians), 
                      between.var = var(groupmedians), within.var = 0.03, 
                      power=0.95,sig.level=0.05,n=NULL)
power_analysis






#Now doing simulating datasets to try and see if we get better performance for stats----
library(SPsimSeq)
set.seed(1)

cox_df_cc_counts <-subset(cox_df, select=(TSPAN6:AC007389.3))
cox_df_cc_counts <- t(cox_df_cc_counts)
cox_df_cc_counts <- as.data.frame(cox_df_cc_counts)

sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = cox_df_cc_counts,
                          group = 1, n.genes = 55317,
                          tot.samples = 501, 
                          pDE = 0.1, lfc.thrld = 0.5, 
                          result.format = "list", verbose = TRUE)



library(seqgendiff)
set.seed(1)




library(splatter)
params <- splatEstimate(cc_tumor_fpkm)
sim <- splatSimulate(params, nGenes = 55186)




#Plotting pvalues----


pvalue_plot <- ggplot(data = my_file, aes(x=Method, y=))

#cv_master_subset$Method <- factor(cv_master_subset$Method, levels = c("MAD", "SDES", "MiRNA", "MAD + SDES", "MiRNA + SDES", "MAD + SDES + MiRNA"))
my_file$MAD.Metric <- formatC(my_file$MAD.Metric, format="f", digits=2)
my_file$SDES.Metric <- formatC(my_file$SDES.Metric, format="f", digits=2)
my_file$MiRNA.Metric <- formatC(my_file$MiRNA.Metric, format="f", digits=2)
my_file$MiRNA...SDES.Metric <- formatC(my_file$MiRNA...SDES.Metric, format = "f", digits = 2)


p<-ggplot(data=my_file, aes(x=c(MAD.Metric,SDES.Metric,MiRNA.Metric, MiRNA...SDES.Metric), y=c(MAD.Metric,SDES.Metric,MiRNA.Metric, MiRNA...SDES.Metric), fill=c(MAD.Metric,SDES.Metric,MiRNA.Metric, MiRNA...SDES.Metric))) +
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(title = element_text(size=35), axis.text = element_text(size = 20), axis.text.x = element_text(size = 20, face = "bold"))+
  labs(title="10-fold CV Performance Across Methods", 
       x="Method", y = "Mean 10-fold CV C-index")+
  scale_y_discrete(expand = expansion(mult = c(0, .1)))+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=MAD.Metric), vjust=0.0, color="black", size=15)

p+ theme(legend.position="none")

