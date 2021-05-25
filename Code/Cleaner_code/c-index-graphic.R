#c-index-graphic.R
#Simply reading in the file and plotting it.
#Loading needed packages----
library(ggplot2)

#Setting the working directory----
setwd("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/")

#Reading in the file----
my_file <- read.csv("Documentation/c-index-across-datasets-formatted-correctly-for-graphs.csv")
res_aov_1 <- aov(C.index~Method, data = my_file)
aov_sum_1 <- summary(res_aov_1)
pvalue_to_plot_1 <- round(aov_sum_1[[1]][["Pr(>F)"]][1], digits = 5)
tukey_aov_1 <- TukeyHSD(res_aov_1)


other_file <- read.csv("Documentation/other-methods-reordered.csv")
other_file <- filter(other_file, Bulk.dataset!="TCGA-COAD + TCGA-READ")

res_aov <- aov(C.index~Method, data = my_file)
aov_sum <- summary(res_aov)
pvalue_to_plot <- round(aov_sum[[1]][["Pr(>F)"]][1], digits = 5)

tukey_aov <- TukeyHSD(res_aov)


#my_file <- filter(my_file, Method=="scDD (genes only)" | Method== "DEsingle (genes only)" | Method== "MiRNA + MAD" | Method== "MiRNA + SDES")
my_file <- filter(my_file, Method=="MiRNA + SDES Metric + Tumor Stage + N Stage" | Method== "scDD (genes + tumor stage + N stage)" | Method== "DEsingle (genes + tumor stage + N stage)" | Method== " MiRNA + MAD + Tumor stage + N Stage")
my_file[1:18, "Method"] <- "CC Singlecell"
my_file[19:27,"Method"] <- "scDD"
my_file[28:36,"Method"] <- "DEsingle"
my_file[37:42, "Method"] <- "CC Singlecell"
my_file <- filter(my_file, Bulk.dataset== "TCGA-COAD" | Bulk.dataset=="TCGA-READ")
my_comparisons <- list(c("DEsingle, MiRNA + MAD"), c("DEsingle, MiRNA + SDES"), c("scDD, MiRNA + MAD"), c("scDD, MiRNA + SDES"))

p<- ggplot(my_file, aes(x=Method, y=C.index, fill=Method))+
  geom_boxplot(data=my_file, aes(x=Method, y=C.index))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  ggtitle("Concordance Index Performance Across Methods for All Datasets & Combinations")+
  xlab("Method")+
  ylab("Concordance Index")

p


p+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  stat_compare_means(method = "anova", label.x = 1, label.y = 0.60)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "CC Singlecell")
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
  #scale_x_discrete(guide = guide_axis(n.dodge = 3))
  





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

