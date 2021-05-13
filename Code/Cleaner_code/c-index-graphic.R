#c-index-graphic.R
#Simply reading in the file and plotting it.
#Loading needed packages----
library(ggplot2)

#Setting the working directory----
setwd("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/")

#Reading in the file----
my_file <- read.csv("Documentation/C-index-across-methods-table.csv")
my_file <- my_file[,1:15]
my_file <- my_file[c(1:4,6:9),]
my_file <- my_file[,c(1,4:15)]
my_file <- select(my_file, !contains(".SE"))
my_file <- select(my_file, !contains(".Genes"))


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

