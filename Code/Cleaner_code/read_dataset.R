#read_dataset.R
#Purpose: TCGA-READ analysis
#Loading needed packages----
library(ggplot2)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)


#Making the READ data set----
read_query <- GDCquery(project       = "TCGA-READ",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the files
GDCdownload(query           = read_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-READ-Dataset")


#Making the SummarizedExperiment object
read_data_se <- GDCprepare(read_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-read-Dataset/")
read_data_df <- as.data.frame(colData(read_data_se))
read_data_df$vital_status <- factor(read_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
read_data_df$vital_status <- as.numeric(as.character(read_data_df$vital_status))

#Bulk data frame for Read merged data frame----
load("Data/TCGA-READ/read_df.RData", verbose = TRUE)

#MiRNA metric----
mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                cancer.type2 = "colon cancer",
                                ts.org = "Human",
                                ts.version = "7.2",
                                max.miR.targets = 10,
                                cancer.up = TRUE,
                                mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                max.mirnas = 1558,
                                write.heatmap.data = FALSE, 
                                print.ts.targets = TRUE,
                                save.mirna.raw.targets = FALSE,
                                save.mirna.genes = TRUE,
                                mirna.gene.filename = "Data/TCGA-READ/mirna_genes_1558_10_targets.csv",
                                mirna.gene.rfile = "Data/TCGA-READ/mirna_genes_1558_10_targets.RData")


#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "Data/TCGA-READ/MiRNA/three_weight_optimized_output_1558_mirna_10_top.csv",
                                      my.title = "CC Singlecell MMS READ 1558 MiRNA 100 Targets Only",
                                      tumor.data = FALSE, n.data = FALSE, cox.df = cox_df, show.pval = TRUE, show.pval.method = FALSE)
patient_risk


load("Data/TCGA-READ/three_weight_optimized_40_mirna_100_colorect.RData", verbose = TRUE)
mad_sdes_mirna_optimized <- integrated_gene_lists

cox_models <- list()
my_cindicies <- c()
counter <- 1

for (x in mad_sdes_mirna_optimized[2]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/top_performing_read.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
}

top_cindex <-max(my_cindicies)
top_index <- which(my_cindicies==top_cindex)
top_index

cox_models$`1`$CV




#For plotting the coefficients----
coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
coef_df_sub <- filter(coef_df, abs(coefs)>0.1)

#coef_df_sub$p.value.trans <- -log10(coef_df_sub$p.value)
colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("CC Singlecell MMS READ Cox Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot


#Just the 0.1 weighted coefficients----
top_cox <- coxph(data = cox_df, formula = Surv(days.to.last.follow.up, vital.status)~TEDDM1+SLC22A13+OR2H1+NUTM1+LCE1C+FSHB+FSD2+CMC4+AIPL1)
summary(top_cox)



