#glio_dataset.R
#Load needed packages----
library(ggplot2)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)

#Glioblastoma patients----
glio <- read.csv("Data/Single-cell-data/Other-cancers/GSE57872_GBM_data_matrix.txt", sep='\t')
rownames(glio) <- glio$X
glio_dg <- dplyr::select(glio, contains("MGH2") | contains("MGH26Tumor") | contains("MGH28Tumor") | contains("MGH31Tumor"))
glio_dg <- abs(glio_dg)


rownames(glio) <- glio$X
glio <- subset(glio, select=c(MGH264_A01:MGH31Tumor))
glio <- abs(glio)



#Glioblastoma TCGA bulk data----
glio_query <- GDCquery(project        = "TCGA-GBM",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = glio_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/TCGA-GBM/TCGA-Glio-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
Glio_data_se <- GDCprepare(glio_query, summarizedExperiment = TRUE, directory = "Data/TCGA-GBM/TCGA-Glio-Dataset/")
Glio_data_df <- as.data.frame(colData(Glio_data_se))
Glio_data_df$vital_status <- factor(Glio_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
Glio_data_df$vital_status <- as.numeric(as.character(Glio_data_df$vital_status))

#Bulk dataframe for glio merged dataframe----
bulk_rna_df <- Glio_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- Glio_data_se@colData@rownames
rownames(bulk_rna_df) <- Glio_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
Glio_data_df_unique <- subset(Glio_data_df, select = unique(colnames(Glio_data_df)))
merged_df <- merge(bulk_rna_df_unique, Glio_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
calculated_days[calculated_days==0]=1
merged_df <- merged_df[complete.cases(merged_df[, "days_to_last_follow_up"]), ]
cox_time <- merged_df$days_to_last_follow_up
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$tumor_grade
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df <- cox_df %>% filter(!is.na(vital.status))
cox_df$days.to.last.follow.up <- as.numeric(cox_df$days.to.last.follow.up)
cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up==0,1,cox_df$days.to.last.follow.up)
#save(cox_df, file = "Data/TCGA-GBM/gbm_df.RData")



#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human",
                                ts.version = "7.2",
                                max.mirnas = 400,
                                max.miR.targets = 10,
                                cancer.up = TRUE,
                                write.heatmap.data = FALSE,
                                print.ts.targets = TRUE,
                                save.mirna.raw.targets = FALSE,
                                save.mirna.genes = TRUE,
                                cancer.type1 = "brain cancer",
                                mirna.gene.rfile = "Data/TCGA-GBM/mirna_genes_400_10_targets.RData",
                                mirna.gene.filename = "Data/TCGA-GBM/mirna_genes_400_10_targets.csv",
                                mirna.remove = "hsa-miR-129-1-3p")



#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "Data/TCGA-GBM/MiRNA/three_weight_optimized_output_400_mirna_10_top.csv",
                                      my.title = "CC Singlecell MMS GBM",
                                      tumor.data = FALSE,
                                      n.data = FALSE,
                                      cox.df = cox_df,
                                      show.pval = TRUE,
                                      show.pval.method = FALSE)
patient_risk








#For plotting the coefficients----
coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
coef_df_sub <- filter(coef_df, abs(coefs)>0.1)


colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("CC Singlecell MMS GBM Cox Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot










