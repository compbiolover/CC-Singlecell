#caod_and_read_dataset.R
#Loading the packages----
library(ggplot2)
library(TCGAbiolinks)
library(tidyverse)

#Loading the data set----
combined_query <- GDCquery(project       = c("TCGA-READ", "TCGA-COAD"),
                           data.category = "Transcriptome Profiling",
                           data.type     = "Gene Expression Quantification",
                           workflow.type = "HTSeq - Counts")

#Downloading
GDCdownload(query           = combined_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/TCGA-COAD_and_TCGA-READ/TCGA-COAD-and-TCGA-READ-Dataset")

combined_data_se <- GDCprepare(combined_query, summarizedExperiment = TRUE, directory = "Data/TCGA-COAD_and_TCGA-READ/TCGA-COAD-and-TCGA-READ-Dataset/")
combined_data_df <- as.data.frame(colData(combined_data_se))
combined_data_df$vital_status <- factor(combined_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
combined_data_df$vital_status <- as.numeric(as.character(combined_data_df$vital_status))



#Bulk data frame for combined merged data frame----
# combined_bulk_rna_df <- combined_data_se@assays@data@listData[["HTSeq - Counts"]]
# colnames(combined_bulk_rna_df) <- combined_data_se@colData@rownames
# rownames(combined_bulk_rna_df) <- combined_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
# combined_bulk_rna_df <- t(combined_bulk_rna_df)
# combined_bulk_rna_df <- as.data.frame(combined_bulk_rna_df)
# bulk_rownames <- rownames(combined_bulk_rna_df)
# combined_bulk_rna_df$barcode <- bulk_rownames
# 
# combined_bulk_rna_df_unique <- subset(combined_bulk_rna_df, select = unique(colnames(combined_bulk_rna_df)))
# combined_data_df_unique <- subset(combined_data_df, select = unique(colnames(combined_data_df)))
# merged_df <- merge(combined_bulk_rna_df_unique, combined_data_df_unique, by = 'barcode')
# rownames(merged_df) <- merged_df$barcode
# merged_df <- merged_df[,2:length(colnames(merged_df))]
# 
# merged_df$days_to_death <- ifelse(is.na(merged_df$days_to_death),0, merged_df$days_to_death)
# merged_df$days_to_last_follow_up <- ifelse(is.na(merged_df$days_to_last_follow_up),0, merged_df$days_to_last_follow_up)
# 
# calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
# calculated_days <- abs(calculated_days)
# calculated_days[calculated_days==0]=1
# cox_time <- calculated_days
# cox_event <- merged_df$vital_status
# cox_tumor <- merged_df$ajcc_pathologic_stage
# cox_tumor_n <- merged_df$ajcc_pathologic_n
# cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
# cox_df$days.to.last.follow.up <- cox_time
# cox_df$vital.status <- cox_event
# cox_df$tumor.stage <- cox_tumor
# cox_df$ajcc.n <- cox_tumor_n
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage IV", replacement = 4)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage III", replacement = 3)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage II", replacement = 2)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage I", replacement = 1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
# cox_df <- filter(cox_df, !tumor.stage=="not reported")
# cox_df <- cox_df[complete.cases(cox_df[, "ajcc.n"]), ]
# cox_df <- cox_df[complete.cases(cox_df[, "vital.status"]), ]
# save(cox_df, file = "Data/TCGA-COAD_and_TCGA-READ/coad_and_read_df.RData")
load("Data/TCGA-COAD_and_TCGA-READ/coad_and_read_df.RData")

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
                                mirna.gene.filename = "Data/TCGA-COAD_and_TCGA-READ/mirna_genes_1558_10_targets.csv",
                                mirna.gene.rfile = "Data/TCGA-COAD_and_TCGA-READ/mirna_genes_1558_10_targets.RData")
#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "Data/TCGA-COAD_and_TCGA-READ/three_weight_optimized_output_1558_mirna_10_top.csv",
                                      my.title = "CC Singlecell MMS COAD + READ 1558 MiRNA 100 Targets",
                                      tumor.data = FALSE, n.data = FALSE, cox.df = cox_df, show.pval = TRUE, show.pval.method = FALSE)
patient_risk


load("", verbose = TRUE)
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
                                  my.filename = "Data/TCGA-COAD_and_TCGA-READ/top_performing_genes.csv") 
  
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


colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("CC Singlecell MMS COAD + READ Cox Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot


#Just the 0.1 weighted coefficients----
top_cox <- coxph(data = cox_df, formula = Surv(days.to.last.follow.up, vital.status)~DEFA7P+AIFM1P1+AC079804.1)
summary(top_cox)


