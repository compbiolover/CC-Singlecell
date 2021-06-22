#Name: main.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: The main script for my project.

#Loading needed packages----
library(ggplot2)
library(grid)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)

#Loading single-cell data----
cc_tumor_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_tumor_all_cells_FPKM.csv")

k562_cells <- read.csv("Data/Single-cell-data/Other-cancers/GSE65525_RAW/GSM1599500_K562_cells.csv")

glio <- read.csv("Data/Single-cell-data/Other-cancers/GSE57872_GBM_data_matrix.txt", sep='\t')
rownames(glio) <- glio$X
glio_dg <- dplyr::select(glio, contains("MGH2") | contains("MGH26Tumor") | contains("MGH28Tumor") | contains("MGH31Tumor"))
glio_dg <- abs(glio_dg)

cc_cell_line_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_Cell_Line_FPKM.csv")

lc_tumor_tpm <- read.csv("Data/Single-cell-data/Other-cancers/GSE69405_PROCESSED_GENE_TPM_ALL.txt", sep = '\t')
lc_tumor_tpm <- lc_tumor_tpm[!base::duplicated(lc_tumor_tpm$gene_name),]
rownames(lc_tumor_tpm) <- lc_tumor_tpm$gene_name

#Loading in gene sets from other methods to compare to our method----
scdd_genes <- read.csv("Data/Exported-data/Csv-files/scdd_res.csv")
rownames(scdd_genes) <- scdd_genes$X
scdd_genes <- scdd_genes[,c("gene", "nonzero.pvalue", "nonzero.pvalue.adj")]


#A method to extract the active genes from the cox models----
active_gene_extractor <- function(cox.num=1){
  Active.Index <- which(as.logical(cox_models[[cox.num]][2]$Coefficients) != 0)
  Active.Coefficients  <- cox_models[[cox.num]][2]$Coefficients[Active.Index]
  active_genes <-rownames(cox_models[[cox.num]][2]$Coefficients)[Active.Index]
  return_df <- data.frame(active_genes=active_genes, betas=Active.Coefficients)
  return_df<- return_df[base::order(return_df$betas, decreasing = TRUE),]
  rownames(return_df) <- 1:length(return_df$betas)
  return(return_df)
}
#DESingle CC patients----
des_genes <- des_results
des_genes <- head(rownames(des_genes), n=1800)



#Log + 1 transforming the cc_tumor sc data to see if it impacts results----
cc_tumor_fpkm_logged <- log1p(cc_tumor_fpkm[,2:length(colnames(cc_tumor_fpkm))])
cc_tumor_fpkm_logged$X <- cc_tumor_fpkm$X

cc_cell_line_fpkm_logged <- log1p(cc_cell_line_fpkm[,2:length(colnames(cc_cell_line_fpkm))])
cc_cell_line_fpkm_logged$X <- cc_cell_line_fpkm$X


#Loading additional bulk colon cancer datasets----
sw620_sw480 <- read.csv("Data/Datasets-for-comparison/GSE112568_fpkm_table_CRC_C_Williams_Addendum.txt", sep = '\t')

#Loading the types of cancer we can pick from for mirna metric----
high_choices <- read.csv("Data/Data-from-Cleaner-code/high-dbDEMC-cancer-types.csv")
low_choices <- read.csv("Data/Data-from-Cleaner-code/low-dbDEMC-cancer-types.csv")

#Pre-processing for the different sc data files
#CC tumor patients----
rownames(cc_tumor_fpkm) <- cc_tumor_fpkm$X
current_colname_split <- strsplit(colnames(t(cc_tumor_fpkm)), "_")
finished_gene_list <- c()
current_list <- current_colname_split
for (x in seq(1:length(current_list))){
  finished_gene_list <- c(finished_gene_list, current_list[[x]][2])
}

cc_tumor_fpkm <- t(cc_tumor_fpkm)
colnames(cc_tumor_fpkm) <- finished_gene_list
cc_colnames <- unique(colnames(cc_tumor_fpkm))
cc_tumor_fpkm <- t(cc_tumor_fpkm)
cc_tumor_fpkm <- subset(cc_tumor_fpkm, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
cc_tumor_fpkm <- t(cc_tumor_fpkm)
cc_tumor_fpkm <- subset(cc_tumor_fpkm, select=cc_colnames)
cc_tumor_fpkm <- t(cc_tumor_fpkm)
cc_tumor_fpkm <- apply(cc_tumor_fpkm, c(1,2), as.numeric)
cc_tumor_fpkm <- as.matrix(cc_tumor_fpkm)


#CC tumor patients logged----
rownames(cc_tumor_fpkm_logged) <- cc_tumor_fpkm_logged$X
current_colname_split <- strsplit(colnames(t(cc_tumor_fpkm_logged)), "_")
finished_gene_list <- c()
current_list <- current_colname_split
for (x in seq(1:length(current_list))){
  finished_gene_list <- c(finished_gene_list, current_list[[x]][2])
}

cc_tumor_fpkm_logged <- t(cc_tumor_fpkm_logged)
colnames(cc_tumor_fpkm_logged) <- finished_gene_list
cc_colnames <- unique(colnames(cc_tumor_fpkm_logged))
cc_tumor_fpkm_logged <- t(cc_tumor_fpkm_logged)
cc_tumor_fpkm_logged <- subset(cc_tumor_fpkm_logged, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
cc_tumor_fpkm_logged <- t(cc_tumor_fpkm_logged)
cc_tumor_fpkm_logged <- subset(cc_tumor_fpkm_logged, select=cc_colnames)
cc_tumor_fpkm_logged <- t(cc_tumor_fpkm_logged)
cc_tumor_fpkm_logged <- apply(cc_tumor_fpkm_logged, c(1,2), as.numeric)
cc_tumor_fpkm_logged <- as.matrix(cc_tumor_fpkm_logged)

#Lung cancer patients----
lc_tumor_tpm <- select(lc_tumor_tpm, contains("LC.PT.45" ) | contains("LC.MBT.15"))
lc_tumor_tpm <- select(lc_tumor_tpm, contains("LC.PT.45_SC") | contains( "LC.MBT.15_SC"))
lc_names <- c("VIM", "VIMP", "VIMP1")

#K562 leukemia cell line----
rownames(k562_cells) <- k562_cells$X
k562_cells <- subset(k562_cells, select=c(X0:X0.238))
colnames(k562_cells) <- seq(from=1,to=239,by=1)

#Glioblastoma patients----
rownames(glio) <- glio$X
glio <- subset(glio, select=c(MGH264_A01:MGH31Tumor))
glio <- abs(glio)

#LC tumor patient dataset analysis----
lc_tumor_tpm <- magic_denoiser(sc.data = lc_tumor_tpm, magic.seed = 123)
cds_output <- cell_dataset_builder(vim.genes = lc_names, cell.data = lc_tumor_tpm$denoised_sc_dataframe, cell.meta = lc_tumor_tpm$cds_gene_names)

mad.genes <- mad_calculator(lc_tumor_tpm$denoised_sc_dataframe)
save(mad.genes, file = "Data/Data-from-Cleaner-code/lc_tumor_tpm_mad.RData")

sde.genes <- switchde_calculator(lc_tumor_tpm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "Data/Data-from-Cleaner-code/lc_tumor_tpm_sde.RData")

mirna.genes <- mirna_calculator(cancer.type1 = "lung cancer", ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, mirna.filename = "Data/Data-from-Cleaner-code/TargetScan_lc_tumor.RData", mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"))
save(mirna.genes, file = "Data/Data-from-Cleaner-code/lc_tumor_tpm_mirna.RData")

mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "Data/Data-from-Cleaner-code/sde_mirna_optimized_lc.RData")


#Lung cancer TCGA bulk data----
lung_query <- GDCquery(project       = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = lung_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-lungluad-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
lung_data_se <- GDCprepare(lung_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-lungluad-Dataset/")
lung_data_df <- as.data.frame(colData(lung_data_se))
lung_data_df$vital_status <- factor(lung_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
lung_data_df$vital_status <- as.numeric(as.character(lung_data_df$vital_status))

#Bulk dataframe for lung merged dataframe----
bulk_rna_df <- lung_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- lung_data_se@colData@rownames
rownames(bulk_rna_df) <- lung_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
lung_data_df_unique <- subset(lung_data_df, select = unique(colnames(lung_data_df)))
merged_df <- merge(bulk_rna_df_unique, lung_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
calculated_days[calculated_days==0]=1
merged_df <- merged_df[complete.cases(merged_df[, "days_to_last_follow_up"]), ]
merged_df$days_to_last_follow_up <- ifelse(merged_df$days_to_last_follow_up==0,1,merged_df$days_to_last_follow_up)
cox_time <- merged_df$days_to_last_follow_up
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$tumor_stage
cox_tumor_n <- merged_df$ajcc_pathologic_n
cox_df <- subset(merged_df, select=c(TSPAN6:V56404))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.n <- cox_tumor_n
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df <- filter(cox_df, !tumor.stage=="not reported")
cox_df <- filter(cox_df, !ajcc.n=="NX")
cox_df <- cox_df[complete.cases(cox_df[, "ajcc.n"]), ]


#Cox model----
cox_models <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

save(cox_models, file = "Data/Data-from-Cleaner-code/cox_models.RData")

#Cox for individual metrics---
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
save(cox_mad, file = "Data/Data-from-Cleaner-code/cox_mad.RData")
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
save(cox_sdes, file = "Data/Data-from-Cleaner-code/cox_sde.RData")
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
save(cox_mirna, file = "Data/Data-from-Cleaner-code/cox_mirna.RData")

#KM curves----
hr_calcs <- list()
counter <- 1
for (x in cox_models[1:11]) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE)
  hr_calcs[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#HR calcs for individual metrics---
hr_mad <- hr_calculator(model.coefs = cox_mad$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE)
hr_sdes <- hr_calculator(model.coefs = cox_sdes$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE)
hr_mirna <- hr_calculator(model.coefs = cox_mirna$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE)

#KM p-values----
hr_pvalues <- list()
counter <- 1
for (x in hr_calcs[7:11]) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#P-values for HR calcs for individual metrics---
pvalue_mad <- km_pvalue_calculator(surv.time = hr_mad$DF$days.to.last.follow.up, surv.status = hr_mad$DF$vital.status, surv.predictors = hr_mad$DF$risk, surv.df = hr_mad$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_sdes <- km_pvalue_calculator(surv.time = hr_sdes$DF$days.to.last.follow.up, surv.status = hr_sdes$DF$vital.status, surv.predictors = hr_sdes$DF$risk, surv.df = hr_sdes$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_mirna <- km_pvalue_calculator(surv.time = hr_mirna$DF$days.to.last.follow.up, surv.status = hr_mirna$DF$vital.status, surv.predictors = hr_mirna$DF$risk, surv.df = hr_mirna$DF, num.sig.figs = 4, scientific.p.value = TRUE)


#KM plots----
hr_plots <- list()
counter <- 1
for (x in hr_calcs) {
  current_calc <- x
  current_pvalue <- hr_pvalues[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "SDES + MiRNA + Tumor Stage + N Stage Lung Cancer Patients Data")
  hr_plots[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}


#KM Plots for HR calcs for individual metrics---
km_mad <- km_plotter(km.fit = hr_mad$KM, data.source = hr_mad$DF, p.value = pvalue_mad, plot.title = "MAD Metric + Tumor Stage + N Stage Alone Lung Cancer Patients")
km_sdes <- km_plotter(km.fit = hr_sdes$KM, data.source = hr_sdes$DF, p.value = pvalue_sdes, plot.title = "SDES Metric + Tumor Stage + N Stage Alone Lung Cancer Patients")
km_mirna <- km_plotter(km.fit = hr_mirna$KM, data.source = hr_mirna$DF, p.value = pvalue_mirna, plot.title = "MiRNA Metric + Tumor Stage + N Stage Alone Lung Cancer Patients")


#CC tumor patients dataset analysis
#Now denoising the sc-data----
cc_tumor_fpkm <- magic_denoiser(sc.data = cc_tumor_fpkm, magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP", "CDH1", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "GRHL2", "CCND1", "CDK4", "CDK6"), cell.data = cc_tumor_fpkm$denoised_sc_dataframe, cell.meta = cc_tumor_fpkm$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
save(mad.genes, file = "Data/Data-from-Cleaner-code/cc_tumor_fpkm_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(cc_tumor_fpkm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "Data/Data-from-Cleaner-code/cc_tumor_fpkm_sde_ccnd1.RData")

#Mirna metric----
mirna_sizes1 <- seq(5, 50, by=10)
mirna_sizes2 <- seq(50, 200, by=50)
mirna_sizes3 <- seq(200, 1555, by=200)
mirna_sizes4 <- seq(1400, 1550, by=100)
mirna_sizes <- c(mirna_sizes1, mirna_sizes2, mirna_sizes3, mirna_sizes4)
mirna_sizes <- unique(mirna_sizes)
mirna_master_list <- list()
counter <- 1
for(y in mirna_sizes){
  print(y)
  current_mirna <- mirna_calculator(cancer.type1 = "colon cancer", cancer.type2 = "colorectal cancer", ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, mirna.filename = NULL, mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"), max.mirnas = y)
  mirna_master_list[[as.character(counter)]] <- current_mirna
  counter <- counter + 1
}



#For lung cancer----
mirna_master_list <- list()
counter <- 1
for(y in mirna_sizes){
  print(y)
  current_mirna <- mirna_calculator(cancer.type1 = "lung cancer", ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, mirna.filename = NULL, mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"), max.mirnas = y, print.ts.targets = TRUE)
  mirna_master_list[[as.character(counter)]] <- current_mirna
  counter <- counter + 1
}





mirna.genes <- mirna_calculator(cancer.type1 = "colon cancer", cancer.type2 = "colorectal cancer", ts.org = "Human", ts.version = "7.2", max.miR.targets = y, cancer.up = TRUE, mirna.filename = "Data/Data-from-Cleaner-code/TargetScan_cc_tumor_patients_all_targets.RData", mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"))
save(mirna_master_list, file = "Data/Data-from-Cleaner-code/cc_tumor_fpkm_first_200_mirna_at_10_targets_data.RData")

#Optimizing the MiRNA + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "Data/Data-from-Cleaner-code/sde_mirna_optimized.RData")

#Optimizing the MAD + MiRNA metric----
mirna_mad_optimized <- two_weight_optimizer(first.metric = mirna.genes, second.metric = mad.genes, my.filename = "~/Desktop/mirna_mad_optimized.RData")

#Loading the bulk READ dataset from TCGA----
read_query <- GDCquery(project       = "TCGA-READ",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

combined_query <- GDCquery(project       = c("TCGA-READ", "TCGA-COAD"),
                           data.category = "Transcriptome Profiling",
                           data.type     = "Gene Expression Quantification",
                           workflow.type = "HTSeq - Counts")

#Downloading the data
GDCdownload(query           = read_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-read-Dataset")

#Combined download
GDCdownload(query           = combined_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-combined-Dataset")

#READ summarizedExperiment object
read_data_se <- GDCprepare(read_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-read-Dataset/")
read_data_df <- as.data.frame(colData(read_data_se))
read_data_df$vital_status <- factor(read_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
read_data_df$vital_status <- as.numeric(as.character(read_data_df$vital_status))



#Combined summarizedExperiment object
combined_data_se <- GDCprepare(combined_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-combined-Dataset/")
combined_data_df <- as.data.frame(colData(combined_data_se))
combined_data_df$vital_status <- factor(combined_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
combined_data_df$vital_status <- as.numeric(as.character(combined_data_df$vital_status))


#Bulk data frame for combined merged data frame----
combined_bulk_rna_df <- combined_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(combined_bulk_rna_df) <- combined_data_se@colData@rownames
rownames(combined_bulk_rna_df) <- combined_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
combined_bulk_rna_df <- t(combined_bulk_rna_df)
combined_bulk_rna_df <- as.data.frame(combined_bulk_rna_df)
bulk_rownames <- rownames(combined_bulk_rna_df)
combined_bulk_rna_df$barcode <- bulk_rownames

combined_bulk_rna_df_unique <- subset(combined_bulk_rna_df, select = unique(colnames(combined_bulk_rna_df)))
combined_data_df_unique <- subset(combined_data_df, select = unique(colnames(combined_data_df)))
merged_df <- merge(combined_bulk_rna_df_unique, combined_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

merged_df$days_to_death <- ifelse(is.na(merged_df$days_to_death),0, merged_df$days_to_death)
merged_df$days_to_last_follow_up <- ifelse(is.na(merged_df$days_to_last_follow_up),0, merged_df$days_to_last_follow_up)

calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
calculated_days <- abs(calculated_days)
calculated_days[calculated_days==0]=1
cox_time <- calculated_days
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$tumor_stage
cox_tumor_n <- merged_df$ajcc_pathologic_n
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.n <- cox_tumor_n
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iv", replacement = 4)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iii", replacement = 3)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge ii", replacement = 2)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge i", replacement = 1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
cox_df <- filter(cox_df, !tumor.stage=="not reported")
cox_df <- cox_df[complete.cases(cox_df[, "ajcc.n"]), ]
cox_df <- cox_df[complete.cases(cox_df[, "vital.status"]), ]



#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
Read_data_se <- GDCprepare(read_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-read-Dataset/")
Read_data_df <- as.data.frame(colData(Read_data_se))
Read_data_df$vital_status <- factor(Read_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
Read_data_df$vital_status <- as.numeric(as.character(Read_data_df$vital_status))

#Bulk data frame for Read merged data frame----
bulk_rna_df <- read_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- read_data_se@colData@rownames
rownames(bulk_rna_df) <- read_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
read_data_df_unique <- subset(read_data_df, select = unique(colnames(read_data_df)))
merged_df <- merge(bulk_rna_df_unique, read_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

merged_df$days_to_death <- ifelse(is.na(merged_df$days_to_death),0, merged_df$days_to_death)
merged_df$days_to_last_follow_up <- ifelse(is.na(merged_df$days_to_last_follow_up),0, merged_df$days_to_last_follow_up)
read_df <- merged_df

merged_df <- read_df

calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
calculated_days <- abs(calculated_days)
calculated_days[calculated_days==0]=1
cox_time <- calculated_days
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$tumor_stage
cox_tumor_n <- merged_df$ajcc_pathologic_n
cox_gender <- merged_df$gender
cox_race <- merged_df$race
cox_eth <- merged_df$ethnicity
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.n <- cox_tumor_n
cox_df$race <- cox_race
cox_df$ethnicity <- cox_eth
cox_df$gender <- cox_gender
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iv", replacement = 4)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iii", replacement = 3)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge ii", replacement = 2)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge i", replacement = 1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
cox_df <- filter(cox_df, !tumor.stage=="not reported")
cox_df <- cox_df[complete.cases(cox_df[, "ajcc.n"]), ]


#Merged data frame for the normal/main colon cancer dataset
#Load the COAD data frame for use with Cox model----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

calculated_days <- merged_df$days.to.death - merged_df$days.to.last.follow.up
calculated_days[calculated_days==0]=1
merged_df$days.to.last.follow.up <- ifelse(is.na(calculated_days), merged_df$days.to.last.follow.up, calculated_days)
merged_df$days.to.last.follow.up <- ifelse(merged_df$days.to.last.follow.up==0, 1, merged_df$days.to.last.follow.up)
cox_time <- merged_df$days.to.last.follow.up
cox_event <- merged_df$vital.status
cox_tumor <- merged_df$tumor.stage
cox_tumor_m <- merged_df$ajcc.pathologic.m
cox_tumor_n <- merged_df$ajcc.pathologic.n
cox_gender <- merged_df$gender
cox_eth <- merged_df$ethnicity
cox_race <- merged_df$race
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.m <- cox_tumor_m
cox_df$ajcc.n <- cox_tumor_n
cox_df$race <- cox_race
cox_df$ethnicity <- cox_eth
cox_df$gender <- cox_gender
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iv", replacement = 4)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iii", replacement = 3)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge ii", replacement = 2)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge i", replacement = 1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
cox_df <- filter(cox_df, !tumor.stage=="not reported")
cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]



#List from the comparison paper----
compared_list <- list("SPOCK1"=0.287634, "VIM"=-0.70744, "C5AR"=0.376991, "WWTR1"=0.204543, "SERPINE1"=0.164739, "EFEMP1"=0.085466, "FSCN1"=0.11085, "FLNA"=-0.08217, "CXCL8"=0.18518, "NOX1"=0.04164)
compared_list <- as.numeric(compared_list)
names(compared_list) <- c("SPOCK1", "VIM", "C5AR1", "WWTR1", "SERPINE1", "EFEMP1", "FSCN1", "FLNA", "CXCL8", "NOX1")

#Cox model----
cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1
gene_sizes <- seq(from=100, to=2500, by=250)
for (y in gene_sizes){
  print(y)
  for (x in mirna_sde_optimized) {
    current_weight <- x
    current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = y, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
    cox_models[[as.character(counter)]] <- current_cox
    counter <- counter + 1
    
    #Storing all of the c-index values in a vector that we can use later to build the plot
    c_finder <-current_cox$CV$index[1]
    current_c <- current_cox$CV$cvm[c_finder]
    current_c <- round(current_c, digits = 4)
    my_cindicies <- c(my_cindicies, current_c)
    my_gene_sizes <- c(my_gene_sizes, y)
    
  }
}




#Just regular cox model, not changing gene size----
cox_models <- list()
my_cindicies <- c()
counter <- 1
for (x in mirna_sde_optimized[6]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight, tumor.stage = TRUE, tumor.n = FALSE, tumor.m = FALSE, regular.cox = TRUE, save.regular.cox.genes = TRUE, my.filename = "Data/Data-from-Cleaner-code/Regular_cox_model_outputs/coad_and_read_regular_cox_genes_tumor_1800.csv") 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
}

all_dfs <- list()
counter <- 1
for(x in cox_models){
  current_cox <- x
  current_df <- km_data_generator(data.source = current_cox$`Active Genes`, data.betas = current_cox$`Active Coefficients`)
  all_dfs[[as.character(counter)]] <- current_df
  counter <- counter + 1
}


#Just cox for the mirna metric at different sizes----
cox_models <- list()
my_cindicies <- c()
counter <- 1
for (x in mirna_master_list[1:14]) {
  current_mirnas <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = length(current_mirnas), cox.predictors = current_mirnas, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  
}



my_df <- data.frame(mirna_num=as.numeric(names(mirna_master_list)), concordance_index=my_cindicies)
my_df$cindex_se <- sqrt(my_df$concordance_index/sqrt(length(my_df$concordance_index)))
write.csv(my_df, file = "Data/Data-from-Cleaner-code/leukimia_patients_different_mirna_size_data_from_optimal_model.csv")
#my_df$mean_c_index <- mean(my_df$concordance_index)
df_to_plot <- data.frame(data=aggregate(x = coad_df$c_index,              
          by = list(coad_df$gene_size),              
          FUN = mean))    


colnames(df_to_plot) <- c("gene_num", "concordance_index")

des_res <- filter(des_res, gene_num<600)
my_df <- filter(my_df, gene_num<600)
all_methods_across_gene_size_df <- merge(des_res,my_df, by = "gene_num")
all_methods_across_gene_size_df <- read.csv("Data/Data-from-Cleaner-code/Graphs-for-c-index-performance/all_methods_across_gene_size_data.csv")
#all_methods_across_gene_size_df <- filter(all_methods_across_gene_size_df, method!="scDD")

coad_df <- filter(my_df, dataset=="TCGA-COAD")
coad_df <- select(coad_df, mirna_num:dataset)

read_df <- filter(my_df, dataset=="TCGA-READ")
read_df <- select(read_df,mirna_num:dataset)


#read only
read_gene_num <-ggplot(read_df, aes(x=mirna_num, y=concordance_index, group=1)) +
    geom_line(colour="#00BFC4")+
    geom_point(colour="#00BFC4")


res_aov_read <- aov(concordance_index~mirna_num, data = read_df)
aov_sum_read <- summary(res_aov_read)
pvalue_to_plot_read <- round(aov_sum_read[[1]][["Pr(>F)"]][1], digits = 5)

read_plot <- read_gene_num + ggtitle("TCGA-READ") +
  xlab("Number of MiRNAs Used in Model") +
  ylab("C-index")+
  ylim(0.40, 0.65)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("ANOVA, p=", pvalue_to_plot_read), x = 500, y = 0.42, hjust = 0, vjust = 1, size=4)
#scale_color_discrete(name = "Dataset", labels = c("TCGA-COAD", "TCGA-READ"))



#coad only
coad_gene_num <- ggplot(coad_df, aes(x=gene_size, y=c_index, group=1)) +
  geom_line(color="#f8766d")+
  geom_point(colour="#f8766d")



res_aov_coad <- aov(concordance_index~gene_size, data = coad_df)
aov_sum_coad <- summary(res_aov_coad)
pvalue_to_plot_coad <- round(aov_sum_coad[[1]][["Pr(>F)"]][1], digits = 5)



coad_plot <- coad_gene_num + ggtitle("TCGA-COAD") +
  xlab("Number of MiRNAs Used in Model") +
  ylab("C-index")+
  ylim(0.40, 0.65)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("ANOVA, p=", pvalue_to_plot_coad), x = 500, y = 0.42, hjust = 0, vjust = 1, size=4)
#scale_color_discrete(name = "Dataset", labels = c("TCGA-COAD", "TCGA-READ"))




#LC only
lc_gene_num <-ggplot(my_df, aes(x=mirna_num, y=concordance_index, group=1)) +
  geom_line(colour="#00BFC4")+
  geom_point(colour="#00BFC4")


res_aov_lc <- aov(concordance_index~mirna_num, data = my_df)
aov_sum_lc <- summary(res_aov_lc)
pvalue_to_plot_lc <- round(aov_sum_lc[[1]][["Pr(>F)"]][1], digits = 5)

lc_plot <- lc_gene_num + ggtitle("TCGA-LUAD") +
  xlab("Number of MiRNAs Used in Model") +
  ylab("C-index")+
  ylim(0.40, 0.66)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("ANOVA, p=", pvalue_to_plot_lc), x = 400, y = 0.42, hjust = 0, vjust = 1, size=4)
#scale_color_discrete(name = "Dataset", labels = c("TCGA-COAD", "TCGA-READ"))


#Cox for DESeq2----
resOrdered_subset_finished <- read.csv("Data/Data-from-Cleaner-code/deseq2_top1800_genes_cc_patients.csv")
cox_deseq2 <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = resOrdered_subset_finished$gene, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)

#Cox for edgeR----
edger_df_done <- read.csv("Data/Data-from-Cleaner-code/finished_edgeR_genes.csv")
clean_edger_names <- gene_vector_cleaner(edger_df_done$X)
edger_df_done$X <- clean_edger_names
colnames(edger_df_done)[1] <- "gene"
cox_edger2 <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = edger_df_done$gene, tumor.stage = FALSE, tumor.n =FALSE, tumor.m = FALSE)


#Leukemia only----
leuk_gene_num <-ggplot(my_df, aes(x=mirna_num, y=concordance_index, group=1)) +
  geom_line(colour="#f8766d")+
  geom_point(colour="#f8766d")


res_aov_leuk <- aov(concordance_index~mirna_num, data = my_df)
aov_sum_leuk <- summary(res_aov_leuk)
pvalue_to_plot_leuk <- round(aov_sum_leuk[[1]][["Pr(>F)"]][1], digits = 5)

leuk_plot <- leuk_gene_num + ggtitle("TCGA-DLBC") +
  xlab("Number of MiRNAs Used in Model") +
  ylab("C-index")+
  ylim(0.40, 0.67)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("ANOVA, p=", pvalue_to_plot_leuk), x = 500, y = 0.42, hjust = 0, vjust = 1, size=4)
#scale_color_discrete(name = "Dataset", labels = c("TCGA-COAD", "TCGA-READ"))



#both DLBC and LUAD on same plot
both_leuk_and_luad_data <- read.csv("Data/Data-from-Cleaner-code/Graphs-for-c-index-performance/Different-mirna-sizes/leuk_and_luad_datasets.csv")
both_other_dataset_plot <- ggplot(both_leuk_and_luad_data, aes(x=mirna_num, y=concordance_index, group=dataset)) +
  geom_line(aes(color=dataset))+
  geom_point(aes(color=dataset))


res_aov_both <- aov(concordance_index~dataset, data = both_leuk_and_luad_data)
aov_sum_both <- summary(res_aov_both)
pvalue_to_plot_both <- round(aov_sum_both[[1]][["Pr(>F)"]][1], digits = 5)

both_plot <- both_other_dataset_plot + ggtitle("TCGA-DLBC vs. TCGA-LUAD") +
  xlab("Number of MiRNAs Used in Model") +
  ylab("C-index")+
  ylim(0.40, 0.67)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("ANOVA, p=", pvalue_to_plot_both), x = 500, y = 0.43, hjust = 0, vjust = 1, size=4)+
  scale_color_discrete(name = "Dataset", labels = c("TCGA-DLBC", "TCGA-LUAD"))



#both COAD and READ on same plot
both_gene_num <- ggplot(my_df, aes(x=mirna_num, y=concordance_index, group=dataset)) +
  geom_line(aes(color=dataset))+
  geom_point(aes(color=dataset))



res_aov_both <- aov(concordance_index~dataset, data = my_df)
aov_sum_both <- summary(res_aov_both)
pvalue_to_plot_both <- round(aov_sum_both[[1]][["Pr(>F)"]][1], digits = 5)



both_plot <- both_gene_num + ggtitle("TCGA-COAD vs. TCGA-READ") +
  xlab("Number of MiRNAs Used in Model") +
  ylab("C-index")+
  ylim(0.40, 0.65)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("ANOVA, p=", pvalue_to_plot_both), x = 1000, y = 0.43, hjust = 0, vjust = 1, size=4)+
  scale_color_discrete(name = "Dataset", labels = c("TCGA-COAD", "TCGA-READ"))







#Figure of the mirna size panels----
mirna_figure <- ggarrange(both_plot, ggarrange(coad_plot, read_plot, ncol = 2, labels = c("B.", "C.")), nrow = 2, labels = "A.")
mirna_figure2 <- ggarrange(both_plot, ggarrange(lc_plot, leuk_plot, ncol = 2, labels = c("B.", "C.")), nrow = 2, labels = "A.")

#Code to plot different pseudotime markers----
pt_different_markers_coad <- read.csv("~/Desktop/different_sdes_markers_for_server/different_sdes_markers_cc_patients.csv")
pt_different_markers_coad$method <- c("snai1", "twist1", "twist1-redo", "zeb1", "vim")


res_aov <- aov(concordance_index~method, data = pt_different_markers)
aov_sum <- summary(res_aov)
pvalue_to_plot <- round(aov_sum[[1]][["Pr(>F)"]][1], digits = 5)


pt_average_df <- merge(pt_different_markers, pt_different_markers_coad, by="method")

pt_average_df$mean_performance <- c(0.58065, 0.5897, 0.56025, 0.5897, 0.56025)


res_aov <- aov(mean_performance~method, data = pt_average_df)
aov_sum <- summary(res_aov)
pvalue_to_plot <- round(aov_sum[[1]][["Pr(>F)"]][1], digits = 5)



p <- ggplot(pt_average_df, aes(x=method, y=mean_performance)) + geom_point()
p+ ggtitle("Concordance Index Performance for Different SDES Markers CC & RC Patients") + xlab("Marker") + ylab("Mean Concordance Index") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), panel.background = element_blank())


for(x in cox_models){
  print(x$CV)
}


#Simple code to get the names of the active genes from the cox model
Active.Index <- which(as.logical(cox_models$`6`$Coefficients) != 0)
Active.Coefficients  <- cox_models$`6`$Coefficients[Active.Index]
active_genes <-rownames(cox_models$`6`$Coefficients)[Active.Index]
print(active_genes)


#Cox model with just the active genes that were previously identified----
active_genes_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = length(active_genes), cox.predictors = active_genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)

#Cox for individual metrics---
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)


#Cox for scDD genes across gene sizes----
scdd_genes <- scdd_res
scdd_genes <- intersect(scdd_genes$gene, colnames(cox_df))


cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1
gene_sizes <- seq(from=100, to=500, by=50)
for (y in gene_sizes) {
  print(y)
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = y, cox.predictors = scdd_genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  my_gene_sizes <- c(my_gene_sizes, y)
}


my_df <- data.frame(gene_num=my_gene_sizes, concordance_index=my_cindicies)
my_df$cindex_se <- sqrt(my_df$concordance_index/sqrt(length(my_df$concordance_index)))
write.csv(my_df, file = "Data/Data-from-Cleaner-code/mirna_sde_cc_patients_gene_size_data.csv")
df_to_plot <- data.frame(data=aggregate(x = my_df$concordance_index,              
                                        by = list(my_df$gene_num),              
                                        FUN = mean))    


colnames(df_to_plot) <- c("gene_num", "concordance_index")


#df_to_plot$gene_num <- factor(df_to_plot$gene_num)
res_aov <- aov(concordance_index~gene_num, data = df_to_plot)
aov_sum <- summary(res_aov)
pvalue_to_plot <- round(aov_sum[[1]][["Pr(>F)"]][1], digits = 5)


p<-ggplot(df_to_plot, aes(x=gene_num, y=concordance_index, group=1)) +
  geom_line(aes(color="red"))+
  geom_point()

p + ggtitle("scDD Concordance Index Across Gene Number") +
  xlab("Number of Genes Used in Model") +
  ylab("Mean Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("P=", pvalue_to_plot), x = 400, y = 0.58, hjust = 0, vjust = 1, size=5)



#Just for testing scDD active genes----
scdd_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 511, cox.predictors = scdd_genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)


#Cox for DESingle----
desingle_genes <- intersect(rownames(des_results), colnames(cox_df))


cox_models <- list()
my_cindicies <- c()
my_gene_sizes <- c()
counter <- 1
gene_sizes <- seq(from=100, to=1800, by=50)
for (y in gene_sizes) {
  print(y)
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = y, cox.predictors = scdd_genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  my_gene_sizes <- c(my_gene_sizes, y)
}


my_df <- data.frame(gene_num=my_gene_sizes, concordance_index=my_cindicies)
my_df$cindex_se <- sqrt(my_df$concordance_index/sqrt(length(my_df$concordance_index)))
write.csv(my_df, file = "Data/Data-from-Cleaner-code/desingle_cc_patients_gene_size_data.csv")
df_to_plot <- data.frame(data=aggregate(x = my_df$concordance_index,              
                                        by = list(my_df$gene_num),              
                                        FUN = mean))    


colnames(df_to_plot) <- c("gene_num", "concordance_index")


#df_to_plot$gene_num <- factor(df_to_plot$gene_num)
res_aov <- aov(concordance_index~gene_num, data = df_to_plot)
aov_sum <- summary(res_aov)
pvalue_to_plot <- round(aov_sum[[1]][["Pr(>F)"]][1], digits = 5)


p<-ggplot(df_to_plot, aes(x=gene_num, y=concordance_index, group=1)) +
  geom_line(aes(color="red"))+
  geom_point()

p + ggtitle("DEsingle Concordance Index Across Gene Number") +
  xlab("Number of Genes Used in Model") +
  ylab("Mean Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))+
  annotate(geom = 'text', label = paste0("P=", pvalue_to_plot), x = 1500, y = 0.58, hjust = 0, vjust = 1, size=5)









desingle_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = desingle_genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)




#Subset to a particular race
cox_df <- subset(cox_df, cox_df$race=="black or african american")




#Filtering the merged df to just stage I and II patients----
load("Data/Exported-data/R-objects/COA_data_df.RData", verbose = TRUE)
COA_data_df_sub <- filter(COA_data_df, tumor_stage=="stage i" | tumor_stage=="stage ia" | tumor_stage=="stage ii" | tumor_stage=="stage iia" | tumor_stage=="stage iib" | tumor_stage=="stage iic")
COA_data_df_sub$tumor_stage <- gsub(COA_data_df_sub$tumor_stage, pattern="a", replacement="")
COA_data_df_sub$tumor_stage <- gsub(COA_data_df_sub$tumor_stage, pattern="b", replacement="")
COA_data_df_sub$tumor_stage <- gsub(COA_data_df_sub$tumor_stage, pattern="c", replacement="")
merged_df_stageI_II <- filter(merged_df, )

#KM curves----
hr_calcs <- list()
counter <- 1
for (x in cox_models[1:11]) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = FALSE)
  hr_calcs[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#Cox for individual metrics----
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE) 
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE) 
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE) 

#HR calcs for individual metrics---
hr_mad <- hr_calculator(model.coefs = cox_mad$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)
hr_sdes <- hr_calculator(model.coefs = cox_sdes$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)
hr_mirna <- hr_calculator(model.coefs = cox_mirna$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)

#HR calcs for just active genes----
hr_active_genes <- hr_calculator(model.coefs = active_genes_cox$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE) 

#HR calcs for scDD----
hr_scdd <- hr_calculator(model.coefs = scdd_cox$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE)

#HR calcs for DESingle----
hr_desingle <- hr_calculator(model.coefs = desingle_cox$Coefficients, data = cox_df, include.cat.data = FALSE, tumor.stage = FALSE, n.stage = FALSE)



#HR calcs for DESeq2----
hr_deseq2 <- hr_calculator(model.coefs = cox_deseq2$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)

#HR calcs for edgeR----
hr_edger <- hr_calculator(model.coefs = cox_edger$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)

#KM p-values----
hr_pvalues <- list()
counter <- 1
for (x in hr_calcs) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#P-values for HR calcs for individual metrics---
pvalue_mad <- km_pvalue_calculator(surv.time = hr_mad$DF$days.to.last.follow.up, surv.status = hr_mad$DF$vital.status, surv.predictors = hr_mad$DF$risk, surv.df = hr_mad$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_sdes <- km_pvalue_calculator(surv.time = hr_sdes$DF$days.to.last.follow.up, surv.status = hr_sdes$DF$vital.status, surv.predictors = hr_sdes$DF$risk, surv.df = hr_sdes$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_mirna <- km_pvalue_calculator(surv.time = hr_mirna$DF$days.to.last.follow.up, surv.status = hr_mirna$DF$vital.status, surv.predictors = hr_mirna$DF$risk, surv.df = hr_mirna$DF, num.sig.figs = 4, scientific.p.value = TRUE)

#P-value for scDD----
pvalue_scdd <- km_pvalue_calculator(surv.time = hr_scdd$DF$days.to.last.follow.up, surv.status = hr_scdd$DF$vital.status, surv.predictors = hr_scdd$DF$risk, surv.df = hr_scdd$DF, num.sig.figs = 4, scientific.p.value = TRUE)
 

#P-value for DESingle----
pvalue_desingle <- km_pvalue_calculator(surv.time = hr_desingle$DF$days.to.last.follow.up, surv.status = hr_desingle$DF$vital.status, surv.predictors = hr_desingle$DF$risk, surv.df = hr_desingle$DF, num.sig.figs = 4, scientific.p.value = TRUE)



#P-value for DESeq2----
pvalue_deseq2 <- km_pvalue_calculator(surv.time = hr_deseq2$DF$days.to.last.follow.up, surv.status = hr_deseq2$DF$vital.status, surv.predictors = hr_deseq2$DF$risk, surv.df = hr_deseq2$DF, num.sig.figs = 4, scientific.p.value = TRUE)

#P-value for edgeR----
pvalue_edger <- km_pvalue_calculator(surv.time = hr_edger$DF$days.to.last.follow.up, surv.status = hr_edger$DF$vital.status, surv.predictors = hr_edger$DF$risk, surv.df = hr_edger$DF, num.sig.figs = 4, scientific.p.value = TRUE)

#P-value for active genes----
pvalue_active <- km_pvalue_calculator(surv.time = hr_active_genes$DF$days.to.last.follow.up, surv.status = hr_active_genes$DF$vital.status, surv.predictors = hr_active_genes$DF$risk, surv.df = hr_active_genes$DF, num.sig.figs = 4, scientific.p.value = TRUE)

#KM plots----
hr_plots <- list()
counter <- 1
for (x in hr_calcs) {
  current_calc <- x
  current_pvalue <- hr_pvalues[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES + Tumor Stage Rectal Cancer Patients")
  hr_plots[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}


#KM Plots for HR calcs for individual metrics---
km_mad <- km_plotter(km.fit = hr_mad$KM, data.source = hr_mad$DF, p.value = pvalue_mad, plot.title = "MAD + Tumor Stage + N Stage Rectal Cancer Patients")
km_sdes <- km_plotter(km.fit = hr_sdes$KM, data.source = hr_sdes$DF, p.value = pvalue_sdes, plot.title = "SDES + Tumor Stage + N Stage Rectal Cancer Patients")
km_mirna <- km_plotter(km.fit = hr_mirna$KM, data.source = hr_mirna$DF, p.value = pvalue_mirna, plot.title = "MiRNA + Tumor Stage + N Stage Rectal Cancer Patients")

#KM plot for scDD----
km_scdd <- km_plotter(km.fit = hr_scdd$KM, data.source = hr_scdd$DF, p.value = pvalue_scdd, plot.title = "scDD Colon Cancer Patients")


#KM plot for DESingle----
km_desingle <- km_plotter(km.fit = hr_desingle$KM, data.source = hr_desingle$DF, p.value = pvalue_desingle, plot.title = "DESingle Rectal Cancer Patients")


#KM plot for DESeq2----
km_deseq2 <- km_plotter(km.fit = hr_deseq2$KM, data.source = hr_deseq2$DF, p.value = pvalue_deseq2, plot.title = "DESeq2 + Tumor Stage + N Stage Colon & Rectal Cancer Patients")

#KM plot for edgeR----
km_edger <- km_plotter(km.fit = hr_edger$KM, data.source = hr_edger$DF, p.value = TRUE, plot.title = "edgeR Rectal Cancer Patients")

#KM plot for active genes----
km_active <- km_plotter(km.fit = hr_active_genes$KM, data.source = hr_active_genes$DF, p.value = pvalue_active, plot.title = "Active Genes Only SDES + MiRNA Leukemia Patients from K562 Cell-line")


#K562 cell-line analysis

#Now denoising the sc-data----
k562_cells <- magic_denoiser(sc.data = k562_cells,magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = "VIM", cell.data = k562_cells$denoised_sc_dataframe, cell.meta = k562_cells$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
save(mad.genes, file = "k562_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(k562_cells$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "k562_sde.RData")

#Mirna metric----
counter <- 1
mirna_master_list <- list()

for (y in mirna_sizes) {
  print(y)
  current_mirna <- mirna_calculator(ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "gastric cancer", mirna.filename = NULL, mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"), max.mirnas = y)
  mirna_master_list[[as.character(counter)]] <- current_mirna
  counter <- counter+1
}

mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "gastric cancer", mirna.filename = "k562_cells_mirnas_dbdemc_up.RData")
save(mirna.genes, file = "k562_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_k562.RData")


#Optimizing the mirna + MAD metric----
mirna_mad_optimized <- two_weight_optimizer(first.metric = mad.genes, second.metric = mirna.genes, my.filename = "Data/Data-from-Cleaner-code/mad_mirna_optimized_k562.RData")


#Leukimia TCGA bulk data----
leuk_query <- GDCquery(project       = "TCGA-DLBC",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = leuk_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-blymph-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
Leuk_data_se <- GDCprepare(leuk_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-blymph-Dataset/")
Leuk_data_df <- as.data.frame(colData(Leuk_data_se))
Leuk_data_df$vital_status <- factor(Leuk_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
Leuk_data_df$vital_status <- as.numeric(as.character(Leuk_data_df$vital_status))

#Bulk dataframe for Leuk merged dataframe----
bulk_rna_df <- Leuk_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- Leuk_data_se@colData@rownames
rownames(bulk_rna_df) <- Leuk_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
Leuk_data_df_unique <- subset(Leuk_data_df, select = unique(colnames(Leuk_data_df)))
merged_df <- merge(bulk_rna_df_unique, Leuk_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

merged_df$days_to_death <- ifelse(is.na(merged_df$days_to_death),0, merged_df$days_to_death)
merged_df$days_to_last_follow_up <- ifelse(is.na(merged_df$days_to_last_follow_up),0, merged_df$days_to_last_follow_up)
calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
calculated_days <- abs(calculated_days)
calculated_days[calculated_days==0]=1
cox_time <- calculated_days
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$tumor_stage
cox_df <- subset(merged_df, select=c(TSPAN6:V56404))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event

#Cox model----
cox_models_k562 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models_k562[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Cox for individual metrics---
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)




#KM curves----
hr_calcs_k562 <- list()
counter <- 1
for (x in cox_models_k562[6:10]) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df)
  hr_calcs_k562[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#HR calcs for individual metrics---
hr_mad <- hr_calculator(model.coefs = cox_mad$Coefficients, data = cox_df)
hr_sdes <- hr_calculator(model.coefs = cox_sdes$Coefficients, data = cox_df)
hr_mirna <- hr_calculator(model.coefs = cox_mirna$Coefficients, data = cox_df)


#KM p-values----
hr_pvalues_k562 <- list()
counter <- 1
for (x in hr_calcs_k562) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_k562[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#P-values for HR calcs for individual metrics---
pvalue_mad <- km_pvalue_calculator(surv.time = hr_mad$DF$days.to.last.follow.up, surv.status = hr_mad$DF$vital.status, surv.predictors = hr_mad$DF$risk, surv.df = hr_mad$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_sdes <- km_pvalue_calculator(surv.time = hr_sdes$DF$days.to.last.follow.up, surv.status = hr_sdes$DF$vital.status, surv.predictors = hr_sdes$DF$risk, surv.df = hr_sdes$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_mirna <- km_pvalue_calculator(surv.time = hr_mirna$DF$days.to.last.follow.up, surv.status = hr_mirna$DF$vital.status, surv.predictors = hr_mirna$DF$risk, surv.df = hr_mirna$DF, num.sig.figs = 4, scientific.p.value = TRUE)


#KM plots----
hr_plots_k562 <- list()
counter <- 1
for (x in hr_calcs_k562) {
  current_calc <- x
  current_pvalue <- hr_pvalues_k562[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES K562 with TCGA B-cell lymphoma")
  hr_plots_k562[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}

#KM Plots for HR calcs for individual metrics---
km_mad <- km_plotter(km.fit = hr_mad$KM, data.source = hr_mad$DF, p.value = pvalue_mad, plot.title = "MAD Metric Alone")
km_sdes <- km_plotter(km.fit = hr_sdes$KM, data.source = hr_sdes$DF, p.value = pvalue_sdes, plot.title = "SDES Metric Alone")
km_mirna <- km_plotter(km.fit = hr_mirna$KM, data.source = hr_mirna$DF, p.value = pvalue_mirna, plot.title = "MiRNA Metric Alone")




#Glio analysis
#Now denoising the sc-data----
glio <- magic_denoiser(sc.data = glio,magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = "VIM", cell.data = glio$denoised_sc_dataframe, cell.meta = glio$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(glio$denoised_sc_dataframe)
save(mad.genes, file = "glio_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(glio$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "glio_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "brain cancer", mirna.filename = "glio_mirna_targets.RData", mirna.remove = "hsa-miR-129-1-3p")
save(mirna.genes, file = "glio_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_glio.RData")

#Glioblastoma TCGA bulk data----
glio_query <- GDCquery(project        = "TCGA-GBM",
                        data.category = "Transcriptome Profiling",
                        data.type     = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = glio_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-Glio-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
Glio_data_se <- GDCprepare(glio_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-Glio-Dataset/")
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
cox_df <- subset(merged_df, select=c(TSPAN6:V56404))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df <- cox_df %>% filter(!is.na(vital.status))
cox_df$days.to.last.follow.up <- as.numeric(cox_df$days.to.last.follow.up)
cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up==0,1,cox_df$days.to.last.follow.up)



#Cox model----
cox_models_glio <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 
  cox_models_glio[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Cox for individual metrics---
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE)

#Cox for just scDD ----
scdd_genes <- scdd_res
scdd_genes <- intersect(scdd_genes$gene, colnames(cox_df))
cox_scdd_glio <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = , tumor.stage = FALSE, tumor.m = FALSE, tumor.n = FALSE) 

#KM curves----
merged_df <- merged_df %>% filter(!is.na(vital_status))
colnames(merged_df)[colnames(merged_df) == "vital_status"] <- "vital.status"
colnames(merged_df)[colnames(merged_df) == "days_to_last_follow_up"] <- "days.to.last.follow.up"

hr_calcs_glio <- list()
counter <- 1
for (x in cox_models_glio) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df)
  hr_calcs_glio[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#HR calcs for individual metrics---
hr_mad <- hr_calculator(model.coefs = cox_mad$Coefficients, data = cox_df)
hr_sdes <- hr_calculator(model.coefs = cox_sdes$Coefficients, data = cox_df)
hr_mirna <- hr_calculator(model.coefs = cox_mirna$Coefficients, data = cox_df)


#KM p-values----
hr_pvalues_glio <- list()
counter <- 1
for (x in hr_calcs_glio) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_glio[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#P-values for HR calcs for individual metrics---
pvalue_mad <- km_pvalue_calculator(surv.time = hr_mad$DF$days.to.last.follow.up, surv.status = hr_mad$DF$vital.status, surv.predictors = hr_mad$DF$risk, surv.df = hr_mad$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_sdes <- km_pvalue_calculator(surv.time = hr_sdes$DF$days.to.last.follow.up, surv.status = hr_sdes$DF$vital.status, surv.predictors = hr_sdes$DF$risk, surv.df = hr_sdes$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_mirna <- km_pvalue_calculator(surv.time = hr_mirna$DF$days.to.last.follow.up, surv.status = hr_mirna$DF$vital.status, surv.predictors = hr_mirna$DF$risk, surv.df = hr_mirna$DF, num.sig.figs = 4, scientific.p.value = TRUE)


#KM plots----
hr_plots_glio <- list()
counter <- 1
for (x in hr_calcs_glio) {
  current_calc <- x
  current_pvalue <- hr_pvalues_glio[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES Glio")
  hr_plots_glio[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}





#A549 cell-line analysis
#Loading needed function----
gene_name_cleaner <- function(data.to.clean=all_tumor_cells_fpkm_denoised_df){
  data.to.clean <-t(data.to.clean)
  current_colname_split <- strsplit(colnames(data.to.clean), "_")
  finished_gene_list <- c()
  current_list <- current_colname_split
  for (y in seq(1:length(current_list))){
    finished_gene_list <- c(finished_gene_list, current_list[[y]][2])
  }
  colnames(data.to.clean) <- finished_gene_list
  return(data.to.clean)
}

#More pre-processing
rownames(cc_cell_line_fpkm) <- cc_cell_line_fpkm$X
cc_cell_line_fpkm<- subset(cc_cell_line_fpkm, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cc_cell_line_fpkm <- gene_name_cleaner(data.to.clean = cc_cell_line_fpkm)
cc_cell_colnames <- unique(colnames(cc_cell_line_fpkm))
cc_cell_line_fpkm <- subset(cc_cell_line_fpkm, select=cc_cell_colnames)
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_line_fpkm <- as.data.frame(cc_cell_line_fpkm)
a549_cells <- select(cc_cell_line_fpkm, contains("A549"))


#Now denoising the sc-data----
a549_cells <- apply(a549_cells, c(1,2), as.numeric)
a549_cells <- magic_denoiser(sc.data = a549_cells, magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = "VIM", cell.data = a549_cells$denoised_sc_dataframe, cell.meta = a549_cells$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(a549_cells$denoised_sc_dataframe)
save(mad.genes, file = "a549_cells_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(a549_cells$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "a549_cells_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "lung cancer", mirna.filename = "TargetScane_a549_cell_mirna.RData", mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"), print.ts.targets = TRUE)
save(mirna.genes, file = "a549_cells_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_a549_cells.RData")

#Lung cancer TCGA bulk data----
lung_query <- GDCquery(project       = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = lung_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-lungluad-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
lung_data_se <- GDCprepare(lung_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-lungluad-Dataset/")
lung_data_df <- as.data.frame(colData(lung_data_se))
lung_data_df$vital_status <- factor(lung_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
lung_data_df$vital_status <- as.numeric(as.character(lung_data_df$vital_status))

#Bulk dataframe for lung merged dataframe----
bulk_rna_df <- lung_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- lung_data_se@colData@rownames
rownames(bulk_rna_df) <- lung_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
lung_data_df_unique <- subset(lung_data_df, select = unique(colnames(lung_data_df)))
merged_df <- merge(bulk_rna_df_unique, lung_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

calculated_days <- merged_df$days_to_death - merged_df$days_to_last_follow_up
calculated_days[calculated_days==0]=1
merged_df <- merged_df[complete.cases(merged_df[, "days_to_last_follow_up"]), ]
merged_df$days_to_last_follow_up <- ifelse(merged_df$days_to_last_follow_up==0,1,merged_df$days_to_last_follow_up)
cox_time <- merged_df$days_to_last_follow_up
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$tumor_stage
cox_tumor_n <- merged_df$ajcc_pathologic_n
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.n <- cox_tumor_n
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df <- filter(cox_df, !tumor.stage=="not reported")
cox_df <- filter(cox_df, !ajcc.n=="NX")
cox_df <- cox_df[complete.cases(cox_df[, "ajcc.n"]), ]

#Cox model----
cox_models_a549 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE) 
  cox_models_a549[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Cox for individual metrics---
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)

#KM curves----
hr_calcs_a549 <- list()
counter <- 1
for (x in cox_models_a549[1:11]) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)
  hr_calcs_a549[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#HR calcs for individual metrics---
hr_mad <- hr_calculator(model.coefs = cox_mad$Coefficients, data = cox_df)
hr_sdes <- hr_calculator(model.coefs = cox_sdes$Coefficients, data = cox_df)
hr_mirna <- hr_calculator(model.coefs = cox_mirna$Coefficients, data = cox_df)


#KM p-values----
hr_pvalues_a549 <- list()
counter <- 1
for (x in hr_calcs_a549) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_a549[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#P-values for HR calcs for individual metrics---
pvalue_mad <- km_pvalue_calculator(surv.time = hr_mad$DF$days.to.last.follow.up, surv.status = hr_mad$DF$vital.status, surv.predictors = hr_mad$DF$risk, surv.df = hr_mad$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_sdes <- km_pvalue_calculator(surv.time = hr_sdes$DF$days.to.last.follow.up, surv.status = hr_sdes$DF$vital.status, surv.predictors = hr_sdes$DF$risk, surv.df = hr_sdes$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_mirna <- km_pvalue_calculator(surv.time = hr_mirna$DF$days.to.last.follow.up, surv.status = hr_mirna$DF$vital.status, surv.predictors = hr_mirna$DF$risk, surv.df = hr_mirna$DF, num.sig.figs = 4, scientific.p.value = TRUE)



#KM plots----
hr_plots_a549 <- list()
counter <- 1
for (x in hr_calcs_a549) {
  current_calc <- x
  current_pvalue <- hr_pvalues_a549[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES + Tumor Stage + N Stage A549 Cell-line with Lung AUD Bulk Dataset")
  hr_plots_a549[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}



#HCT116 cell-line analysis----
#More pre-processing----
rownames(cc_cell_line_fpkm) <- cc_cell_line_fpkm$X
cc_cell_line_fpkm<- subset(cc_cell_line_fpkm, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cc_cell_line_fpkm <- gene_name_cleaner(data.to.clean = cc_cell_line_fpkm)
cc_cell_colnames <- unique(colnames(cc_cell_line_fpkm))
cc_cell_line_fpkm <- subset(cc_cell_line_fpkm, select=cc_cell_colnames)
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_line_fpkm <- as.data.frame(cc_cell_line_fpkm)
hct116_cells <- select(cc_cell_line_fpkm, contains("HCT116"))


#Logged data----
rownames(cc_cell_line_fpkm_logged) <- cc_cell_line_fpkm_logged$X
cc_cell_line_fpkm_logged<- subset(cc_cell_line_fpkm_logged, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cc_cell_line_fpkm_logged <- gene_name_cleaner(data.to.clean = cc_cell_line_fpkm_logged)
cc_cell_colnames <- unique(colnames(cc_cell_line_fpkm_logged))
cc_cell_line_fpkm_logged <- subset(cc_cell_line_fpkm_logged, select=cc_cell_colnames)
cc_cell_line_fpkm_logged <- t(cc_cell_line_fpkm_logged)
cc_cell_line_fpkm_logged <- as.data.frame(cc_cell_line_fpkm_logged)
hct116_cells <- select(cc_cell_line_fpkm_logged, contains("HCT116"))


#Now denoising the sc-data----
hct116_cells <- apply(hct116_cells, c(1,2), as.numeric)
hct116_cells <- magic_denoiser(sc.data = hct116_cells, magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = NULL, cell.data = hct116_cells$denoised_sc_dataframe, cell.meta = hct116_cells$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(hct116_cells$denoised_sc_dataframe)
save(mad.genes, file = "Data/Data-from-Cleaner-code/hct116_cells_mad_logged.RData")

#Switchde metric----
sde.genes <- switchde_calculator(hct116_cells$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "Data/Data-from-Cleaner-code/hct116_cells_sde_logged.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "colon cancer", cancer.type2 = "colorectal cancer", mirna.filename = "TargetScane_hct116_cell_mirna_logged.RData", print.ts.targets = TRUE, mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"))
save(mirna.genes, file = "Data/Data-from-Cleaner-code/hct116_cells_mirna_logged.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_hct116_cells_logged.RData")


#Optimizing the mirna + MAD metric----
mirna_mad_optimized <- two_weight_optimizer(first.metric = mad.genes, second.metric = mirna.genes, my.filename = "mad_mirna_optimized_hct116_cells.RData")


load("Data/Exported-data/R-objects/merged_df_replaced.RData", verbose = TRUE)
calculated_days <- merged_df$days.to.death - merged_df$days.to.last.follow.up
calculated_days[calculated_days==0]=1
merged_df$days.to.last.follow.up <- ifelse(is.na(calculated_days), merged_df$days.to.last.follow.up, calculated_days)
merged_df$days.to.last.follow.up <- ifelse(merged_df$days.to.last.follow.up==0, 1, merged_df$days.to.last.follow.up)


calculated_days <- merged_df$days.to.death - merged_df$days.to.last.follow.up
calculated_days[calculated_days==0]=1
merged_df$days.to.last.follow.up <- ifelse(is.na(calculated_days), merged_df$days.to.last.follow.up, calculated_days)
merged_df$days.to.last.follow.up <- ifelse(merged_df$days.to.last.follow.up==0, 1, merged_df$days.to.last.follow.up)
cox_time <- merged_df$days.to.last.follow.up
cox_event <- merged_df$vital.status
cox_tumor <- merged_df$tumor.stage
cox_tumor_n <- merged_df$ajcc.pathologic.n
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
cox_df <- filter(cox_df, !tumor.stage=="not reported")

#Cox model----
cox_models_hct116 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors =current_weight, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE) 
  cox_models_hct116[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Cox for individual metrics---
cox_mad <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mad.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)
cox_sdes <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = sde.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)
cox_mirna <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = mirna.genes, tumor.stage = TRUE, tumor.m = TRUE, tumor.n = FALSE)


#KM curves----
hr_calcs_hct116 <- list()
counter <- 1
for (x in cox_models_hct116[1:11]) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df, include.cat.data = TRUE, tumor.stage = TRUE, n.stage = TRUE)
  hr_calcs_hct116[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#HR calcs for individual metrics---
hr_mad <- hr_calculator(model.coefs = cox_mad$Coefficients, data = cox_df)
hr_sdes <- hr_calculator(model.coefs = cox_sdes$Coefficients, data = cox_df)
hr_mirna <- hr_calculator(model.coefs = cox_mirna$Coefficients, data = cox_df)


#KM p-values----
hr_pvalues_hct116 <- list()
counter <- 1
for (x in hr_calcs_hct116) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_hct116[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#P-values for HR calcs for individual metrics---
pvalue_mad <- km_pvalue_calculator(surv.time = hr_mad$DF$days.to.last.follow.up, surv.status = hr_mad$DF$vital.status, surv.predictors = hr_mad$DF$risk, surv.df = hr_mad$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_sdes <- km_pvalue_calculator(surv.time = hr_sdes$DF$days.to.last.follow.up, surv.status = hr_sdes$DF$vital.status, surv.predictors = hr_sdes$DF$risk, surv.df = hr_sdes$DF, num.sig.figs = 4, scientific.p.value = TRUE)
pvalue_mirna <- km_pvalue_calculator(surv.time = hr_mirna$DF$days.to.last.follow.up, surv.status = hr_mirna$DF$vital.status, surv.predictors = hr_mirna$DF$risk, surv.df = hr_mirna$DF, num.sig.figs = 4, scientific.p.value = TRUE)


#KM plots----
hr_plots_hct116 <- list()
counter <- 1
for (x in hr_calcs_hct116) {
  current_calc <- x
  current_pvalue <- hr_pvalues_hct116[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES + Tumor Stage + N Stage HCT116 Cell-line Log Transformed Data")
  hr_plots_hct116[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}

#KM Plots for HR calcs for individual metrics---
km_mad <- km_plotter(km.fit = hr_mad$KM, data.source = hr_mad$DF, p.value = pvalue_mad, plot.title = "MAD + Tumor Stage + N Stage HCT116 Log Transformed Data Metric Alone")
km_sdes <- km_plotter(km.fit = hr_sdes$KM, data.source = hr_sdes$DF, p.value = pvalue_sdes, plot.title = "SDES + Tumor Stage + N Stage HCT116 Log Transformed Data Metric Alone")
km_mirna <- km_plotter(km.fit = hr_mirna$KM, data.source = hr_mirna$DF, p.value = pvalue_mirna, plot.title = "MiRNA + Tumor Stage + N Stage HCT116 Log Transformed Data Metric Alone")


#H1437 cell-line analysis---
#Loading needed function----
gene_name_cleaner <- function(data.to.clean=all_tumor_cells_fpkm_denoised_df){
  data.to.clean <-t(data.to.clean)
  current_colname_split <- strsplit(colnames(data.to.clean), "_")
  finished_gene_list <- c()
  current_list <- current_colname_split
  for (y in seq(1:length(current_list))){
    finished_gene_list <- c(finished_gene_list, current_list[[y]][2])
  }
  colnames(data.to.clean) <- finished_gene_list
  return(data.to.clean)
}

#More pre-processing
rownames(cc_cell_line_fpkm) <- cc_cell_line_fpkm$X
cc_cell_line_fpkm<- subset(cc_cell_line_fpkm, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cc_cell_line_fpkm <- gene_name_cleaner(data.to.clean = cc_cell_line_fpkm)
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_rownames <- unique(rownames(cc_cell_line_fpkm))
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_line_fpkm <- subset(cc_cell_line_fpkm, select=cc_cell_rownames)
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_line_fpkm <- as.data.frame(cc_cell_line_fpkm)
h1437_cells <- select(cc_cell_line_fpkm, contains("H1437"))


#Now denoising the sc-data----
h1437_cells <- apply(h1437_cells, c(1,2), as.numeric)
h1437_cells <- magic_denoiser(sc.data = h1437_cells, magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP"), cell.data = h1437_cells$denoised_sc_dataframe, cell.meta = h1437_cells$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(h1437_cells$denoised_sc_dataframe)
save(mad.genes, file = "h1437_cells_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(h1437_cells$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "h1437_cells_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "lung cancer", mirna.filename = "TargetScane_h1437_cell_mirna.RData", mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"), print.ts.targets = TRUE)
save(mirna.genes, file = "h1437_cells_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_h1437_cells.RData")

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/df_for_train_test_split.RData")

#Cox model----
cox_models_h1437 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = df_for_train_test_split, gene.num = 900, cox.predictors = current_weight) 
  cox_models_h1437[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#KM curves----
hr_calcs_h1437 <- list()
counter <- 1
for (x in cox_models_h1437) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
  hr_calcs_h1437[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_h1437 <- list()
counter <- 1
for (x in hr_calcs_h1437) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_h1437[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots_h1437 <- list()
counter <- 1
for (x in hr_calcs_h1437) {
  current_calc <- x
  current_pvalue <- hr_pvalues_h1437[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES H1437 Cell-line")
  hr_plots_h1437[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}





#CC cell-line all of the cell-lines together----
rownames(cc_cell_line_fpkm) <- cc_cell_line_fpkm$X
cc_cell_line_fpkm <- gene_name_cleaner(data.to.clean = cc_cell_line_fpkm)
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_line_fpkm<- subset(cc_cell_line_fpkm, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cc_cell_line_fpkm <- as.data.frame(cc_cell_line_fpkm)

#Now denoising the sc-data----
cc_cell_line_fpkm <- apply(cc_cell_line_fpkm, c(1,2), as.numeric)
cc_cell_line_fpkm <- magic_denoiser(sc.data = cc_cell_line_fpkm, magic.seed = 123)


#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP"), cell.data = cc_cell_line_fpkm$denoised_sc_dataframe, cell.meta = cc_cell_line_fpkm$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(cc_cell_line_fpkm$denoised_sc_dataframe)
save(mad.genes, file = "cc_cell_line_fpkm_all_cell_lines_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(cc_cell_line_fpkm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "cc_cell_line_fpkm_all_cell_lines_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "leukemia", cancer.type2 = "lung cancer", cancer.type3 = "colon cancer", mirna.remove = c("hsa-miR-129-2-3p","hsa-miR-129-1-3p"))
save(mirna.genes, file = "cc_cell_line_fpkm_all_cell_lines_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_cc_cell_line_fpkm_all_cell_lines.RData")

#Cox model----
cox_models_cc_all_cell_lines <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = df_for_train_test_split, gene.num = 900, cox.predictors = current_weight) 
  cox_models_cc_all_cell_lines[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#KM curves----
hr_calcs_cc_all_cell_lines <- list()
counter <- 1
for (x in cox_models_cc_all_cell_lines) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
  hr_calcs_cc_all_cell_lines[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_cc_all_cell_lines <- list()
counter <- 1
for (x in hr_calcs_cc_all_cell_lines) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_cc_all_cell_lines[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots_cc_all_cell_lines <- list()
counter <- 1
for (x in hr_calcs_cc_all_cell_lines) {
  current_calc <- x
  current_pvalue <- hr_pvalues_cc_all_cell_lines[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES All Cell-lines")
  hr_plots_cc_all_cell_lines[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}

