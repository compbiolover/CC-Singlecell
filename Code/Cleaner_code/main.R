#Name: main.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: The main script for my project.

#Loading single-cell data----
cc_tumor_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_tumor_all_cells_FPKM.csv")
k562_cells <- read.csv("Data/Single-cell-data/Other-cancers/GSE65525_RAW/GSM1599500_K562_cells.csv")
glio <- read.csv("Data/Single-cell-data/Other-cancers/GSE57872_GBM_data_matrix.txt", sep='\t')
cc_cell_line_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_Cell_Line_FPKM.csv")

#Loading additional bulk colon cancer datasets----
sw620_sw480 <- read.csv("Data/Datasets-for-comparison/GSE112568_fpkm_table_CRC_C_Williams_Addendum.txt", sep = '\t')

#Loading the types of cancer we can pick from for mirna metric----
high_choices <- read.csv("Data/Data-from-Cleaner-code/high-dbDEMC-cancer-types.csv")
low_choices <- read.csv("Data/Data-from-Cleaner-code/low-dbDEMC-cancer-types.csv")


#Preporcessing the various bulk datasets for comparision----


#Pre-processing for the different sc data files----

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

#K562 gastric cell line----
rownames(k562_cells) <- k562_cells$X
k562_cells <- subset(k562_cells, select=c(X0:X0.238))
colnames(k562_cells) <- seq(from=1,to=239,by=1)

#Glioblastoma patients----
rownames(glio) <- glio$X
glio <- subset(glio, select=c(MGH264_A01:MGH31Tumor))
glio <- abs(glio)

#CC tumor patients dataset analysis-----
#Now denoising the sc-data----
cc_tumor_fpkm <- magic_denoiser(sc.data = cc_tumor_fpkm,magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP"), cell.data = cc_tumor_fpkm$denoised_sc_dataframe, cell.meta = cc_tumor_fpkm$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
save(mad.genes, file = "cc_tumor_fpkm_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(cc_tumor_fpkm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "cc_tumor_fpkm_sde.RData")

#Mirna metric----
#mirna.genes <- mirna_calculator(cancer.type1 = "colon cancer", cancer.type2 = "colorectal cancer", ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, mirna.filename = "TargetScan_cc_tumor_patients.RData")
#save(mirna.genes, file = "Data/Data-from-Cleaner-code/cc_tumor_fpkm_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized.RData")

#Load the merged_df data frame for use with Cox model----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

calculated_days <- merged_df$days.to.death - merged_df$days.to.last.follow.up
calculated_days[calculated_days==0]=1
merged_df$days.to.last.follow.up <- ifelse(is.na(calculated_days), merged_df$days.to.last.follow.up, calculated_days)
merged_df$days.to.last.follow.up <- ifelse(merged_df$days.to.last.follow.up==0, 1, merged_df$days.to.last.follow.up)
cox_time <- merged_df$days.to.last.follow.up
cox_event <- merged_df$vital.status
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event

#Cox model----
cox_models <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#KM curves----
hr_calcs <- list()
counter <- 1
for (x in cox_models) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
  hr_calcs[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues <- list()
counter <- 1
for (x in hr_calcs) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots <- list()
counter <- 1
for (x in hr_calcs) {
  current_calc <- x
  current_pvalue <- hr_pvalues[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES")
  hr_plots[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}





#K562 cell-line analysis----
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
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "gastric cancer", mirna.filename = "k562_cells_mirnas_dbdemc_up.RData")
save(mirna.genes, file = "k562_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_k562.RData")

#Leukimia TCGA bulk data----
leuk_query <- GDCquery(project       = "TCGA-LAML",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = leuk_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-leuk-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
Leuk_data_se <- GDCprepare(leuk_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-leuk-Dataset/")
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
cox_time <- merged_df$days_to_last_follow_up
cox_event <- merged_df$vital_status
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$days.to.last.follow.up <- as.numeric(cox_df$days.to.last.follow.up)
cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up==0,1,cox_df$days.to.last.follow.up)



#Cox model----
cox_models_k562 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight) 
  cox_models_k562[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#KM curves----
hr_calcs_k562 <- list()
counter <- 1
for (x in cox_models_k562) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df)
  hr_calcs_k562[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_k562 <- list()
counter <- 1
for (x in hr_calcs_k562) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_k562[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots_k562 <- list()
counter <- 1
for (x in hr_calcs_k562) {
  current_calc <- x
  current_pvalue <- hr_pvalues_k562[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES K562")
  hr_plots_k562[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}



#Glio analysis----
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
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
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
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight) 
  cox_models_glio[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}



#KM curves----
merged_df <- merged_df %>% filter(!is.na(vital_status))
colnames(merged_df)[colnames(merged_df) == "vital_status"] <- "vital.status"
colnames(merged_df)[colnames(merged_df) == "days_to_last_follow_up"] <- "days.to.last.follow.up"

hr_calcs_glio <- list()
counter <- 1
for (x in cox_models_glio) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
  hr_calcs_glio[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_glio <- list()
counter <- 1
for (x in hr_calcs_glio) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_glio[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

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





#A549 cell-line analysis----
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
            directory       = "Data/Bulk-data/TCGA-lungaud-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
lung_data_se <- GDCprepare(lung_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-lungaud-Dataset/")
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
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event

#Cox model----
cox_models_a549 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors = current_weight) 
  cox_models_a549[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#KM curves----
hr_calcs_a549 <- list()
counter <- 1
for (x in cox_models_a549) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = cox_df)
  hr_calcs_a549[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_a549 <- list()
counter <- 1
for (x in hr_calcs_a549) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_a549[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots_a549 <- list()
counter <- 1
for (x in hr_calcs_a549) {
  current_calc <- x
  current_pvalue <- hr_pvalues_a549[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES A549 Cell-line with Lung AUD Bulk Dataset")
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


#Now denoising the sc-data----
hct116_cells <- apply(hct116_cells, c(1,2), as.numeric)
hct116_cells <- magic_denoiser(sc.data = hct116_cells, magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = "VIM", cell.data = hct116_cells$denoised_sc_dataframe, cell.meta = hct116_cells$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(hct116_cells$denoised_sc_dataframe)
save(mad.genes, file = "hct116_cells_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(hct116_cells$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "hct116_cells_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", max.miR.targets = 10, cancer.up = TRUE, cancer.type1 = "colon cancer", cancer.type2 = "colorectal cancer", mirna.filename = "TargetScane_hct116_cell_mirna.RData", print.ts.targets = TRUE, mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"))
save(mirna.genes, file = "hct116_cells_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_hct116_cells.RData")

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#Cox model----
cox_models_hct116 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800, cox.predictors =mirna.genes) 
  cox_models_hct116[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#KM curves----
hr_calcs_hct116 <- list()
counter <- 1
for (x in cox_models_hct116) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
  hr_calcs_hct116[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_hct116 <- list()
counter <- 1
for (x in hr_calcs_hct116) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_hct116[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots_hct116 <- list()
counter <- 1
for (x in hr_calcs_hct116) {
  current_calc <- x
  current_pvalue <- hr_pvalues_hct116[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA HCT116 Cell-line")
  hr_plots_hct116[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}




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

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/df_for_train_test_split.RData")

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

