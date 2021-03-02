#Name: main.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: The main script for my project.

#Loading single-cell data----
cc_tumor_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_tumor_all_cells_FPKM.csv")
k562_cells <- read.csv("Data/Single-cell-data/Other-cancers/GSE65525_RAW/GSM1599500_K562_cells.csv")
glio <- read.csv("Data/Single-cell-data/Other-cancers/GSE57872_GBM_data_matrix.txt", sep='\t')
cc_cell_line_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_Cell_Line_FPKM.csv")

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
cds_output <- cell_dataset_builder(vim.genes = "VIM", cell.data = cc_tumor_fpkm$denoised_sc_dataframe, cell.meta = cc_tumor_fpkm$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
save(mad.genes, file = "cc_tumor_fpkm_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(cc_tumor_fpkm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "cc_tumor_fpkm_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(cancer.type = "colon cancer", ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10)
save(mirna.genes, file = "cc_tumor_fpkm_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized.RData")

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/df_for_train_test_split.RData")

#Cox model----
cox_models <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = df_for_train_test_split, gene.num = 900, cox.predictors = current_weight) 
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

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
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10)
save(mirna.genes, file = "k562_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_k562.RData")

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/df_for_train_test_split.RData")

#Cox model----
cox_models_k562 <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = df_for_train_test_split, gene.num = 900, cox.predictors = current_weight) 
  cox_models_k562[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#KM curves----
hr_calcs_k562 <- list()
counter <- 1
for (x in cox_models_k562) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
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
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10)
save(mirna.genes, file = "glio_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_glio.RData")

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/df_for_train_test_split.RData")

#Cox model----
cox_models_glio <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = df_for_train_test_split, gene.num = 900, cox.predictors = current_weight) 
  cox_models_glio[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#KM curves----
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




#CC cell-line analysis----
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
cc_cell_line_fpkm <- gene_name_cleaner(data.to.clean = cc_cell_line_fpkm)
cc_cell_line_fpkm <- t(cc_cell_line_fpkm)
cc_cell_line_fpkm<- subset(cc_cell_line_fpkm, select=c(RHA015__A549__turquoise:RHC2506__H1_B2__brown))
cc_cell_line_fpkm <- as.data.frame(cc_cell_line_fpkm)
h1_and_gm_cells <- select(cc_cell_line_fpkm, contains("H1") | contains("GM"))
h1_and_gm_cells <- select(h1_and_gm_cells, contains("RHG") | contains("RHC"))

#Now denoising the sc-data----
h1_and_gm_cells <- apply(h1_and_gm_cells, c(1,2), as.numeric)
h1_and_gm_cells <- magic_denoiser(sc.data = h1_and_gm_cells, magic.seed = 123)

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = "VIM", cell.data = h1_and_gm_cells$denoised_sc_dataframe, cell.meta = h1_and_gm_cells$cds_gene_names)

#MAD metric----
mad.genes <- mad_calculator(h1_and_gm_cells$denoised_sc_dataframe)
save(mad.genes, file = "h1_and_gm_cells_mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(h1_and_gm_cells$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "h1_and_gm_cells_sde.RData")

#Mirna metric----
mirna.genes <- mirna_calculator(ts.org = "Human", ts.version = "7.2", ts.num = 900, max.miR.targets = 10)
save(mirna.genes, file = "h1_and_gm_cells_mirna.RData")

#Optimizing the mirna + SDE metric----
mirna_sde_optimized <- two_weight_optimizer(first.metric = sde.genes, second.metric = mirna.genes, my.filename = "sde_mirna_optimized_h1_and_gm_cells.RData")

#Load the df_for_train_test_split data frame for use with Cox model----
load("Data/Exported-data/R-objects/df_for_train_test_split.RData")

#Cox model----
cox_models_h1_and_gm <- list()
counter <- 1
for (x in mirna_sde_optimized) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = df_for_train_test_split, gene.num = 900, cox.predictors = current_weight) 
  cox_models_h1_and_gm[[as.character(counter)]] <- current_cox
  counter <- counter + 1
}

#Loading data.frame needed for KM curves----
load("Data/Exported-data/R-objects/merged_df_replaced.RData")

#KM curves----
hr_calcs_h1_and_gm <- list()
counter <- 1
for (x in cox_models_h1_and_gm) {
  current_cox <- x
  current_hr <- hr_calculator(model.coefs = current_cox$Coefficients, data = merged_df)
  hr_calcs_h1_and_gm[[as.character(counter)]] <- current_hr
  counter <- counter + 1
}

#KM p-values----
hr_pvalues_h1_and_gm <- list()
counter <- 1
for (x in hr_calcs_h1_and_gm) {
  current_hr <- x
  current_pvalue <- km_pvalue_calculator(surv.time = current_hr$DF$days.to.last.follow.up, surv.status = current_hr$DF$vital.status, surv.predictors = current_hr$DF$risk, surv.df = current_hr$DF, num.sig.figs = 4, scientific.p.value = TRUE)
  hr_pvalues_h1_and_gm[[as.character(counter)]] <- current_pvalue
  counter <- counter + 1
}

#KM plots----
hr_plots_h1_and_gm <- list()
counter <- 1
for (x in hr_calcs_h1_and_gm) {
  current_calc <- x
  current_pvalue <- hr_pvalues_h1_and_gm[counter]
  current_plot <- km_plotter(km.fit = current_calc$KM, data.source = current_calc$DF, p.value = current_pvalue, plot.title = "MiRNA + SDES H1 and GM CC Cell-line")
  hr_plots_h1_and_gm[[as.character(counter)]] <- current_plot
  counter <- counter + 1
}
