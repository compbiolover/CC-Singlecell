#lusc_dataset.R
#Loading needed packages----
library(ggplot2);packageVersion("ggplot2")
library(grid);packageVersion("grid")
library(SummarizedExperiment);packageVersion("SummarizedExperiment")
library(TCGAbiolinks);packageVersion("TCGAbiolinks")
library(tidyverse);packageVersion("tidyverse")

#Single-cell data----
lc_tumor_tpm <- read.csv("Data/Single-cell-data/Other-cancers/GSE69405_PROCESSED_GENE_TPM_ALL.txt", sep = '\t')
lc_tumor_tpm <- lc_tumor_tpm[!base::duplicated(lc_tumor_tpm$gene_name),]
rownames(lc_tumor_tpm) <- lc_tumor_tpm$gene_name

#Pre-processing
lc_tumor_tpm <- select(lc_tumor_tpm, contains("LC.PT.45" ) | contains("LC.MBT.15"))
lc_tumor_tpm <- select(lc_tumor_tpm, contains("LC.PT.45_SC") | contains( "LC.MBT.15_SC"))
lc_names <- c("VIM", "VIMP", "VIMP1")

#Denoising
lc_tumor_tpm <- magic_denoiser(sc.data = lc_tumor_tpm, magic.seed = 123)
cds_output <- cell_dataset_builder(vim.genes = lc_names, cell.data = lc_tumor_tpm$denoised_sc_dataframe, cell.meta = lc_tumor_tpm$cds_gene_names)

#MAD metric
mad.genes <- mad_calculator(lc_tumor_tpm$denoised_sc_dataframe)
save(mad.genes, file = "Data/TCGA-LUSC/mad.RData")

#SDES metric
sde.genes <- switchde_calculator(lc_tumor_tpm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "Data/TCGA-LUSC/sde.RData")

#MiRNA metric
mirnas <- seq(1, 810, by = 50)
#For constant miRNA number while varying the target number
for(y in mirnas){
  for(x in mirnas[1:length(mirnas)]){
    mirna.genes <- mirna_calculator(cancer.type1 = "lung cancer",
                                    max.miR.targets = x,
                                    cancer.up = TRUE,
                                    mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"),
                                    max.mirnas = x,
                                    save.mirna.genes = TRUE,
                                    mirna.gene.filename = paste0("Data/TCGA-LUSC/MiRNA/Global_Search/mirna_genes_global_search_mirna_fill_in_", y,"_", x, "_targets.csv"),
                                    mirna.gene.rfile = paste0("Data/TCGA-LUSC/Mirna/Global_Search/mirna_genes_global_search_mirna_fill_in_", y,"_", x, "_targets.RData"))
    
  }
}


#For changing miRNA and miRNA targets
for(x in mirnas[1:length(mirnas)]){
  mirna.genes <- mirna_calculator(cancer.type1 = "lung cancer",
                                  max.miR.targets = x,
                                  cancer.up = TRUE,
                                  mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"),
                                  max.mirnas = x,
                                  save.mirna.genes = TRUE,
                                  mirna.gene.filename = paste0("Data/TCGA-LUSC/MiRNA/Global_Search/mirna_genes_global_search_mirna_", x,"_", x, "_targets.csv"),
                                  mirna.gene.rfile = paste0("Data/TCGA-LUSC/MiRNA/Global_Search/mirna_genes_global_search_mirna_", x,"_", x, "_targets.RData"))
  
}



#Single miRNA
mirna.genes <- mirna_calculator(cancer.type1 = "lung cancer",
                                ts.org = "Human",
                                ts.version = "7.2",
                                max.miR.targets = 800,
                                max.mirnas = 800,
                                cancer.up = TRUE,
                                save.mirna.raw.targets = FALSE,
                                save.mirna.genes = TRUE,
                                print.ts.targets = TRUE,
                                write.heatmap.data = FALSE,
                                mirna.gene.filename = "Data/TCGA-LUSC/mirna_genes_800_800_targets.csv",
                                mirna.gene.rfile = "Data/TCGA-LUSC/mirna_genes_800_800_targets.RData",
                                mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p"))


#Lung cancer TCGA bulk data----
load("Data/TCGA-LUSC/lusc_df.RData", verbose = TRUE)

load("Data/TCGA-LUSC/three_weight_optimized_20_mirna_10.RData", verbose = TRUE)
mad_sdes_mirna_optimized <- integrated_gene_lists

cox_models <- list()
my_cindicies <- c()
counter <- 1

for (x in mad_sdes_mirna_optimized[113]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "Data/TCGA-LUSC/top_lusc.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
}


cox_models$`1`$CV



#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "Data/TCGA-LUSC/three_weight_optimized_output_20_mirna_10_top.csv",
                                      my.title = "CC Singlecell MMS LUSC",
                                      tumor.data = FALSE,
                                      n.data = FALSE,
                                      cox.df = cox_df,
                                      show.pval = TRUE,
                                      show.pval.method = FALSE)
patient_risk




#For plotting the coefficients----
coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
coef_df_sub <- filter(coef_df, abs(coefs)>0.0001)
colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("CC Singlecell MMS LUSC Cox Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot




