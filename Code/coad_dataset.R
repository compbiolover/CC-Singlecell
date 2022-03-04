#Name: coad_dataset.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: All of the analysis for TCGA-COAD analysis

#Load needed packages----
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library(TCGAbiolinks)
library(tidyverse)
library(viridis)

#Loading single-cell data----
cc_tumor_fpkm <- read.csv("Data/Single-cell-data/FPKM/GSE81861_CRC_tumor_all_cells_FPKM.csv")

#Processing the single-cell data to clean up the gene names----
rownames(cc_tumor_fpkm) <- cc_tumor_fpkm$X
current_colname_split <- strsplit(colnames(t(cc_tumor_fpkm)), "_")
finished_gene_list <- c()
for (x in seq(1:length(current_colname_split))){
  finished_gene_list <- c(finished_gene_list, current_colname_split[[x]][2])
}

cc_tumor_fpkm <- t(cc_tumor_fpkm)
colnames(cc_tumor_fpkm) <- finished_gene_list
cc_tumor_fpkm <- t(cc_tumor_fpkm)
cc_tumor_fpkm <- subset(cc_tumor_fpkm,
                        select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
cc_tumor_fpkm <- apply(cc_tumor_fpkm, c(1,2), as.numeric)
cc_tumor_fpkm <- as.matrix(cc_tumor_fpkm)


#Denoising the single-cell data----
cc_tumor_fpkm <- magic_denoiser(sc.data = cc_tumor_fpkm, magic.seed = 123)
save(cc_tumor_fpkm, file = "Data/TCGA-COAD/denoised-single-cell-data.RData")

#Now getting pseudotime info from Moncocle3----
cds_output <- cell_dataset_builder(vim.genes = c("VIM", "VIMP", "CDH1", "CDH2", "SNAI1",
                                                 "SNAI2", "TWIST1", "ZEB1", "GRHL2", "CCND1",
                                                 "CDK4", "CDK6"),
                                   cell.data = cc_tumor_fpkm$denoised_sc_dataframe,
                                   cell.meta = cc_tumor_fpkm$cds_gene_names)


#MAD metric----
mad.genes <- mad_calculator(cc_tumor_fpkm$denoised_sc_dataframe)
save(mad.genes, file = "Data/TCGA-COAD/mad.RData")

#Switchde metric----
sde.genes <- switchde_calculator(cc_tumor_fpkm$denoised_sc_dataframe, pseudo.time = cds_output$Pseudotime)
save(sde.genes, file = "Data/TCGA-COAD/sde.RData")

#MiRNA metric----
#For loop for testing several different mirna number or mirna target size
#mirnas <- seq(1,1557, by=2)
#For cancer upregulated miRNAS
#mirnas <- seq(1, 1560, by = 50)
#For cancer downregulaed miRNas
#mirnas <- seq(1, 82, by = 2)

for(x in mirnas[1:length(mirnas)]){
  mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                  cancer.type2 = "colon cancer",
                                  max.miR.targets = x,
                                  cancer.up = TRUE,
                                  mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                  max.mirnas = x,
                                  ts.org = "Human",
                                  ts.version = "7.2",
                                  print.ts.targets = TRUE,
                                  save.mirna.genes = TRUE,
                                  mirna.gene.filename = paste0("Data/TCGA-COAD/Mirna/Global_Search/mirna_genes_global_search_mirna_", x,"_", x, "_targets.csv"),
                                  mirna.gene.rfile = paste0("Data/TCGA-COAD/Mirna/Global_Search/mirna_genes_global_search_mirna_", x,"_", x, "_targets.RData"))
  
}





#Looping through several global search terms now that we know it works well
#Global target search
mirna_targets <- seq(10, 1010, by=100)
mirna_num <- seq(10, 100, 10)

#Local target around max
# mirna_targets <- seq(890,930, by = 5)
# mirna_targets <- mirna_targets[-5]
mirna_num <- 500

#Local target around second highest value to max
# mirna_targets <- seq(90, 130, by = 5)
# mirna_targets <- mirna_targets[-5]
# mirna_num <- 400

for(mt in mirna_targets[1]){
  mirna.genes <- mirna_calculator(cancer.type1 = "colon cancer",
                                  cancer.type2 = "colorectal cancer",
                                  max.miR.targets = mt,
                                  cancer.up = TRUE,
                                  mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                  max.mirnas = mirna_num,
                                  save.mirna.genes = TRUE,
                                  mirna.gene.rfile = paste0("~/Desktop/",mirna_num,"_",mt,"_targets.rds"))
  
  
  
  
}

#For individual miRNAs
mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                cancer.type2 = "colon cancer",
                                max.miR.targets = 110,
                                ts.version = "7.2",
                                cancer.up = TRUE,
                                mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p", "hsa-miR-873-5p", "hsa-miR-505-3p", "hsa-miR-496", "hsa-miR-483-3p", "hsa-miR-455-3p", "hsa-miR-411-5p"),
                                max.mirnas = 100,
                                save.mirna.genes = TRUE,
                                mirna.ranking.name = "~/Desktop/100-mirna-110-targets-mirna-ranking-test.rds")


#Optimizing the weights of the three metric linear model----
mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
                                                   second.metric = mirna.ranking,
                                                   third.metric = sde.genes,
                                                   my.filename = "~/Desktop/Optimization-test.RData")


#Loading the merged data frame----
coad_query <- GDCquery(project       = "TCGA-COAD",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the files
GDCdownload(query           = coad_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-COAD-Test-Dataset")


#Making the SummarizedExperiment object
coad_data_se <- GDCprepare(coad_query, summarizedExperiment = TRUE,
                           directory = "Data/Bulk-data/TCGA-COAD-Test-Dataset/")
coad_data_df <- as.data.frame(colData(coad_data_se))
coad_data_df$vital_status <- factor(coad_data_df$vital_status,
                                    levels = c("Alive", "Dead"),
                                    labels = c(0,1))
coad_data_df$vital_status <- as.numeric(as.character(coad_data_df$vital_status))


bulk_rna_df <- coad_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- coad_data_se@colData@rownames
rownames(bulk_rna_df) <- coad_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df,
                             select = unique(colnames(bulk_rna_df)))
coad_data_df_unique <- subset(coad_data_df,
                              select = unique(colnames(coad_data_df)))
merged_df <- merge(bulk_rna_df_unique, coad_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]



merged_df$days_to_last_follow_up <- ifelse(merged_df$vital_status==1,
                                           merged_df$days_to_death, 
                                           merged_df$days_to_last_follow_up)

merged_df <- filter(merged_df, days_to_last_follow_up != "NA")


cox_time <- merged_df$days_to_last_follow_up
cox_event <- merged_df$vital_status
cox_tumor <- merged_df$ajcc_pathologic_stage
cox_tumor_n <- merged_df$ajcc_pathologic_n
cox_tumor_m <- merged_df$ajcc_pathologic_m
cox_gender <- merged_df$gender
cox_eth <- merged_df$ethnicity
cox_race <- merged_df$race
cox_type <- merged_df$definition
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.m <- cox_tumor_m
cox_df$ajcc.n <- cox_tumor_n
cox_df$race <- cox_race
cox_df$ethnicity <- cox_eth
cox_df$gender <- cox_gender
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iv", replacement = 4)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iii", replacement = 3)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage ii", replacement = 2)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage i", replacement = 1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
cox_df$sample.type <- cox_type
cox_df <- filter(cox_df, !tumor.stage=="not reported")
cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]
cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up < 1, 1,
                 cox_df$days.to.last.follow.up)
saveRDS(cox_df, "Data/TCGA-COAD/coad_df_finished_v2.rds")
#save(cox_df, file = "Data/TCGA-COAD/coad_df_fpkm.RData")

# load("Data/Exported-data/R-objects/coad_df.RData")
# calculated_days <- merged_df$days.to.death - merged_df$days.to.last.follow.up
# calculated_days[calculated_days==0]=1
# merged_df$days.to.last.follow.up <- ifelse(is.na(calculated_days), merged_df$days.to.last.follow.up, calculated_days)
# merged_df$days.to.last.follow.up <- ifelse(merged_df$days.to.last.follow.up==0, 1, merged_df$days.to.last.follow.up)
# cox_time <- merged_df$days.to.last.follow.up
# cox_event <- merged_df$vital.status
# cox_tumor <- merged_df$tumor.stage
# cox_tumor_m <- merged_df$ajcc.pathologic.m
# cox_tumor_n <- merged_df$ajcc.pathologic.n
# cox_gender <- merged_df$gender
# cox_eth <- merged_df$ethnicity
# cox_race <- merged_df$race
# cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.3))
# cox_df$days.to.last.follow.up <- cox_time
# cox_df$vital.status <- cox_event
# cox_df$tumor.stage <- cox_tumor
# cox_df$ajcc.m <- cox_tumor_m
# cox_df$ajcc.n <- cox_tumor_n
# cox_df$race <- cox_race
# cox_df$ethnicity <- cox_eth
# cox_df$gender <- cox_gender
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="a", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="b", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="c", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iv", replacement = 4)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge iii", replacement = 3)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge ii", replacement = 2)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "stge i", replacement = 1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
# cox_df <- filter(cox_df, !tumor.stage=="not reported")
# cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]
#save(cox_df, file = "Data/TCGA-COAD/coad_df.RData")
load("Data/TCGA-COAD/coad_df.RData", verbose = TRUE)


#Cox Model Weight Optimization----
load("~/Desktop/Optimization_500_810_targets.RData", verbose = TRUE)
mad_sdes_mirna_optimized <- integrated_gene_lists

cox_models <- vector(mode = "list", length = 12)
my_cindicies <- vector(mode = "numeric", length = 1)
my_active_coefs <- vector(mode = "character", length = 1)
counter <- 1


  
  
for (x in mad_sdes_mirna_optimized[27]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1,
                                  cox.df = cox_df,
                                  gene.num = 1800,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/test-mean-active-genes.csv") 
  
  cox_models[[as.character(counter)]] <- current_cox
  counter <- counter + 1
  
  #Storing all of the c-index values in a vector that we can use later to build the plot
  c_finder <-current_cox$CV$index[1]
  current_c <- current_cox$CV$cvm[c_finder]
  current_c <- round(current_c, digits = 4)
  my_cindicies <- c(my_cindicies, current_c)
  current_coefs <- length(current_cox$`Active Coefficients`)
  my_active_coefs <- c(my_active_coefs, current_coefs)
}

# top_cindex <-max(my_cindicies)
top_index <- which.max(my_cindicies)
top_index

cox_models$`1`$CV





#Just individual cox model
coad_cox<- cox_model_fitter(my.seed = 1,
                                cox.df = cox_df,
                                gene.num = 1800,
                                cox.predictors = integrated_gene_lists[[27]],
                                tumor.stage = FALSE,
                                tumor.n = FALSE,
                                tumor.m = FALSE,
                                regular.cox = TRUE,
                                save.regular.cox.genes = TRUE,
                                my.filename = "~/Desktop/coad_top_genes.csv") 


#Generating random genes that are equal to our ideal CC Singlecell MMS amount----
random_genes <- sample(colnames(cox_df), size = 1100, replace = FALSE)

#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "~/Desktop/test-mean-active-genes.csv",
                                      my.title = "COAD",
                                      tumor.data = FALSE, n.data = FALSE, cox.df = cox_df, show.pval = TRUE, show.pval.method = FALSE)
patient_risk




#For plotting the coefficients----
coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
coef_df_sub <- filter(coef_df, abs(coefs)>0.5)


#coef_df_sub$p.value.trans <- -log10(coef_df_sub$p.value)
colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("READ Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"))+
  scale_color_viridis_d(direction = -1)+
  scale_fill_viridis_d(direction = -1)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot






#For plotting c-index across gene size----
coad_gene_size_df <- read.csv("Data/TCGA-COAD/cc_singlecell_mms_coad_gene_size_data.csv")
coad_gene_size_df$dataset <- "TCGA-COAD"
#coad_read_gene_size_df <- read.csv("Data/TCGA-COAD_and_TCGA-READ/cc_singlecell_mms_gene_size_data.csv")
#coad_read_gene_size_df$dataset <- "TCGA-COAD + TCGA-READ"
read_gene_size_df <- read.csv("Data/TCGA-READ/cc_singlecell_mms_read_gene_size_data.csv")
finished_size_df <- read.csv("Data/TCGA-COAD/finished_gene_size_data.csv")
read_gene_size_df$dataset <- "TCGA-READ"
coad_size_plot2 <- ggplot(data = df_to_plot, aes(x=Gene_num, y=C_index))+
  geom_line()+
  geom_point()+
  labs(x = "Gene Number",
       y = "Concordance Index",
       title = "Gene Size vs. Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

coad_size_plot2 

coad_read_size_plot <- ggplot(data = coad_read_gene_size_df, aes(x=gene_num, y=concordance_index))+
  geom_line()+
  geom_point()+
  labs(x = "Gene Number",
       y = "Concordance Index",
       title = "TCGA-COAD + TCGA-READ Gene Size vs. Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

coad_read_size_plot


read_size_plot <- ggplot(data = read_gene_size_df, aes(x=gene_num, y=concordance_index))+
  geom_line()+
  geom_point()+
  labs(x = "Gene Number",
       y = "Concordance Index",
       title = "TCGA-READ Gene Size vs. Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))


read_size_plot

combined_size_df <- rbind(coad_gene_size_df, read_gene_size_df)
combined_gene_size_plot <- ggplot(data = finished_size_df, aes(x=Gene_size, y=C_index, color=Dataset))+
  geom_line()+
  geom_point(size = 4.5)+
  labs(x = "Gene Number",
       y = "Mean Concordance Index",
       title = "Gene Size vs. Concordance Index",
       color = "Dataset")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom")+
  scale_color_viridis_d(direction = -1)

combined_gene_size_plot



#Gene size vs. active coefs
active_coefs_gene_size_plot <- ggplot(data = combined_active_coefs_df, aes(x=gene_sizes, y=active_coefs, color=Dataset))+
  geom_line()+
  geom_point(size = 4.5)+
  labs(x = "Gene Number",
       y = "Mean Active Coefficients",
       title = "Gene Size vs. Active Coefficients",
       color = "Dataset")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom")+
  scale_color_viridis_d(direction = -1)

active_coefs_gene_size_plot




#Stacked bar graph of top coefficient c-index vs total c-index----
stacked_data <- read.csv("Data/TCGA-COAD/top-predictor-cindex-vs-total-cindex.csv")
stacked_data <- filter(stacked_data, Dataset != "TCGA-COAD + TCGA-READ")
stacked_data$percent_covered <- c(93.6, 98.6)
#stacked_data_des <- discretize(stacked_data$Concordance_index_top_only)



#MiRNA size heat map----
mirna_heatmap_df <- read.csv("Data/TCGA-COAD/MiRNA/heatmap_data.csv")


#Reading all the files in to get the C-index
c_index_getter <- function(filename){
  finished_filename <- paste0("Data/TCGA-COAD/Mirna/Global_Search/Outputs/Top-cindices/",filename)
  my_file <- read.csv(finished_filename)
  current_c <- my_file$c_index[1]
  return(current_c)
}

files <- list.files("Data/TCGA-COAD/Mirna/Global_Search/Outputs/Top-cindices/")
all_cindicies <-sapply(files, c_index_getter)
all_cindicies <- unname(all_cindicies)

#Appending our c-indicies to the larger file for global search
mirna_heatmap_df$C_index <- all_cindicies


#Making the heat map
heatmap_coad <- ggplot(data = mirna_heatmap_df, aes(x=Mirna, y=Mirna_targets, fill=C_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(C_index, 4)), color = "white")+ 
  coord_fixed()+
  labs(x ="# of MiRNAs",
       y = "# of MiRNA Targets",
       title = "Old COAD Global Search",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "bottom",
        legend.key.width = unit(2.5,"cm"))

#Changing to color-blind friendly palette
heatmap_coad + scale_fill_viridis_c()

#Individual metrics----
individual_metrics_df <- read.csv(file = "Data/TCGA-COAD/indivdual_metrics_cindex_updated.csv")
individual_metrics_df$Method <- factor(individual_metrics_df$Method, levels = c("CC Singlecell MM", "CC Singlecell MS", "CC Singlecell MMS", "MAD", "MiRNA", "SDE", "Random Genes"))
individual_metrics_read_df <-read.csv("Data/TCGA-READ/individual_vs_cc_singlecell_mms_comparison_updated.csv")
individual_metrics_read_df$Method <- factor(individual_metrics_read_df$Method, levels = c("CC Singlecell MM", "CC Singlecell MS", "CC Singlecell MMS", "MAD", "MiRNA", "SDE", "Random Genes"))
individual_graph_coad <-ggplot(data = individual_metrics_df, aes(x=Method, y=C_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-COAD",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 40, family = "sans", angle = 90),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none",
        legend.key.width = unit(2.5, "cm"))

individual_graph_coad + coord_cartesian(ylim = c(0.5,0.71))+ scale_fill_manual(values = c("#00AFBB","#00AFBB", "#00AFBB", "#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07")) 

  


#Other methods comparison----
coad_comparison_df <- read.csv("Data/TCGA-COAD/cindices_vs_other_methods.csv")
coad_comparison_df <- coad_comparison_df[1:5,1:5]
colnames(coad_comparison_df) <- c("Num","Method", "Bulk Dataset", "Method 2", "C_index" )
comparison_graph_coad <-ggplot(data = coad_comparison_df, aes(x=Method, y=C_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-COAD",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none",
        legend.key.width = unit(2.5, "cm"))

comparison_graph_coad + coord_cartesian(ylim = c(0.5,0.71)) + scale_fill_viridis_d(direction = -1) 


#READ comparison
read_comparison_df <- read.csv("Data/TCGA-READ/methods-comparison-read-updated.csv")
colnames(read_comparison_df) <- c("Num","Method", "Bulk Dataset", "Method 2", "C_index" )
comparison_graph_read <-ggplot(data = read_comparison_df, aes(x=Method, y=C_index, fill=Method))+
  geom_bar(stat = "identity")+
  labs(title = "TCGA-READ",
       x = "Method",
       y = "Concordance Index",
       fill = "Concordance Index")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none",
        legend.key.width = unit(2.5, "cm"))

comparison_graph_read + coord_cartesian(ylim = c(0.5,0.85)) + scale_fill_viridis_d(direction = -1) 


#Combining the plots together
combined_comparison_graph <- ggarrange(comparison_graph_coad, comparison_graph_read, ncol = 1, nrow = 2, labels = c("A.", "B."))+ scale_fill_viridis_d(direction = -1)

#Active gene set size graph----
active_genes_plot <- ggplot(data = gene_sizes_df, aes(x=Method, y=Gene_num, fill=Method))+
  geom_bar(stat = "identity")+
  theme_bw()+
  ggtitle("COAD")+
  ylab("Active Gene Set Sizes")+
  xlab("Methods")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"))+
  scale_fill_viridis_d(direction = -1)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


active_genes_plot


