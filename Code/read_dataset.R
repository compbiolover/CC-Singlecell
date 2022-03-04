#Name: read_dataset.R
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
            directory       = "Data/TCGA-READ/TCGA-READ-COUNTS-Dataset")


#Making the SummarizedExperiment object
read_data_se <- GDCprepare(read_query, summarizedExperiment = TRUE,
                           directory = "Data/TCGA-READ/TCGA-READ-COUNTS-Dataset/")
read_data_df <- as.data.frame(colData(read_data_se))
read_data_df$vital_status <- factor(read_data_df$vital_status,
                                    levels = c("Alive", "Dead"),
                                    labels = c(0,1))
read_data_df$vital_status <- as.numeric(as.character(read_data_df$vital_status))


bulk_rna_df <- read_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- read_data_se@colData@rownames
rownames(bulk_rna_df) <- read_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

bulk_rna_df_unique <- subset(bulk_rna_df,
                             select = unique(colnames(bulk_rna_df)))
read_data_df_unique <- subset(read_data_df,
                              select = unique(colnames(read_data_df)))
merged_df <- merge(bulk_rna_df_unique, read_data_df_unique, by = 'barcode')
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
saveRDS(cox_df, "Data/TCGA-READ/read_df_finished_v2.rds")
# read_data_se <- GDCprepare(read_query, summarizedExperiment = TRUE, directory = "Data/TCGA-READ/TCGA-READ-COUNTS-Dataset/")
# read_data_df <- as.data.frame(colData(read_data_se))
# read_data_df$vital_status <- factor(read_data_df$vital_status, levels = c("Alive", "Dead"), labels = c(0,1))
# read_data_df$vital_status <- as.numeric(as.character(read_data_df$vital_status))
# 
# bulk_rna_df <- read_data_se@assays@data@listData[["HTSeq - Counts"]]
# colnames(bulk_rna_df) <- read_data_se@colData@rownames
# rownames(bulk_rna_df) <- read_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
# bulk_rna_df <- t(bulk_rna_df)
# bulk_rna_df <- as.data.frame(bulk_rna_df)
# bulk_rownames <- rownames(bulk_rna_df)
# bulk_rna_df$barcode <- bulk_rownames
# 
# bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
# read_data_df_unique <- subset(read_data_df, select = unique(colnames(read_data_df)))
# merged_df <- merge(bulk_rna_df_unique, read_data_df_unique, by = 'barcode')
# rownames(merged_df) <- merged_df$barcode
# merged_df <- merged_df[,2:length(colnames(merged_df))]
# 
# 
# merged_df <- merged_df[complete.cases(merged_df[, "days_to_last_follow_up"]), ]
# merged_df$days_to_last_follow_up <- ifelse(merged_df$days_to_last_follow_up==0,1,merged_df$days_to_last_follow_up)
# cox_time <- merged_df$days_to_last_follow_up
# cox_event <- merged_df$vital_status
# cox_tumor <- merged_df$ajcc_pathologic_stage
# cox_tumor_n <- merged_df$ajcc_pathologic_n
# cox_tumor_m <- merged_df$ajcc_pathologic_m
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
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iv", replacement = 4)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iii", replacement = 3)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage ii", replacement = 2)
# cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage i", replacement = 1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
# cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
# cox_df <- filter(cox_df, !tumor.stage=="not reported")
# cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]
# #save(cox_df, file = "Data/TCGA-READ/read_df_fpkm.RData")

#Bulk data frame for Read merged data frame----
load("Data/TCGA-READ/read_df.RData", verbose = TRUE)

#MiRNA metric----
mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                max.miR.targets = 40,
                                cancer.up = TRUE,
                                mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                max.mirnas = 100,
                                save.mirna.genes = TRUE,
                                mirna.gene.filename = "Data/TCGA-READ/mirna_genes_40_100_targets.csv",
                                mirna.gene.rfile = "Data/TCGA-READ/mirna_genes_40_100_targets.RData")


#Different miRNA sizes
#mirnas <- seq(1,1074, by=2)
mirnas <- seq(1, 1074, by = 50)

#For constant miRNA number while varying the target number
for(y in mirnas){
  for(x in mirnas[1:length(mirnas)]){
    mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                    max.miR.targets = x,
                                    cancer.up = TRUE,
                                    mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                    max.mirnas = x,
                                    save.mirna.genes = TRUE,
                                    mirna.gene.filename = paste0("Data/TCGA-READ/MiRNA/Global_Search/mirna_genes_global_search_mirna_fill_in_", y,"_", x, "_targets.csv"),
                                    mirna.gene.rfile = paste0("Data/TCGA-READ/Mirna/Global_Search/mirna_genes_global_search_mirna_fill_in_", y,"_", x, "_targets.RData"))
    
  }
}

#For changing miRNA and miRNA targets
for(x in mirnas[1:length(mirnas)]){
  mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                  max.miR.targets = x,
                                  cancer.up = TRUE,
                                  mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                  max.mirnas = x,
                                  save.mirna.genes = TRUE,
                                  mirna.gene.filename = paste0("Data/TCGA-READ/MiRNA/Global_Search/mirna_genes_global_search_mirna_", x,"_", x, "_targets.csv"),
                                  mirna.gene.rfile = paste0("Data/TCGA-READ/MiRNA/Global_Search/mirna_genes_global_search_mirna_", x,"_", x, "_targets.RData"))
  
}




#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "Data/TCGA-READ/MiRNA/three_weight_optimized_output_1558_mirna_10_top.csv",
                                      my.title = "CC Singlecell MMS READ 1558 MiRNA 100 Targets Only",
                                      tumor.data = FALSE, n.data = FALSE, cox.df = cox_df, show.pval = TRUE, show.pval.method = FALSE)
patient_risk


load("Data/TCGA-READ/MiRNA/three_weight_optimized_40_mirna_100_colorect.RData", verbose = TRUE)
mad_sdes_mirna_optimized <- integrated_gene_lists

cox_models <- list()
my_cindicies <- c()
counter <- 1

for (x in mad_sdes_mirna_optimized[2]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 10,
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

#MiRNA size heat map----
mirna_heatmap_df <- read.csv("Data/TCGA-READ/MiRNA/Global_Search/mirna-gene-size-vs-gene-target-data-read.csv")

#Reading all the files in to get the C-index
c_index_getter <- function(filename){
  finished_filename <- paste0("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-READ/Mirna/Global_Search/Finished-outputs/",filename)
  my_file <- read.csv(finished_filename)
  current_c <- my_file$c_index[1]
  return(current_c)
}

files <- list.files("Data/TCGA-READ/Mirna/Global_Search/Finished-outputs/")
all_cindicies <-sapply(files, c_index_getter)
all_cindicies <- unname(all_cindicies)

#Appending our c-indicies to the larger file
mirna_heatmap_df$C_index <- all_cindicies

#Making the heat map
ggplot(data = mirna_heatmap_df, aes(x=Mirna, y=Mirna_targets, fill=C_index))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  coord_fixed()+
  labs(x ="# of MiRNAs",
       y = "# of MiRNA Targets",
       title = "READ Global Search",
       fill = "Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")



#Example of miRNA metric score----
example_metric_graph <- ggplot(data = example_mirna_metric_df, aes(x=Gene, y=score, fill=Gene))+
  geom_bar(stat = "identity")+
  labs(x = "Genes",
       y = "Score",
       title = "miRNA Metric Score")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 30, family = "sans"),
        axis.text.y = element_text(size = 40, family = "sans"),
        legend.text = element_text(size = 25, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "none")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_fill_viridis_d(direction = -1)

#Example miRNA heatmap scoring matrix
#Heatmap data file
example_mirna_heatmap <- read.csv("~/Desktop/example_mirna_metric_matrix.csv")
example_mirna_heatmap_sub <- example_mirna_heatmap[1:80,1:20]
colnames(example_mirna_heatmap_sub)[1] <- "Genes"
example_mirna_heatmap_sub_long <- gather(example_mirna_heatmap_sub, "miRNA", "score", 2:20)

example_heatmap_graph <- ggplot(data = example_mirna_heatmap_sub_long, aes(x=Genes, y = miRNA, fill=factor(score)))+
  geom_tile()+
  coord_fixed()+
  labs(x="Genes",
       y= "miRNAs",
       title = "miRNA Metric Interaction Map",
       fill = "Score")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40, family = "sans"),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 40, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 40, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 20, family = "sans"),
        axis.text.y = element_text(size = 30, family = "sans"),
        legend.text = element_text(size = 30, family = "sans"),
        legend.title = element_text(size = 40, family = "sans"),
        legend.position = "right")+
  scale_x_discrete(guide = guide_axis(angle = 45, check.overlap = TRUE))+
  scale_y_discrete(guide = guide_axis(angle = 45, check.overlap = TRUE))+
  scale_fill_viridis_d()
example_heatmap_graph

