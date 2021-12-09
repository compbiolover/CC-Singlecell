#Name: coad_dataset.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: All of the analysis for TCGA-COAD analysis

#Load needed packages----
library(akima)
library(arules)
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library(TCGAbiolinks)
library(tidyverse)

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
cc_tumor_fpkm <- subset(cc_tumor_fpkm, select=c(RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55))
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
mirna_targets <- seq(374,1550, by=2)
#mirnas <- seq(374,1550, by=2)
mirnas <- seq(1, 1550, by = 50)

for(x in mirnas[31:length(mirnas)]){
  mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                  cancer.type2 = "colon cancer",
                                  ts.org = "Human",
                                  ts.version = "7.2",
                                  max.miR.targets = x,
                                  cancer.up = TRUE,
                                  mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                  max.mirnas = x,
                                  write.heatmap.data = FALSE, 
                                  print.ts.targets = TRUE,
                                  save.mirna.raw.targets = FALSE,
                                  save.mirna.genes = TRUE,
                                  mirna.gene.filename = paste0("Data/TCGA-COAD/mirna_genes_global_search_", x,"_", x, "_targets.csv"),
                                  mirna.gene.rfile = paste0("Data/TCGA-COAD/mirna_genes_global_search_", x,"_", x, "_targets.RData"))
  
}






mirna.genes <- mirna_calculator(cancer.type1 = "colorectal cancer",
                                cancer.type2 = "colon cancer",
                                ts.org = "Human",
                                ts.version = "7.2",
                                max.miR.targets = 200,
                                cancer.up = TRUE,
                                mirna.remove = c("hsa-miR-129-2-3p", "hsa-miR-129-1-3p", "hsa-miR-454-3p", "hsa-miR-365a-5p"),
                                max.mirnas = 1558,
                                write.heatmap.data = FALSE, 
                                print.ts.targets = TRUE,
                                save.mirna.raw.targets = FALSE,
                                save.mirna.genes = TRUE,
                                mirna.gene.filename = "Data/TCGA-COAD/mirna_genes_1558_200_targets.csv",
                                mirna.gene.rfile = "Data/TCGA-COAD/mirna_genes_1558_200_targets.RData")


#Optimizing the weights of the three metric linear model----
mad_sdes_mirna_optimized <- three_weight_optimizer(first.metric = mad.genes,
                                                   second.metric = mirna.genes,
                                                   third.metric = sde.genes,
                                                   my.filename = "Data/TCGA-COAD/three_weight_optimized_40_mirna.RData")


#Loading the merged data frame----
#load("Data/Exported-data/R-objects/coad_df.RData")
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
load("Data/TCGA-COAD/Mirna-metric/three_weight_optimized_800_mirna.RData", verbose = TRUE)
mad_sdes_mirna_optimized <- integrated_gene_lists

cox_models <- list()
my_cindicies <- c()
counter <- 1

for (x in mad_sdes_mirna_optimized[34]) {
  current_weight <- x
  current_cox <- cox_model_fitter(my.seed = 1, cox.df = cox_df, gene.num = 1800,
                                  cox.predictors = current_weight,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = TRUE,
                                  save.regular.cox.genes = TRUE,
                                  my.filename = "~/Desktop/test.csv") 
  
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

#KM risk calculation----
patient_risk <- risk_score_calculator(my.file = "Data/TCGA-COAD/three_weight_optimized_output_1558_mirna_200_top.csv",
                                      my.title = "CC Singlecell MMS COAD 1558 MiRNA 100 Targets Only",
                                      tumor.data = FALSE, n.data = FALSE, cox.df = cox_df, show.pval = TRUE, show.pval.method = FALSE)
patient_risk




#For plotting the coefficients----
coef_df <- data.frame(coefs=cox_models$`1`$`Active Coefficients`, labels=cox_models$`1`$`Active Genes`)
coef_df_sub <- filter(coef_df, abs(coefs)>0.1)


coef_df_sub$p.value.trans <- -log10(coef_df_sub$p.value)
colnames(coef_df_sub)[2] <- "Gene"
coef_plot <- ggplot(data = coef_df_sub, aes(x=Gene, y=coefs, color=Gene, fill=Gene))+
  geom_col()+
  theme_bw()+
  ggtitle("CC Singlecell MMS COAD Cox Coefficients")+
  ylab("Coefficients")+
  xlab("Gene")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_flip()


coef_plot

#Finding the combination with the max value for mirna for this dataset
coad_mirna_df <- read.csv(file = "Documentation/mirna_data.csv")
coad_mirna_df <- filter(coad_mirna_df, Dataset=="TCGA-COAD")
coad_mirna_df <- filter(coad_mirna_df, MiRNA == 20 | MiRNA == 40 | MiRNA == 80 | MiRNA == 200 | MiRNA == 400 | MiRNA == 800 | MiRNA == 1000 | MiRNA == 1200 | MiRNA == 1500 | MiRNA == 1558)
coad_mirna_df <- filter(coad_mirna_df, Concordance_index==max(coad_mirna_df$Concordance_index))

#Making heatmap of optimal mirna value
# ggplot(coad_mirna_df, aes(MiRNA, MiRNA_targets)) +
#   geom_tile(aes(fill = Concordance_index), colour = "white") +
#   scale_fill_gradient(low = "white", high = "red")+
#   theme_minimal()+
#   ylab("MiRNA Targets")+
#   labs(fill="Concordance Index")+
#   ggtitle("TCGA-COAD")+
#   coord_fixed()



resolution <- 0.1
a <- interp(x=coad_mirna_df$MiRNA, y=coad_mirna_df$MiRNA_targets, z=coad_mirna_df$Concordance_index, 
            xo=seq(min(coad_mirna_df$MiRNA),max(coad_mirna_df$MiRNA_targets),by=resolution),
            yo=seq(min(coad_mirna_df$MiRNA_targets),max(coad_mirna_df$MiRNA_targets),by=resolution), duplicate="mean")

x.coordinate <- a$x
y.coordinate <- a$y
my.table <- as.data.frame(a$z) 

names(my.table) <- y.coordinate
my.table[is.na(my.table)] <- 0
my.table$x <- x.coordinate
longenergy <- melt(as.data.table(my.table),id="x")


ggheat <- ggplot(longenergy,aes(x=variable,y=x, fill=value)) +
  geom_tile() + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + 
  xlab("# of MiRNA Targets") + ylab("# of MiRNAs") 


ggheat


#Regular cox with the top coefficients from the penalized cox model
#Code to make the construction of the predictors efficient 
reg_cox <- function(cox.predictors = NULL, my_data=cox_df)
if(is.numeric(cox.predictors)==TRUE){
  my_predictors <- names(cox.predictors)
}else{
  my_predictors <- rownames(cox.predictors)
}

if(is.character(cox.predictors)==TRUE){
  my_predictors <- cox.predictors
}


my_predictors <- head(my_predictors, n=gene.num)
my_predictors <- sapply(my_predictors, gsub, pattern="-",replacement=".")
my_predictors <- unlist(my_predictors)
colname_changes <- sapply(colnames(cox.df), gsub, pattern="-",replacement=".")
colname_changes <- sapply(colnames(cox.df), gsub, pattern="_",replacement=".")
colname_changes <- sapply(colnames(cox.df), gsub, pattern="/",replacement=".")
colname_changes <- unlist(colname_changes)
colnames(cox.df) <- colname_changes
my_predictors <- intersect(my_predictors, colnames(cox.df))
my_predictors <- paste("~", paste(my_predictors[1:length(my_predictors)], collapse = "+"))

res.cox <- coxph(Surv(days.to.last.follow.up, vital.status) ~ SPANXN1 + RFPL3 + PAGE2B + MTRNR2L13 + FAM230A + CST9 + CMC4 + C2orf73 + ARMS2, data =  cox_df)
summary(res.cox)



#For plotting c-index across gene size----
coad_gene_size_df <- read.csv("Data/TCGA-COAD/cc_singlecell_mms_coad_gene_size_data.csv")
coad_gene_size_df$dataset <- "TCGA-COAD"
coad_read_gene_size_df <- read.csv("Data/TCGA-COAD_and_TCGA-READ/cc_singlecell_mms_gene_size_data.csv")
coad_read_gene_size_df$dataset <- "TCGA-COAD + TCGA-READ"
read_gene_size_df <- read.csv("Data/TCGA-READ/cc_singlecell_mms_read_gene_size_data.csv")
read_gene_size_df$dataset <- "TCGA-READ"
coad_size_plot <- ggplot(data = gene_size_df, aes(x=gene_num, y=concordance_index))+
  geom_line()+
  geom_point()+
  labs(x = "Gene Number",
       y = "Concordance Index",
       title = "TCGA-COAD Gene Size vs. Concordance Index")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

coad_size_plot

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

combined_size_df <- rbind(coad_gene_size_df, coad_read_gene_size_df, read_gene_size_df)
combined_gene_size_plot <- ggplot(data = combined_size_df, aes(x=gene_num, y=concordance_index, color=dataset))+
  geom_line()+
  geom_point()+
  labs(x = "Gene Number",
       y = "Concordance Index",
       title = "Gene Size vs. Concordance Index",
       color = "Dataset")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"))

combined_gene_size_plot




#MiRNA plot----
mirna_data <- read.csv("Documentation/mirna_data.csv")
mirna_data <- filter(mirna_data, Dataset=="TCGA-COAD")
#mirna_data <- filter(mirna_data, Dataset=="TCGA-COAD" | Dataset=="TCGA-READ" | Dataset== "TCGA-COAD + TCGA-READ")
mirna_plot <- ggplot(data = mirna_data, aes(x=MiRNA, y=Concordance_index, color=Dataset))+
  geom_line()+
  geom_point()+
  labs(x = "MiRNA Number",
       y = "Concordance Index",
       color = "Dataset",
       title = "TCGA-COAD")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

mirna_plot

#READ mirna plot
mirna_data <- read.csv("Documentation/mirna_data.csv")
mirna_data <- filter(mirna_data, Dataset=="TCGA-READ")
#mirna_data <- filter(mirna_data, Dataset=="TCGA-COAD" | Dataset=="TCGA-READ" | Dataset== "TCGA-COAD + TCGA-READ")
read_mirna_plot <- ggplot(data = mirna_data, aes(x=MiRNA, y=Concordance_index, color=Dataset))+
  geom_line(color="blue")+
  geom_point(color="blue")+
  labs(x = "MiRNA Number",
       y = "Concordance Index",
       color = "Dataset",
       title = "TCGA-READ")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

read_mirna_plot


#TCGA-COAD + TCGA-READ
mirna_data <- read.csv("Documentation/mirna_data.csv")
mirna_data <- filter(mirna_data, Dataset=="TCGA-COAD + TCGA-READ")

coad_read_mirna_plot <- ggplot(data = mirna_data, aes(x=MiRNA, y=Concordance_index, color=Dataset))+
  geom_line(color="green")+
  geom_point(color="green")+
  labs(x = "MiRNA Number",
       y = "Concordance Index",
       color = "Dataset",
       title = "TCGA-COAD + TCGA-READ")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

coad_read_mirna_plot


#MiRNA target graphs
mirna_target_data <- read.csv("Documentation/mirna_data.csv")
mirna_target_data <- filter(mirna_target_data, Dataset=="TCGA-COAD" & Colorect_only_mirnas.=="N")


coad_mirna_target_plot <- ggplot(data = mirna_target_data, aes(x=MiRNA_targets, y=Concordance_index, fill=MiRNA))+
  geom_col()+
  # geom_bar(stat = "identity")+
  # geom_line(color="blue")+
  # geom_point(color="blue")+
  labs(x = "MiRNA Target Number",
       y = "Concordance Index",
       color = "Dataset",
       title = "TCGA-COAD")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"))

coad_mirna_target_plot

ggplot(data = mirna_target_data, aes(x=MiRNA_targets, y=MiRNA, fill=Concordance_index))+
  geom_tile()

#Just the 0.1 weighted coefficients----
top_cox <- coxph(data = cox_df, formula = Surv(days.to.last.follow.up, vital.status)~SPANXN1+RFPL3+PAGE2B+MTRNR2L13+FAM230A+CST9+CMC4+C2orf73+ARMS2)
summary(top_cox)

#Stacked bar graph of top coefficient c-index vs total c-index----
stacked_data <- read.csv("Data/TCGA-COAD/top-predictor-cindex-vs-total-cindex.csv")
stacked_data_des <- discretize(stacked_data$Concordance_index_top_only)
ggplot(data = stacked_data, aes(x=Dataset, y=Concordance_index_all, fill=Concordance_index_top_only))+
  geom_bar(stat = "identity", position = "fill")+
  labs(x="Dataset",
       y="Concordance Index (total)",
       title = "Top Predictors vs. All Predictors")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))

#Heatmap of the top coefficient c-index vs total c-index----
ggplot(data = stacked_data, aes(x=Dataset, y=Concordance_index_all, fill=Concordance_index_top_only))+
  geom_tile()+
  labs(y="Concordance Index (all predictors)",
       title = "Top Predictors vs. All Predictors",
       fill = "Concordance Index (top predictors)")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")


