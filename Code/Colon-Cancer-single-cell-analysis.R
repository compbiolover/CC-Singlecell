#Author: Andrew Willems <awillems@vols.utk.edu>.
#Purpose: To take scRNA-seq colon cancer data and begin analysis.
#R version: 4.0.2.

#Loading needed R packages. ----
library(caret);packageVersion("caret")
library(dplyr);packageVersion("dplyr")
library(glmnet);packageVersion("glmnet")
library(hoardeR);packageVersion("hoardeR")
library(Matrix);packageVersion("Matrix")
library(mlr);packageVersion("mlr")
library(mlr3);packageVersion("mlr3")
library(mlr3proba);packageVersion("mlr3proba")
library(monocle3);packageVersion("monocle3")
library(PTC);packageVersion("PTC")
library(pROC);packageVersion("pROC")
library(Rmagic);packageVersion("Rmagic")
library(SummarizedExperiment);packageVersion("SummarizedExperiment")
library(stringr);packageVersion("stringr")
library(survival);packageVersion("survival")
library(survminer);packageVersion("survminer")
library(switchde);packageVersion("switchde")
library(targetscan.Hs.eg.db);packageVersion("targetscan.Hs.eg.db")
library(TCGAbiolinks);packageVersion("TCGAbiolinks")
library(tidyverse);packageVersion("tidyverse")

#Setting the directory to where our data files are.----
setwd("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/")

#Copying in the functions that are helpful from the Li et al. paper. ----

# geneRank <- function(ranking1 = NULL,
#                     ranking2  = NULL,
#                     a1        = 1,
#                     a2        = 0){
#   gns = c(names(ranking1),names(ranking2))
#   gns = unique(gns)
#   res = rep(0, length(gns))
#   names(res) = gns
#   for(i in names(res)) {
#     res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2
#   }
#   #res=res/sum(res) 
#   res = res[order(res, decreasing = T)]
#   res
# }

# getRank <- function(ranking = NULL,
#                    gn      = NULL){
#   if (gn %in% names(ranking)) {
#     return(ranking[gn])
#   }
#   else return(0.0)
# }

# Ranking genes by three measurements
geneRank <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, a1 = 1,
                     a2 = 0, a3 = 0){
  gns = c(names(ranking1),names(ranking2),names(ranking3))
  gns = unique(gns)
  res = rep(0, length(gns))
  names(res) = gns
  for(i in names(res)) {
    res[i] = getRank(ranking1, i)*a1+getRank(ranking2, i)*a2+getRank(ranking3, i)*a3
  }
  #res=res/sum(res)
  #res = res[order(res, decreasing = T, method = "radix")]
  res
}

getRank <- function(ranking = NULL, gn = NULL){
  if (gn %in% names(ranking)) {
    return(ranking[gn])
  }
  else return(0.0)
}
# grid.search <- function(ranking1 = NULL,
#                        ranking2 = NULL,
#                        N        = 50){
#   res = NULL
#   for(a1 in seq(0,1)){
#     for(a2 in seq(0,1-a1)){
#       ranking = geneRank(ranking1, ranking2,a1,a2)
#       temp = benchdb(ranking, N)
#       View(temp)
#       View(res)
#       if(is.null(res)) res = temp
#       else res = cbind(res, temp)
#     }
#   }
#   res
#   View(res)
# }


grid.search <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, N = 50){
  res = NULL
  for(a1 in seq(0,1,0.1)){
    for(a2 in seq(0,1-a1,0.1)){
      a3 = 1 - (a1+a2)
      ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
      temp = benchdb(ranking, N)
      if(is.null(res)) res = temp
      else res = cbind(res, temp)
    }
  }
  res
}


makeCVDataSets <- function(dd  = NULL,
                          pd   = NULL,
                          list = NULL){
  idx = which(rownames(dd) %in% list)
  #print(length(idx))
  expd=t(dd[idx,])
  t1 = as.numeric(pd[,1])
  t1[is.na(t1)] <- 0
  t1[t1<0] <- 0
  t2 = as.numeric(pd[,2])
  pd = data.frame(cbind(t1,t2))
  my_pd = pd
  View(my_pd)
  colnames(pd) = c("time", "status")
  if(nrow(expd) == 1) {data = cbind(t(expd), pd)
  colnames(data) = c("G1","time", "status")
  }
  else data = cbind(expd, pd)
  CV_method_data = data
  #CV_method_data$status = as.logical(CV_method_data$status)
  CV_method_data = CV_method_data[-c(52, 378), ]
  View(CV_method_data)
  datasets <- makeSurvTask(data = CV_method_data, target = c("time", "status"))
  View(datasets)
  return(datasets)
}

doPHModel <- function(train.task = NULL,
                     test.task  = NULL,
                     n.iter     = 100) {
  res.aggr = matrix(data = NA, nrow = n.iter, ncol = 2)
  colnames(res.aggr) = c("ci", "lgrk")
  for (i in 1:n.iter) {
    tryCatch({
      print(i)
      #set.seed(seed[i])
      lrn <- makeLearner(cl = "surv.coxph")
      print(lrn)
      rdesc  <-  makeResampleDesc(method = "CV", iters = 10L)
      print(rdesc)
      res = resample(lrn, train.task, rdesc, measures = list(ci, lgrk),show.info = FALSE)
      print(res)
      temp = res$aggr
      ####test join from ci
      # resp = res$pred
      # a = Hmisc::rcorr.cens(x = getPredictionResponse(resp),
      #                       S = getPredictionTruth(resp))
      # res.list[[i]] = list("ci" = a[["C Index"]], "mean.aggr" = temp)
      res.aggr[i,"ci"] = temp[1]
      res.aggr[i,"lgrk"] = temp[2]
      #hrph.vec[i] = temp[3]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  means = colMeans(res.aggr,na.rm = TRUE)
  return (means)
}

run_topN <- function(mRNA    = NULL,
                    pd      = NULL,
                    ranking = NULL,
                    N       = 50){
  my_ranking = names(ranking) ###notice this line
  index = which(my_ranking %in% rownames(mRNA))
  t_rank = my_ranking[index]
  ddata = mRNA[t_rank,]
  genes = t_rank[1:N]
  dataSets = makeCVDataSets(ddata,pd,genes)
  res = doPHModel(dataSets)
  res
}

benchdb <- function(ranking = NULL,
                   N       = 50) {
  pd = my_survival_data[,c(3,1)] # rf
  A = run_topN(mRNA = RNASeq_data, pd = pd, ranking = ranking,N = N)
}

#Loading the single-cell data files. ----

#all_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_NM_all_cells_FPKM.csv")
#epi_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_NM_epithelial_cells_FPKM.csv")
all_tumor_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_tumor_all_cells_FPKM.csv")
#epi_cells_tumor_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_tumor_epithelial_cells_FPKM.csv")
#geo_ega_id_match <- read.csv("Data/Single-cell-data/GSE81861_GEO_EGA_ID_match.csv")



#Loading the bulk RNA seq files. ----
#Define query.
colon_query <- GDCquery(project       = "TCGA-COAD",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

#Downloading the data.
GDCdownload(query           = colon_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/Bulk-data/TCGA-Colon-Cancer-Dataset")

#Making the summarizedExperiment object and then removing all entries that lacked days_to_last_follow_up information
COA_data_se <- GDCprepare(colon_query, summarizedExperiment = TRUE, directory = "Data/Bulk-data/TCGA-Colon-Cancer-Dataset/")
colon_ind <- is.na(COA_data_se$days_to_last_follow_up)
COA_data_se$days_to_last_follow_up[colon_ind] <- COA_data_se$days_to_death[colon_ind]



#Turning the summarizedExperiment object into a data frame and re-factoring the labels on the vital_status column from Alive/Dead to 0/1.
#We then turn them into the boolean values TRUE/FALSE for use in later calculations. 
COA_data_df <- as.data.frame(colData(COA_data_se))
COA_data_df$vital_status <- factor(COA_data_df$vital_status, levels = c("Alive", "Dead"),
                                         labels = c(0,1))

COA_data_df$vital_status <- as.numeric(as.character(COA_data_df$vital_status))
COA_survival_data <- subset(COA_data_df, select=c("vital_status", "days_to_last_follow_up"))
COA_survival_data$vital_status <- factor(COA_survival_data$vital_status, levels = c(0,1),
                                         labels = c(FALSE,TRUE))

# Building the bulk RNA dataframe to merge with clinical data to run Cox PH on
#We do this by extracting various components of the summarizedExperiment object that was made earlier.
bulk_rna_df <- COA_data_se@assays@data@listData[["HTSeq - Counts"]]
colnames(bulk_rna_df) <- COA_data_se@colData@rownames
rownames(bulk_rna_df) <- COA_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames

#Making a pseudotime of the all tumor cells from the scRNA-seq dataset (data is already pre-processed).
#We do this through the magic() function. We use the seed parameter set to '123' to make the results reproducible.
#We first turn the single-cell data into a data frame and get the rownames of the frame (gene names).
#We then subset the data frame to remove the unnecessary column that contains the rownames in a column now that 
#they are stored in the actual rownames. The input of the magic() function requires the data to be in matrix
#format. We then transpose the results and do some basic filtering. The filtering involves filtering to include 
#rows that are greater than 0.3 and 0. This mean value is then compared to 0.2 and must be greater than 0.2 to be
#kept in the dataset. We then convert this output back to a dataframe and include a new name for the first colum
#Which includes the short names of the genes. Once this is finished then we finally convert the dataframe back to
#a matrix for future calculations. The end product is a matrix that contains 375 cells (columns) and 23,479 genes (rows).
all_tumor_cells_fpkm_df <- as.data.frame(all_tumor_cells_fpkm)
non_denoised_gene_names <- all_tumor_cells_fpkm_df$X
non_denoised_gene_names <- as.data.frame(non_denoised_gene_names)
colnames(non_denoised_gene_names)[1] <- "gene_short_name"
rownames(all_tumor_cells_fpkm_df) <- non_denoised_gene_names$gene_short_name
all_tumor_cells_fpkm_denoised <- subset(all_tumor_cells_fpkm_df, select=RHC3546__Tcell__.C6E879:RHC6041__Macrophage__.FFFF55)
all_tumor_cells_fpkm_denoised <- as.matrix(all_tumor_cells_fpkm_denoised)
all_tumor_cells_fpkm_denoised <- magic(t(all_tumor_cells_fpkm_denoised), seed = 123)
all_tumor_cells_fpkm_denoised <- t(all_tumor_cells_fpkm_denoised[["result"]])
all_tumor_cells_fpkm_denoised <- all_tumor_cells_fpkm_denoised[rowMeans(all_tumor_cells_fpkm_denoised) > 0.3 & rowMeans(all_tumor_cells_fpkm_denoised > 0) > 0.2,]
all_tumor_cells_fpkm_denoised <-  as.data.frame(all_tumor_cells_fpkm_denoised)
tumor_gene_names <- rownames(all_tumor_cells_fpkm_denoised)
tumor_gene_names <- as.data.frame(tumor_gene_names)
colnames(tumor_gene_names)[1] <- "gene_short_name"
rownames(tumor_gene_names) <- tumor_gene_names$gene_short_name
all_tumor_cells_fpkm_denoised <- as.matrix(all_tumor_cells_fpkm_denoised)

#Doing some pre-processing of the rownames to make stuff cleaner----
current_rowname_split <- strsplit(rownames(all_tumor_cells_fpkm_denoised), "_")

finished_gene_list <- c()
current_list <- current_rowname_split
for (y in seq(1:length(current_list))){
  #print(current_list[[y]][2])
  finished_gene_list <- c(finished_gene_list, current_list[[y]][2])
}

#rownames(all_tumor_cells_fpkm_denoised) <- rownames(finished_gene_list)
# finished_gene_list <- unique(finished_gene_list)
# finished_gene_list <- as.data.frame(finished_gene_list)
# colnames(finished_gene_list)[1] <- "gene_short_name"
# rownames(finished_gene_list) <- finished_gene_list$gene_short_name

#all_tumor_cells_fpkm_denoised <- filter(all_tumor_cells_fpkm_denoised, rownames(all_tumor_cells_fpkm_denoised) %in% rownames(finished_gene_list))

#Moncole3 steps ----
#Here We first are making some vectors of VIM genes (mesencymal progression marker) and CDH genes (epithelial marker) that I can use to
#determine the starting point of the pseudotime calculation carried out by the Moncocle3 package. We make a new cell data set and run
#through the typical steps to get the pseudotime calculated by Monocle3. This includes steps such as dimensionality reduction through
#PCA and UMAP to help cluster the cells and learn their order orientation in space. I plot the VIM and CDH1 expression levels to
#understand where to set the start point of the pseudotime calculation. I then plot the graph with the pseudotime calculated
#to get a view of what it looks like. 
vim_genes <- c("chr10:17256237-17279592_VIM_ENSG00000026025.9", "chr15:101811021-101817705_VIMP_ENSG00000131871.10", "chr6:126923501-126924795_VIMP1_ENSG00000220548.3", "chr10:17256237-17279592_VIM-AS1_ENSG00000229124.2")
cdh1_genes <- c("chr5:141232937-141258811_PCDH1_ENSG00000156453.9", "chr16:68771127-68869451_CDH1_ENSG00000039068.14")
cell_data_set <- new_cell_data_set(all_tumor_cells_fpkm_denoised,gene_metadata=tumor_gene_names)
cds <- preprocess_cds(cell_data_set, num_dim=100, method="PCA")
cds <- reduce_dimension(cds)
plot_cells(cds)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by="partition")
plot_cells(cds, genes=vim_genes, cell_size=0.50)
plot_cells(cds, genes=cdh1_genes, cell_size=0.50)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds                 = cds,
           color_cells_by      = "pseudotime",
           label_cell_groups   = FALSE,
           label_leaves        = FALSE,
           label_branch_points = FALSE,
           graph_label_size    = 1.5)


#Integrating the new pseudotime data into the count matrix. Here I extract just the pseudotime calculations from the cds
#object to a use for the swithde() calculation later. The output of this code is a dataframe that has 2 columns (pseudotime
#value and sample name[cell name]) and 375 #rows (cells with their associated pseudotime).
pseudotime_data <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudotime_data <- as.data.frame(pseudotime_data)
colnames(pseudotime_data) <- "Pseudotime"
pseudotime_contents <- rownames(pseudotime_data)
pseudotime_data$Samples <- pseudotime_contents

#Different metrics ----
#3.1 The median absolute deviation metric (MAD) ----
mads <- apply(all_tumor_cells_fpkm_denoised,1,mad)
index <- order(mads, decreasing=TRUE)
mad.ranking<- mads[index]
mad.ranking<-abs(mad.ranking)/sum(abs(mad.ranking))
save(mad.ranking, file = "Data/Exported-data/R-objects/mad.ranking.RData")

#3.2 The switchde (SDE) metric ----
sdes <- switchde(all_tumor_cells_fpkm_denoised, as.numeric(pseudotime_data$Pseudotime), verbose = TRUE)
sde.filtered <- filter(sdes, qval < 0.05)
index <- order(abs(sde.filtered$k), decreasing = T)
vim.sdes.rank <- sde.filtered[index,]
vim.sdes.ranking <- vim.sdes.rank$k
names(vim.sdes.ranking) <- vim.sdes.rank$gene
vim.sdes.ranking<-abs(vim.sdes.ranking)/sum(abs(vim.sdes.ranking))
save(vim.sdes.ranking, file = "Data/Exported-data/R-objects/vim.sdes.ranking.RData")

#3.3 The miRNA metric ----
#Loading the microRNA data from Mary----
my_miRNA <- read.table(file="Data/miRNA-data/miRNA_list.txt", sep="\t", header=FALSE)
my_miRNA <- my_miRNA[1]
colnames(my_miRNA) <- "miRNAs"
my_miRNA <- my_miRNA$miRNAs
my_miRNA <- as.vector(my_miRNA)
my_miRNA <- trimws(my_miRNA)


#miRNAs from dbDEMC version 2.0----
dbDEMC_high <- read.csv(file = "Data/miRNA-data/List-of-dbDEMC-2-0-miRNAs/dbDEMC-2.0-high.txt", sep = '\t')
dbDEMC_low <- read.csv(file = "Data/miRNA-data/List-of-dbDEMC-2-0-miRNAs/dbDEMC-2.0-low.txt", sep = '\t')
colnames(dbDEMC_low)[1] <- "miRNA.ID"

#Filtering to just the miRNAs associated with colorectal or colon cancer. 
dbDEMC_high <- filter(dbDEMC_high, Cancer.Type=="colorectal cancer" | Cancer.Type=="colon cancer")
dbDEMC_high_miRNAs <- subset(dbDEMC_high, select = miRBase.Update.ID)
dbDEMC_high_miRNAs <- as.vector(dbDEMC_high_miRNAs)

#miRNAs from miRmap----
miRmap_transcripts <- read.csv(file = "Data/miRNA-data/MiRMap-data/mirmap201301e_homsap_transcripts.csv", sep = ',')
miRmap_mirnas <- read.csv(file = "Data/miRNA-data/MiRMap-data/mirmap201301e_homsap_mirnas.csv", sep = ',')

#f <- function(x, pos) print(pos)
#miRmap_targets <- read_csv_chunked(file = "Data/miRNA-data/MiRMap-data/mirmap201301e_homsap_targets.csv", chunk_size = 1000, callback = DataFrameCallback$new(f))


#Common miRNAs between databases----
#Intersection between the miRNAs of the two databases
common_mirnas <- intersect(miRmap_mirnas$mature_name, dbDEMC_high_miRNAs$miRBase.Update.ID)

#Now submitting these miRNAs to TargetScan to get genes to make a gene list for the third metric----
my_num <- 1
miRNA_targets <- list()
for (m in common_mirnas[1:300]) {
  print(m)
  current_target <- targetScan(mirna=common_mirnas[my_num], species="Human", release="7.2", maxOut= NULL)
  miRNA_name <- m
  miRNA_name_final <- rep(miRNA_name, times=length(current_target$Ortholog))
  current_target <- cbind(current_target,miRNA_name_final)
  miRNA_targets[[m]] <- current_target
  my_num <- my_num + 1
  print(my_num)
}

#Simplifying the output of the targetscan commands to just the Gene name
#and the mirna columns
counter <- 1
total_list <- list()
for (i in miRNA_targets){
  current_df <- i
  gene_list <- current_df$Ortholog
  mirna_list <- current_df$miRNA_name_final
  simple_list <- cbind(gene_list,mirna_list)
  total_list[[counter]] <- simple_list
  counter <- counter + 1
}

counter <- 1
all_miRs_for_score <- list()
all_genes_for_score <- list()
for (l in total_list){
  current_df <- l
  current_miRs <- current_df[,'mirna_list']
  all_miRs_for_score[[counter]] <- unique(current_miRs)
  current_genes <- current_df[,'gene_list']
  all_genes_for_score[[counter]] <- current_genes
  counter <- counter + 1
}

all_miRs_for_score <- unlist(all_miRs_for_score)
all_genes_for_score <-unlist(all_genes_for_score)
all_genes_for_score_unique <- unique(all_genes_for_score)

miRNA_score <- matrix(data = 0, nrow = length(all_genes_for_score_unique), ncol = length(all_miRs_for_score), dimnames = list(all_genes_for_score_unique, all_miRs_for_score))
miRNA_score <- as.data.frame(miRNA_score)

#Checking to see for each miRNA (colname) if it interacts with a particular row.
#If it does it get a plus one to that cell. If it does not it moves to next cell
for (i in rownames(miRNA_score)){
  for (x in miRNA_targets){
    current_df <- x
    if (i %in% current_df$Ortholog){
      miRNA_to_add <- unique(current_df$miRNA_name_final)
      miRNA_score[i,miRNA_to_add]<- miRNA_score[i, miRNA_to_add] + 1
    }
    
  }
}

#Now calculating the rowsums of each gene for total number of miRNA interactions----
mirna_gene_list <- rowSums(miRNA_score)
mirna_gene_list <- as.data.frame(mirna_gene_list)
colnames(mirna_gene_list)[1] <- "Counts"
mirna_gene_list <- arrange(mirna_gene_list, desc(Counts))
mirna_gene_list <- as.vector(mirna_gene_list)
mirna.ranking<-abs(mirna_gene_list)/sum(abs(mirna_gene_list))
colnames(mirna.ranking)[1] <- "Score"

#Converting the gene symbols to Ensembl gene ids to have a common gene naming----
#system with other metrics
library(EnsDb.Hsapiens.v79)

# 2. Convert from gene.symbol to ensembl.gene
geneSymbols <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')

geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

mirna.ranking <- as.vector(mirna.ranking)
save(mirna.ranking, file = "Data/Exported-data/R-objects/mirna.ranking.RData")
#write.csv(mirna.ranking, file = "Data/Exported-data/mirna-ranking.csv")


#Loading in a large set of well known EMT genes----
emt_genes <- read.csv(file = "Data/EMT-gene-data/emt-hallmark-genes.csv", sep = ',')
common_genes <- intersect(head(emt_genes$Genes, n=20), rownames(mirna.ranking))

#Loading in a defined set of colon cancer genes----
cc_genes <- read.csv(file = "Data/Colon-cancer-markers/cc-markers.csv", sep = ',')
common_cc_genes <- intersect(cc_genes$Markers, head(emt_genes$Genes, n =80))


#All of this code is on hiatus until a further date----
# #Reading in the miRNA similarity score information
# library(readxl)
# miRNA_similarity_names <- read_excel("Data/miRNA-data/miRNA-similarity-data/microRNA name.xls", col_names = FALSE)
# colnames(miRNA_similarity_names) <- "miRNAs"
# rownames(miRNA_similarity_names) <- miRNA_similarity_names$miRNAs
# miRNA_similarity <- read.csv("Data/miRNA-data/miRNA-similarity-data/miRNA similarity matrix.txt", sep = '\t')
# for (x in seq(1:length(colnames(miRNA_similarity)))){
#   colnames(miRNA_similarity)[x] <- rownames(miRNA_similarity_names)[x]
# }
# 
# #This next section will be dedicated to getting the disease similarity part setup
# 
# 
# # Using the hoardeR package targetScan() function to get the miRNA targets.
# my_num <- 1
# remove <- c("miR-451", "miR-21", "miR-221", "miR-154-3p/487-3p", "miR-487a")
# miRNA_targets <- list()
# my_miRNA[1] <- "let-7-5p/98-5p"
# my_miRNA[4] <- "miR-15a-3p"
# my_miRNA[5] <- "miR-15-5p"
# my_miRNA[12] <- "miR-30-5p"
# my_miRNA[14] <- "miR-25-3p/32-5p/92-3p/363-3p/367-3p"
# my_miRNA[16] <- "miR-34-5p/449-5p"
# my_miRNA[41] <- "miR-15-5p/16-5p/195-5p/424-5p/497-5p"
# my_miRNA[45] <- "miR-204-5p/211-5p"
# my_miRNA[52] <- "miR-302-3p/372-3p/373-3p/520-3p"
# my_miRNA[56] <- "miR-302-3p/372-3p/373-3p/520-3p"
# my_miRNA[77] <- "miR-96-5p/1271-5p"
# my_miRNA[80] <- "miR-96-5p/1271-5p"
# my_miRNA[86] <- "miR-21-5p/590-5p"
# 
# 
# my_miRNA <- setdiff(my_miRNA, remove)
# my_miRNA <- trimws(my_miRNA)
# 
# for (m in my_miRNA[1:121]) {
#   print(m)
#   current_target <- targetScan(mirna=my_miRNA[my_num], species="Human", release="7.2", maxOut= 50)
#   miRNA_name <- m
#   miRNA_name_final <- rep(miRNA_name, times=length(current_target$Ortholog))
#   current_target <- cbind(current_target,miRNA_name_final)
#   miRNA_targets[[m]] <- current_target
#   my_num <- my_num + 1
#   print(my_num)
# }
# all_miRNA_targets <- do.call(rbind, miRNA_targets)
# all_miRNA_targets <- arrange(all_miRNA_targets, desc(consSites))
# miRNA_gene_intersect <- intersect(all_miRNA_targets$Ortholog, genes_in_bulk_RNA)
# miRNA_gene_intersect

#Now doing the grid search of optimal values for the linear, integrated model----
#Currently just two metrics (MAD and SDE).
weights <- seq(from = 0, to=1, by=0.1)
df_index <- 1
integrated_gene_lists <- list()
a3 <- 0.1

for (x in weights) {
  print(x)
  current_ranking <- geneRank(ranking1 = mad.ranking, ranking2 = vim.sdes.ranking, ranking3 = mirna.ranking,  a1=x, a2=1-(x+a3), a3= 0.1)
  current_ranking <- as.data.frame(current_ranking)
  integrated_gene_lists[[df_index]] <- current_ranking
  df_index <- df_index + 1
}
#Subsetting the gene expression data frame to just the top N genes found from the integrated ranking at different combinations/values of weights ----
gene_lists_to_test <- list()
all_tumor_cells_fpkm_denoised_df <- as.data.frame(all_tumor_cells_fpkm_denoised)
all_tumor_cells_fpkm_denoised_df <- t(all_tumor_cells_fpkm_denoised_df)
for (x in seq(1:length(integrated_gene_lists))){
  print(x)
  current_gene_list <- as.data.frame(integrated_gene_lists[[x]])
  current_gene_list$GeneName <- rownames(current_gene_list)
  gene_names_to_test <- head(current_gene_list, n = 900)
  gene_lists_to_test[[x]] <- gene_names_to_test
}

genes_of_interest <- list()
for (x in seq(1:length(integrated_gene_lists))){
  current_list <- gene_lists_to_test[[x]]
  current_list <- as.data.frame(current_list)
  current_genes <- subset(all_tumor_cells_fpkm_denoised_df, select = current_list$GeneName) 
  genes_of_interest[[x]] <- current_genes
}

all_split_colnames <- list()
for (x in seq(1:length(integrated_gene_lists))){
  for (y in seq(1:length(colnames(genes_of_interest[[x]])))){
    current_df <- genes_of_interest[[x]]
    current_df <- as.data.frame(current_df)
    current_colname_split <- strsplit(colnames(current_df), "_")
    all_split_colnames[[x]] <- current_colname_split
  }
  
}

finished_sets <- list()
for (x in seq(1:length(integrated_gene_lists))) {
  finished_gene_list <- c()
  current_list <- all_split_colnames[[x]]
  for (y in seq(1:length(current_list))){
    #print(current_list[[y]][2])
    finished_gene_list <- c(finished_gene_list, current_list[[y]][2])
  }
  finished_gene_list <- unique(finished_gene_list)
  finished_sets[[x]] <-finished_gene_list
}

#Making the dataframe that contains just the genes by samples that we need from bulk RNA-seq data
gene_expression_info <- COA_data_se@assays@data@listData[["HTSeq - Counts"]]
rownames(gene_expression_info) <- COA_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
colnames(gene_expression_info) <- COA_data_se@colData@rownames
gene_expression_info <- t(gene_expression_info)


all_intersections <- list()
for (x in seq(1:length(finished_sets))){
  current_set <- finished_sets[[x]]
  current_intersect <- intersect(colnames(gene_expression_info), current_set)
  all_intersections[[x]] <- current_intersect
}



#Merging the two dataframes together into a larger dataframe that we can use for the Cox PH
bulk_rna_df_unique <- subset(bulk_rna_df, select = unique(colnames(bulk_rna_df)))
COA_data_df_unique <- subset(COA_data_df, select = unique(colnames(COA_data_df)))
merged_df <- merge(bulk_rna_df_unique, COA_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]

#Making the dataframe with the survival info
survival_df <- subset(merged_df, select=vital_status)
time <- as.numeric(merged_df$days_to_last_follow_up)
survival_df$time <- cbind(time)
survival_df <- na.omit(survival_df)



#69 gene list--Concordance = 0.751 (se = 0.025)
#58 gene list--Concordance = 0.735 (se = 0.024)
#48 gene list--Concordance = 0.705 (se = 0.025)
#38 gene list--Concordance = 0.70  (se = 0.025)
#28 gene list--Concordance = 0.668 (se = 0.026)
#18 gene list--Concordance = 0.618 (se = 0.028)
#8 gene list--Concordance  = 0.533 (se = 0.029)

#c2 <- colnames(merged_df)[2:40]


all_intersections_cleaned <- list()
for (x in seq(1:length(all_intersections))){
  genes_in_bulk_RNA <- all_intersections[[x]]
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="-",replacement=".")
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="_", replacement=".")
  genes_in_bulk_RNA <- sapply(genes_in_bulk_RNA, gsub, pattern="/", replacement=".")
  all_intersections_cleaned[[x]] <- genes_in_bulk_RNA
}

rows_to_remove <- setdiff(rownames(merged_df), rownames(survival_df))
merged_df <- merged_df[!(row.names(merged_df) %in% rows_to_remove), ]


colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern="-",replacement=".")
colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern="_", replacement=".")
colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern="/", replacement=".")

#Now doing data splitting for training and testing sets----
df_for_train_test_split <- merge(merged_df, survival_df, by="row.names")
df_for_train_test_split <- subset(df_for_train_test_split, select=-c(Row.names, vital_status, time))
df_for_train_test_split <- filter(.data = df_for_train_test_split, days.to.last.follow.up!=0)
my_time <- df_for_train_test_split$days.to.last.follow.up
my_status <- df_for_train_test_split$vital.status
df_for_train_test_split <- subset(df_for_train_test_split, select=c(TSPAN6:AC007389.3))
df_for_train_test_split$days.to.last.follow.up <- my_time
df_for_train_test_split$vital.status <- my_status
set.seed(1)
index <- createDataPartition(df_for_train_test_split$vital.status, p = 0.6, list = F)
merged_train <- df_for_train_test_split[index, ] # 60%
merged_test <- df_for_train_test_split[-index, ] # 40%


#Cox models----
cox_models <- list()
f_objects <- list()
lambdas <- list()
c_indicies <- list()
prediction_res <- list()
all_coefs <- list()
all_active_coefs <- list()

for (x in seq(1:length(all_intersections_cleaned))){
  current_formula_data <- all_intersections_cleaned[[x]]
  current_formula_data <- as.vector(current_formula_data)
  current_formula_data <- current_formula_data[-107]
  
  f <- as.formula(paste("~", paste(current_formula_data[1:length(current_formula_data)], collapse = "+")))
  
  
  
  my_x <- model.matrix(f, df_for_train_test_split)
  my_y <- Surv(time = df_for_train_test_split$days.to.last.follow.up, event = df_for_train_test_split$vital.status)
  cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = 10, type.measure = "C", maxit=100000, family="cox")
  current_lambdas <- cv_fit$lambda
  current_c_indicies <- cv_fit$cvm
  
  
  
  cox_models[[x]] <- cv_fit
  f_objects[[x]] <- f
  lambdas[[x]] <- current_lambdas
  c_indicies[[x]] <- current_c_indicies

  
  #Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family =  "cox", maxit = 100000)
  Coefficients <- coef(fit, s = cv_fit$lambda.min)
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  all_coefs[[x]] <- Coefficients
  all_active_coefs[[x]] <- Active.Coefficients
  
}

for(x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  print(current_model)
}

all_c_indexes_finished <- list()
all_lambdas <-list()
min_lambda <- list()
for (x in seq(1:length(cox_models))){
  current_model <- cox_models[[x]]
  all_c_indexes_finished[[x]] <- current_model$cvm
  all_lambdas[[x]] <- current_model$lambda
  min_lambda[[x]] <- current_model$lambda.min
}


my_newdata <- merged_test
my_newdata <- model.matrix(f_pred,data = my_newdata)
my_model_pred <- predict(cv_fit, s = cv_fit$lambda.min, newx =my_newdata, type.measure="C")

my_assesment <- assess.glmnet(object = cv_fit, newx = my_newdata, newy = Surv(time = merged_test$days.to.last.follow.up, event = merged_test$vital.status),
                              family = "cox", s=cv_fit$lambda.min, type.measure="C")

my.filename <- "10-fold-cv-active-coefs-for-all-11-cox-models-60-40-split.csv"
all_cvs <- do.call(rbind, all_c_indexes_finished)
all_lambdas_to_save <- do.call(rbind, all_lambdas)
all_min_lambdas <- do.call(rbind, min_lambda)
all_coefs_merged <- do.call(rbind, all_coefs)
all_active_coefs_merged <- do.call(rbind, all_active_coefs)
#test_preds <- do.call(rbind, prediction_res)
#colnames(test_preds)[1] <- "Predicted C-index on Test data"
# write.csv(all_active_coefs_merged, file = my.filename, row.names = TRUE, quote = FALSE)



#Doing heatmap of weights for linear model with the color shading indicating what the c-index value is----
for (x in cox_models){
  print(x)
}


#With random genes selected of the same size as our N significant ones ----
#N times, repeatedly sample rows from the data and write to a csv file
total_files <- 15000
df_index <- 1
concordance_value <- list()
concordance_value_se <- list()
logtest_value <- list()
waldtest_value <- list()
sctest_value <- list()
current_gene_list <- list()
gene_set_size <- 111

#Making a large pool of randomly drawn genes from the bulk RNA-seq data
sampled_subset <- sample(colnames(merged_df[,2:56567]), size=total_files*gene_set_size, replace =TRUE)
sampled_subset <- unique(sampled_subset)

while (df_index <= total_files) {
  #Now I am randomly drawing genes from that generated subset. 
  current_subset <- sample(sampled_subset, size = gene_set_size, replace = FALSE)
  
  
  # Write the data to a csv file
  filename_to_use <- paste0("Sample_Data", "_", df_index, ".csv")
  
  # write.csv(x = current_subset,
  #           file = filename_to_use)
  
  random_genes_in_bulk_data <- intersect(current_subset, colnames(gene_expression_info))
  random_genes_in_bulk_data <- sapply(random_genes_in_bulk_data, gsub, pattern="-",replacement=".")
  random_genes_in_bulk_data <- sapply(random_genes_in_bulk_data, gsub, pattern="_", replacement=".")
  random_genes_in_bulk_data <- sapply(random_genes_in_bulk_data, gsub, pattern="/", replacement=".")
  random_genes_in_bulk_data <- as.vector(random_genes_in_bulk_data)
  numeric_starters<-grep('^[0-9]',random_genes_in_bulk_data, value = TRUE)
  numeric_starters<-gsub("^[[:digit:]]", "AW", numeric_starters)
  numeric_starters_index <-grep('^[0-9]',random_genes_in_bulk_data, value = FALSE) 
  random_genes_in_bulk_data[numeric_starters_index] <- numeric_starters
  colnames(merged_df)<-gsub("-",".",colnames(merged_df))
  colnames(merged_df)<-gsub("_",".",colnames(merged_df))
  colnames(merged_df)<-gsub("/",".",colnames(merged_df))
  numeric_starters<-grep('^[0-9]',colnames(merged_df), value = TRUE)
  numeric_starters<-gsub("^[[:digit:]]", "AW", numeric_starters)
  numeric_starters_index <-grep('^[0-9]',colnames(merged_df), value = FALSE) 
  colnames(merged_df)[numeric_starters_index] <- numeric_starters
  
  f2 <- as.formula(paste("surv_object_se2 ~ ",
                        paste(random_genes_in_bulk_data[1:length(random_genes_in_bulk_data)], collapse= "+")))
  
  
  res.cox6 <-coxph(f2, data=merged_df)
  summary_of_random_sigs <- summary(res.cox6)
  
  
  
  current_concordance_value <- as.data.frame(res.cox6$concordance["concordance"])
  current_concordance_value_se <- as.data.frame(res.cox6$concordance["std"])
  current_logtest_value <- as.data.frame(summary_of_random_sigs$logtest["pvalue"])
  current_waldtest_value <- as.data.frame(summary_of_random_sigs$waldtest["pvalue"])
  current_sctest_value <- as.data.frame(summary_of_random_sigs$sctest["pvalue"])
  current_gene_value <- as.data.frame(random_genes_in_bulk_data)
  
  concordance_value[[df_index]] <- current_concordance_value
  concordance_value_se[[df_index]] <- current_concordance_value_se
  logtest_value[[df_index]] <- current_logtest_value
  waldtest_value[[df_index]] <- current_waldtest_value
  sctest_value[[df_index]] <- current_sctest_value
  current_gene_list[[df_index]] <- current_gene_value
  
  
  #Writing the files out to .csvs
  # write.table(subset_of_random_sigs$coefficients, "p-values-for-N-random-sig-genes.csv", append = TRUE, sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)
  # write.table(subset_of_random_sigs$concordance, "concordance-value-for-N-random-sig-genes.csv", append = TRUE, sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)
  # write.table(subset_of_random_sigs$logtest, "logtest-for-N-random-sig-genes.csv", append = TRUE, sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)
  # write.table(subset_of_random_sigs$sctest, "sctest-for-N-random-sig-genes.csv", append = TRUE, sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)
  # write.table(subset_of_random_sigs$waldtest, "waldtest-for-N-random-gene-sigs.csv", append = TRUE, sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  df_index <- df_index + 1
}

cox_results_concordance_vector <- do.call(rbind, concordance_value)
cox_results_concordance_se_vector <- do.call(rbind, concordance_value_se)
cox_results_logtest_vector <- do.call(rbind, logtest_value)
cox_results_waldtest_vector <- do.call(rbind, waldtest_value)
cox_results_sctest_vector <- do.call(rbind, sctest_value)
final_concordance_df <- cbind(cox_results_concordance_vector, cox_results_concordance_se_vector, cox_results_logtest_vector, cox_results_waldtest_vector, cox_results_sctest_vector)
colnames(final_concordance_df) <- c("Concordance Value (C)", "Concordance Value SE", "Log Test P-value", "Wald Test P-value", "Sc Test P-value")
rownames(final_concordance_df) <- c(1:total_files)


#Plotting for random simulations----
hist(final_concordance_df$`Concordance Value (C)`, xlab = "Concordance Value (C)", main = "Concordance Value (C) Histogram for 150 Gene Set 15000 Simulations")
boxplot(final_concordance_df$`Concordance Value (C)`)
plot(final_concordance_df$`Concordance Value (C)`, final_concordance_df$`Wald Test P-value`, xlab="Concordance Value (C)", ylab="Wald Test P-value", main="Concordance Values vs. Wald Test P-values")
plot(final_concordance_df$`Concordance Value (C)`, final_concordance_df$`Log Test P-value`, xlab="Concordance Value (C)", ylab="Log Test P-value", main="Concordance Values vs. Log Test P-values")
plot(final_concordance_df$`Concordance Value (C)`, final_concordance_df$`Sc Test P-value`, xlab="Concordance Value (C)", ylab="Sc Test P-value", main="Concordance Values vs. Sc Test P-values")

gene_list_vector <- do.call(rbind, current_gene_list)
unique_gene_list_vector <- unique(gene_list_vector)
percent_same_across_runs <- abs(nrow(gene_list_vector)-nrow(unique_gene_list_vector))/(nrow(gene_list_vector))*100
percent_same_across_runs

genes_in_bulk_RNA_vector <- as.vector(genes_in_bulk_RNA)
intersect_my_genes_with_random <- intersect(unique_gene_list_vector$random_genes_in_bulk_data, genes_in_bulk_RNA_vector)
num_in_both <- length(intersect_my_genes_with_random)
num_in_my_list <- length(genes_in_bulk_RNA_vector)
percent_shared <- (num_in_both/(num_in_my_list-1))* 100
percent_shared


#Comparison to other methods ----

#1. Xu et al. method ----
Xu_gene_sigs <- c("HES5", "ZNF417", "GLRA2", "OR8D2", "HOXA7", "FABP6", "MUSK", "HTR6", "GRIP2", "KLRK1", "VEGFA", "AKAP12", "RHEB", "NCRNA00152", "PMEPA1")

Xu_shared_with_my_set <- intersect(Xu_gene_sigs, genes_in_bulk_RNA_vector)

Xu_shared_with_my_bulk_data <- intersect(Xu_gene_sigs, colnames(merged_df))

#Xu_df <- subset(merged_df, colnames(merged_df) %in% Xu_shared_with_my_bulk_data)





f3 <- as.formula(paste("surv_object_se2 ~ ",
                       paste(Xu_shared_with_my_bulk_data[1:length(Xu_shared_with_my_bulk_data)], collapse= "+")))


res.cox.xu <- coxph(f3, data=merged_df)
summary_of_xu_sigs <- summary(res.cox.xu)


#2. Barrier et al. method ----
Barrier_nm_sig_genes <- read.table(file = "Barrier-paper-nm-genes.txt", sep = "\t")
Barrier_t_sig_genes <- read.table(file = "Barrier-paper-tumor-genes.txt", sep = "\t")

Barrier_nm_sig_genes <- subset(Barrier_nm_sig_genes, select="V1")
Barrier_t_sig_genes <- subset(Barrier_t_sig_genes, select="V1")

# Converting PROBEIDs to Gene name and symbols
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = Barrier_nm_sig_genes$V1, uniqueRows=TRUE)

Barrier_nm_sig_genes_gn <- annotLookup$external_gene_name

Barrier_nm_shared_with_my_set <- intersect(Barrier_nm_sig_genes_gn, genes_in_bulk_RNA_vector)
Barrier_nm_shared_with_bulk_data <- intersect(Barrier_nm_sig_genes_gn, colnames(merged_df))

f4 <- as.formula(paste("surv_object_se2 ~ ",
                       paste(Barrier_nm_shared_with_bulk_data[1:length(Barrier_nm_shared_with_bulk_data)], collapse= "+")))


res.cox.barrier <- coxph(f4, data=merged_df)
summary_of_barrier_sigs <- summary(res.cox.barrier)



require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = Barrier_t_sig_genes$V1, uniqueRows=TRUE)

Barrier_t_sig_genes_gn <-annotLookup$external_gene_name

Barrier_t_shared_with_my_set <- intersect(Barrier_t_sig_genes_gn, genes_in_bulk_RNA_vector)
Barrier_t_shared_with_bulk_data <- intersect(Barrier_t_sig_genes_gn, colnames(merged_df))

f5 <- as.formula(paste("surv_object_se2 ~ ",
                       paste(Barrier_t_shared_with_bulk_data[1:length(Barrier_t_shared_with_bulk_data)], collapse= "+")))


res.cox.barrier_tumor <- coxph(f5, data=merged_df)
summary_of_barrier_sigs_tumor <- summary(res.cox.barrier_tumor)
  
#3. Goa et al. gene sets ----
goa_analyzer <- function(my.file,my.column){
  subset_list <- list()
  goa_genes <- read.table(file = my.file, sep = "\t", header = TRUE)
  goa_subset <- subset(goa_genes, select=my.column)
  goa_subset <- goa_subset[1:30,1]
  
  goa_shared_with_my_set <- intersect(goa_subset, genes_in_bulk_RNA_vector)
  goa_shared_with_bulk_data <- intersect(goa_subset, colnames(merged_df))
  subset_list[[1]] <- goa_shared_with_my_set
  subset_list[[2]] <- goa_shared_with_bulk_data
  
  fsub <- as.formula(paste("surv_object_se2 ~ ",
                         paste(goa_shared_with_bulk_data[1:length(goa_shared_with_bulk_data)], collapse= "+")))
  
  
  res.cox <- coxph(fsub, data=merged_df)
  summary_of_res_cox <- summary(res.cox)
  subset_list[[3]] <- summary_of_res_cox
  
  return(subset_list)
}

apop <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Apoptosis")
cc <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Cell_Cycle")
cd <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Cell_Death")
cm <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Cell_Motility")
dr <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "DNA_Repair")
ir <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Immune_Response")
phos1 <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Phosphorylation_1")
phos2 <- goa_analyzer(my.file = "Gene-sets-from-Gao-paper.txt", my.column = "Phosphorylation_2")

#4. Watanabe gene set (ref 26) ----

#5. Kwon et al. ----
kwon_raw_data <- read.table(file = "Kwon-gene-set.txt", sep = "\t", header = TRUE)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://may2009.archive.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup_kwon <- getBM(
  mart=mart,
  attributes=c(
    #"refseq_mrna",
    #"refseq_mrna_predicted",
    "refseq_dna",
    "refseq_dna_predicted",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_id"),
    #"external_gene_name"),
  filter = "refseq_dna", 
  #filter = "refseq_mrna",
  values = kwon_raw_data$Gene_Symbol, uniqueRows=TRUE)

kwon_gene_symbols_old <- annotLookup_kwon$external_gene_name
kwon_gene_symbols_old <- unique(kwon_gene_symbols)
kwon_gene_symbols_old <- as.data.frame(kwon_gene_symbols)
colnames(kwon_gene_symbols_old) <- "Gene_Symbol"

#Back to normal ----
results_from_other_papers <- list()

coxph_batch <- function(my.file){
  subset_list <- list()
  my_subset <- read.table(file = my.file, sep = "\t", header = TRUE)
  #my_subset <- my.file
  shared_with_my_set <- intersect(my_subset$Gene_Symbol, genes_in_bulk_RNA_vector)
  shared_with_bulk_data <- intersect(my_subset$Gene_Symbol, colnames(merged_df))
  subset_list[[1]] <- shared_with_my_set
  subset_list[[2]] <- shared_with_bulk_data
  
  fsub <- as.formula(paste("surv_object_se2 ~ ",
                           paste(shared_with_bulk_data[1:length(shared_with_bulk_data)], collapse= "+")))

  
  res.cox <- coxph(fsub, data=merged_df)
  summary_of_res_cox <- summary(res.cox)
  subset_list[[3]] <- summary_of_res_cox
  
  return(subset_list)
  
}
results_adder <- function(my.items,my.list, list.length, my.index=1){
  for (my_index in list.length){
    my.list[[my_index]] <- my.items[my_index]
    my_index <- my_index + 1
  }
  return(my.list)
}

watanabe <- coxph_batch(my.file = "Watanabe-gene-set2.txt")
bertucci <- coxph_batch(my.file = "Bertucci-gene-set.txt")
arango <- coxph_batch(my.file = "Arango-gene-set.txt")
nguyen <- coxph_batch(my.file = "Nguyen-gene-set.txt")
lenehan <- coxph_batch(my.file = "Lenehan-gene-set.txt") #Can't run because no overlap with my bulk RNA-seq set
kwon <- coxph_batch(my.file = kwon_gene_symbols)

results_from_other_papers[["watanabe"]] <- watanabe
results_from_other_papers[["bertucci"]] <- bertucci
results_from_other_papers[["arango"]] <- arango
results_from_other_papers[["kwon"]] <- kwon
results_from_other_papers[["barrier"]] <- summary_of_barrier_sigs
results_from_other_papers[["barrier_tumor"]] <- summary_of_barrier_sigs_tumor
results_from_other_papers[["xu"]] <- summary_of_xu_sigs
results_from_other_papers[["goa_apop"]] <- apop
results_from_other_papers[["goa_cc"]] <- cc
results_from_other_papers[["goa_cd"]] <- cd
results_from_other_papers[["goa_cm"]] <- cm
results_from_other_papers[["goa_dr"]] <- dr
results_from_other_papers[["goa_ir"]] <- ir
results_from_other_papers[["goa_phos1"]] <- phos1
results_from_other_papers[["goa_phos2"]] <- phos2
