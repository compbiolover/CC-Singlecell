#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Code to test several scRNA-seq methods to my method

#Loading needed packages----
library(magrittr);packageVersion("magrittr")
library(MAST);packageVersion("MAST")


#Loading our single cell data----
all_tumor_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_tumor_all_cells_FPKM.csv")
all_tumor_cells_fpkm$Condition <- "Cancer"
colnames(all_tumor_cells_fpkm)[1] <- "Genes"
all_nm_cells_fpkm <- read.csv("Data/Single-cell-data/GSE81861_CRC_NM_all_cells_FPKM.csv")
all_nm_cells_fpkm$Condition <- "Normal_Mucosa"
colnames(all_nm_cells_fpkm)[1] <- "Genes"
combined_data <- all_nm_cells_fpkm %>% full_join(all_tumor_cells_fpkm, by="Genes")
#combined_data <- as.matrix(combined_data)




#Looking at tutorial data to understand how to use MAST----
scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)

#Trying to load my own data into MAST format----
my_fData <- as.data.frame(dimnames(all_tumor_cells_fpkm)[[1]])
my_cData <- as.data.frame(dimnames(all_tumor_cells_fpkm)[[2]])

my_data_raw <- FromMatrix(exprsArray   = all_tumor_cells_fpkm,
                          fData        = my_fData, 
                          cData        = my_cData,
                          check_sanity = FALSE)
