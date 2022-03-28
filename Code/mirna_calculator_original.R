#Name: mirna_calculator.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently build miRNA metric outputs

#mirna_calculator----
mirna_calculator <- function(ts.org                      ="Human", 
                             ts.version                  ="7.2",
                             max.miR.targets             =10,
                             cancer.up                   =TRUE,
                             cancer.type1                ="colon cancer",
                             cancer.type2                ="rectal cancer",
                             print.ts.targets            =TRUE,
                             mirna.remove                ="hsa-miR-129-1-3p",
                             max.mirnas                  =1559,
                             save.mirna.genes            =TRUE,
                             mirna.gene.rfile            ="~/Desktop/my_mirnas.rds"){
  
  #Loading required package----
  require(hoardeR)
  require(tidyverse)
  
  
  #miRNAs from miRmap----
  miRmap_mirnas <- read.csv(file = "Data/miRNA-data/mirmap_mirnas.csv",
                            sep = ',')
  
  if(cancer.up==TRUE){
    dbDEMC_high <- read.csv(file = "Data/miRNA-data/dbDEMC-2.0-high.txt",
                            sep = '\t')
    
    #Filtering to just the miRNAs associated with a particular type of cancer. 
    dbDEMC_high <- filter(dbDEMC_high, 
                          Cancer.Type==cancer.type1 | Cancer.Type==cancer.type2)
    dbDEMC_high_miRNAs <- subset(dbDEMC_high, select = miRBase.Update.ID)
    if ("unknown" %in% dbDEMC_high_miRNAs$miRBase.Update.ID==TRUE){
      dbDEMC_high_miRNAs <- filter(dbDEMC_high_miRNAs,
                                   miRBase.Update.ID != "unknown")
      dbDEMC_high_miRNAs <- as.vector(dbDEMC_high_miRNAs)
    }else{
      dbDEMC_high_miRNAs <- as.vector(dbDEMC_high_miRNAs)
    }
    
    #Common mirnas
    common_mirnas <- intersect(miRmap_mirnas$mature_name,
                               dbDEMC_high_miRNAs$miRBase.Update.ID)
  }else{
    dbDEMC_low <- read.csv(file = "Data/miRNA-data/dbDEMC-2.0-low.txt",
                           sep = '\t')
    colnames(dbDEMC_low)[1] <- "miRNA.ID"
    
    #Filtering to just the miRNAs associated with a particular type of cancer. 
    dbDEMC_low <- filter(dbDEMC_low, 
                         Cancer_Type==cancer.type1 | Cancer_Type==cancer.type2)
    dbDEMC_low_miRNAs <- subset(dbDEMC_low, select = miRBase_update)
    dbDEMC_low_miRNAs <- as.vector(dbDEMC_low_miRNAs)
    
    #Common miRNAs between databases
    common_mirnas <- intersect(miRmap_mirnas$mature_name,
                               dbDEMC_low_miRNAs$miRBase_update)
  }
  
  #Now submitting these miRNAs to TargetScan to get genes to make a gene list
  #for the third metric testing to see if all of the miRNAs exist in TargetScan
  #before getting all submitted at once
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna.remove]
  if(length(common_mirnas)>= length(common_mirnas[1:max.mirnas])){
    common_mirnas <- common_mirnas[1:max.mirnas]
  }else{
    print("There are fewer target mirnas available than your input.
          Using the largest number of common mirnas for this submission
          to TargetScan")
    common_mirnas <- common_mirnas[1:length(common_mirnas)]
    print("The number of common mirnas is")
    print(length(common_mirnas))
  }
  
  if(print.ts.targets==TRUE){
    print(length(common_mirnas))
  }
  my_num <- 1
  miRNA_targets <- vector(mode = "list", length = length(common_mirnas))
  for (m in common_mirnas[1:length(common_mirnas)]) {
    percent_done <- (my_num/length(common_mirnas))*100
    #print(percent_done)
    #print(m)
    current_target <- targetScan(mirna=common_mirnas[my_num],
                                 species=ts.org,
                                 release=ts.version,
                                 maxOut= max.miR.targets)
    miRNA_name <- m
    miRNA_name_final <- rep(miRNA_name, times=length(current_target$Ortholog))
    current_target <- cbind(current_target,miRNA_name_final)
    miRNA_targets[[m]] <- current_target
    #print(my_num)
    my_num <- my_num + 1
  }
  
  
  
  #Simplifying the output of the TargetScan commands to 
  #just the Gene name and the miRNA columns
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
  
  miRNA_score <- matrix(data = 0, nrow = length(all_genes_for_score_unique),
                        ncol = length(all_miRs_for_score),
                        dimnames = list(all_genes_for_score_unique,
                                        all_miRs_for_score))
  miRNA_score <- as.data.frame(miRNA_score)
  
  #Checking to see for each miRNA (column name) if it interacts with a particular row.
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
  
  
  #Now calculating the row sums of each gene for total number of
  #the miRNA interactions----
  mirna_gene_list <- rowSums(miRNA_score)
  mirna.ranking<-abs(mirna_gene_list)/sum(abs(mirna_gene_list))
  mirna.ranking <- sort(mirna.ranking, decreasing = TRUE)
  if(save.mirna.genes==TRUE){
    saveRDS(mirna.ranking, file = mirna.gene.rfile)
  }
  #Return object----
  return(mirna.ranking)
}
