#Name: mirna_calculator_new.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently build mirna metric
#outputs

#mirna_calculator----
mirna_calculator <- function(ts.org                      ="Human", 
                             ts.version                  ="7.2",
                             max.miR.targets             =10,
                             cancer.up                   =TRUE,
                             cancer.type1                ="colon cancer",
                             cancer.type2                ="colorectal cancer",
                             mirna.remove                ="hsa-miR-129-1-3p",
                             max.mirnas                  =1559,
                             save.mirna.genes            =TRUE,
                             mirna.gene.rfile            ="~/Desktop/my_mirnas.rds",
                             mirna.ranking               ="~/Desktop/mirna.ranking.rds"){
  
  #Loading required package----
  require(hoardeR)
  require(tidyverse)
  
  #Our return list 
  miRNA_returns <- vector(mode = "list", length = 4)
  
  #miRNAs from miRmap----
  miRmap_mirnas <- read.csv(file = "Data/miRNA-data/mirmap_mirnas.csv", sep = ',')
  
  if(cancer.up==TRUE){
    dbDEMC_high <- read.csv(file = "Data/miRNA-data/dbDEMC-2.0-high.txt", sep = '\t')
    
    
    
    #Filtering to just the miRNAs associated with a particular type of cancer. 
    dbDEMC_high <- filter(dbDEMC_high, Cancer.Type==cancer.type1 | Cancer.Type==cancer.type2)
    dbDEMC_high_miRNAs <- subset(dbDEMC_high, select = miRBase.Update.ID)
    if ("unknown" %in% dbDEMC_high_miRNAs$miRBase.Update.ID==TRUE){
      dbDEMC_high_miRNAs <- filter(dbDEMC_high_miRNAs, miRBase.Update.ID != "unknown")
      dbDEMC_high_miRNAs <- as.vector(dbDEMC_high_miRNAs)
    }else{
      dbDEMC_high_miRNAs <- as.vector(dbDEMC_high_miRNAs)
    }
    
    #Common miRNAs
    common_mirnas <- intersect(miRmap_mirnas$mature_name, dbDEMC_high_miRNAs$miRBase.Update.ID)
    
  }else{
    dbDEMC_low <- read.csv(file = "Data/miRNA-data/dbDEMC-2.0-low.txt", sep = '\t')
    colnames(dbDEMC_low)[1] <- "miRNA.ID"
    
    #Filtering to just the miRNAs associated with a particular type of cancer. 
    dbDEMC_low <- filter(dbDEMC_low, Cancer_Type==cancer.type1 | Cancer_Type==cancer.type2)
    dbDEMC_low_miRNAs <- subset(dbDEMC_low, select = miRBase_update)
    dbDEMC_low_miRNAs <- as.vector(dbDEMC_low_miRNAs)
    
    #Common miRNAs between databases
    common_mirnas <- intersect(miRmap_mirnas$mature_name, dbDEMC_low_miRNAs$miRBase_update)
  }
  
  #Now submitting these miRNAs to targetScan to get genes to make a gene list for the third metric----
  #testing to see if all of the miRNAs exist in targetScan before getting all submitted at once
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna.remove]
  if(length(common_mirnas)>= length(common_mirnas[1:max(max.mirnas)])){
    common_mirnas <- common_mirnas[1:max.mirnas]
  }else{
    print("There are fewer target miRNAs available than your input.\n Using the largest number of common miRNAs for this submission to TargetScan")
    common_mirnas <- common_mirnas[1:length(common_mirnas)]
    print(paste0("The number of common mirnas is:", print(length(common_mirnas))))
  }
  
  #Using the targetScan() function of the hoardeR package to get our gene
  #gene targets for our specific miRNAs
    miRNA_targets <- targetScan(mirna=common_mirnas,
                                species=ts.org,
                                release=ts.version,
                                maxOut= max.miR.targets)
  
  #Saving the original unaltered targetScan output
  miRNA_returns[[1]] <- miRNA_targets
    
  #Subsetting the returning output from targetScan  
  # miRNA_targets <- lapply(miRNA_targets, "[", c(1,3,4))
  
  #Our numeric conversion functions for lapply()
  numeric_converter1 <- function(df){
  within(df, consSites <- as.numeric(consSites))
  }

  numeric_converter2 <- function(df){
    within(df, poorlySites <- as.numeric(poorlySites))
  }


  #Turning our specific columns in all of our list's data frames to
  #numeric data types
  miRNA_targets <- lapply(miRNA_targets, numeric_converter1)
  miRNA_targets <- lapply(miRNA_targets, numeric_converter2)


  # #Using a loop to go through each data frame in our list of data frames and
  # #removing any values that are equal to NA and replacing them with 0 to 
  # #make sure our scoring calculation don't break. 
  miRNA_targets_mod <- list()
  counter <- 1

  for(x in miRNA_targets){
    current_df <- x
    current_df[is.na(current_df)] <- 0
    miRNA_targets_mod[[counter]] <- current_df
    counter <- counter + 1
  }



  #We create a column in all of our list's data frames that is the fraction
  #of well conserved binding sites over total binding sites and remove any
  #entries that lead to NaN. We then create another column in all of the list's
  #data frames that sums all of the binding sites (well conserved and poorly
  #conserved). Finally, we create a score column that sums all of the numeric
  #values in each of the list's data frames and add an additional plus 1 for
  #the fact that this gene is considered a target of this miRNA by targetScan().
  miRNA_targets_mod <- lapply(miRNA_targets_mod, function(x) {x$mirna_frac <- x$consSites/(x$consSites + x$poorlySites);return(x)})
  miRNA_targets_mod <- lapply(miRNA_targets_mod, na.omit)
  miRNA_targets_mod <- lapply(miRNA_targets_mod, function(x) {x$total_binding_sites <-(x$consSites + x$poorlySites);return(x)})
  miRNA_targets_mod <- lapply(miRNA_targets_mod, function(x) {x$score <- rowSums(x[3:6]);return(x)})
  miRNA_targets_mod <- lapply(miRNA_targets_mod, function(x) {x$score <- x$score +1;return(x)})

  #Getting the names of the miRNAs from the original set of lists for our
  #modified lists
  names(miRNA_targets_mod) <- names(miRNA_targets)


  #Adding the miRNA-gene target matrices to our return list
  miRNA_returns[[2]] <- miRNA_targets_mod

  #Creating our score matrix by binding all of the data frames in our
  #list together and then sorting them from greatest to least score
  miRNA_score <- dplyr::bind_rows(miRNA_targets_mod)
  miRNA_score <- miRNA_score[order(-miRNA_score$score), ]

  #Function for standardizing our score from 0 to 1
  range_standridizing <- function(x){(x-min(x))/(max(x)-min(x))}
  miRNA_score$score <- range_standridizing(miRNA_score$score)
  
  #Now removing duplicate gene rows
  miRNA_score <- miRNA_score[!duplicated(miRNA_score[c('Ortholog')]), ]

  #Adding our score matrix to our return list and returning all our objects
  miRNA_returns[[3]] <- miRNA_score
  
  #Taking our finished data frame and constructing a named numeric from it 
  mirna.ranking <- miRNA_score$score
  names(mirna.ranking) <- miRNA_score$Ortholog
  
  #Adding our named numeric to the return object
  miRNA_returns[[4]] <- mirna.ranking
  
  #Saving the ranked miRNAs to a .rds file for use later
  if(save.mirna.genes==TRUE){
    saveRDS(mirna.ranking, file = mirna.gene.rfile)
  }
  
  return(miRNA_returns)
  
}


