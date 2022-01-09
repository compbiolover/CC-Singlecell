#Name: magic_denoiser.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Denoise single-celll RNA seq data.
#The genes must be on the rows and the cells
#on the columns.

#magic_denoiser----
magic_denoiser <- function(sc.data      =sc, 
                           magic.seed   =123,
                           magic.knn    =10,
                           magic.solver ='approximate',
                           magic.cores  =1,
                           magic.verbose=1,
                           filt.1       =0.3, 
                           filt.2       =0, 
                           filt.3       =0.2){
  #Loading required packages----
  require(ggplot2)
  require(Rmagic)
  
  #A list to store our return objects in----
  magic_return_list <- list()
  
  #The actual denoising of the single-cell data----
  sc.data <- t(sc.data)
  denoised_sc <- magic(sc.data, seed = magic.seed, solver=magic.solver, knn=magic.knn, n.jobs=magic.cores, verbose=magic.verbose)
  denoised_sc <- t(denoised_sc[["result"]])
  denoised_sc <- denoised_sc[rowMeans(denoised_sc) > filt.1 & rowMeans(denoised_sc > filt.2) > filt.3,]
  denoised_sc <-  as.data.frame(denoised_sc)
  
  #Getting the gene names of the single-cell data and saving it----
  #as a data.frame for use in the cell data set object in Monocle3
  cds_gene_names <- rownames(denoised_sc)
  cds_gene_names <- as.data.frame(cds_gene_names)
  colnames(cds_gene_names)[1] <- "gene_short_name"
  rownames(cds_gene_names) <- cds_gene_names$gene_short_name
  denoised_sc <- as.matrix(denoised_sc)
  
  #Our return list----
  magic_return_list[["denoised_sc_dataframe"]] <- denoised_sc
  magic_return_list[["cds_gene_names"]] <- cds_gene_names
  return(magic_return_list)
}
