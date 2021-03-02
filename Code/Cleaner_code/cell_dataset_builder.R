#Name: cell_dataset_builder.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently build cell dataset objects
#and get pseudotime information. 

#cell_dataset_builder----
cell_dataset_builder <- function(vim.genes=vim,
                                 cell.data=denoised_sc,
                                 cell.meta=gene_metadata){
  #Required packages----
  require(monocle3)
  
  #Return list----
  cds_return_list <- list()
  
  #Making the CDS object----
  vim_genes <- vim.genes
  cell_data_set <- new_cell_data_set(cell.data, gene_metadata = cell.meta)
  cds <- preprocess_cds(cell_data_set, num_dim=100, method="PCA")
  cds <- reduce_dimension(cds)
  plot_cells(cds)
  cds <- cluster_cells(cds)
  plot_cells(cds, color_cells_by="partition")
  vim_plot <- plot_cells(cds, genes=vim_genes, cell_size=0.50)
  cds <- learn_graph(cds)
  cds <- order_cells(cds)
  pt_graph <-plot_cells(cds                 = cds,
             color_cells_by      = "pseudotime",
             label_cell_groups   = FALSE,
             label_leaves        = FALSE,
             label_branch_points = FALSE,
             graph_label_size    = 1.5)
  
  
  #Integrating the new pseudotime data into the count matrix.----
  #Here I extract just the pseudotime calculations from the cds
  #object to a use for the swithde() calculation later. 
  #The output of this code is a dataframe that has 2 columns (pseudotime
  #value and sample name[cell name]) and 375 rows
  #(cells with their associated pseudotime).
  pseudotime_data <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pseudotime_data <- as.data.frame(pseudotime_data)
  colnames(pseudotime_data) <- "Pseudotime"
  pseudotime_contents <- rownames(pseudotime_data)
  pseudotime_data$Samples <- pseudotime_contents
  
  #Returning our objects----
  cds_return_list[["CDS"]] <- cds
  cds_return_list[["Pseudotime"]] <- pseudotime_data
  cds_return_list[["VIM Graph"]] <- vim_plot
  cds_return_list[["PT Graph"]] <- pt_graph
  return(cds_return_list)
}

