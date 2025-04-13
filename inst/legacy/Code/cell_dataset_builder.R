# Name: cell_dataset_builder.R----
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Efficiently build cell data set objects
# and get pseudotime information.

cell_dataset_builder <- function(vim.genes = vim,
                                 cell.data = denoised_sc,
                                 cell.meta = gene_metadata,
                                 norm.flag = c("log"),
                                 graph.root = 1.5,
                                 show.traj = TRUE,
                                 point.size = 2.5,
                                 plot.cells = TRUE,
                                 my.root = "Y_19",
                                 save.gene.expr = FALSE,
                                 my.monocle.graph = "~/Desktop/pseudotime_graph.svg",
                                 my.monocle.graph.genes = "~/Desktop/genes_graph.svg",
                                 my.moncole.plot.type = "svg",
                                 my.monocle.plot.dpi = 600,
                                 my.monocle.width = 7.9,
                                 my.monocle.height = 7.9,
                                 my.monocle.unit = "in",
                                 use.function = TRUE,
                                 my.cds.filename = "~/Desktop/cds_output.rds",
                                 my.pt.data.filename = "~/Desktop/pt_data.csv") {
  # Required packages----
  require(ggplot2)
  require(monocle3)

  # Return list----
  cds_return_list <- list()

  # a helper function to identify the root principal points:
  get_earliest_principal_node <- function(cds) {
    initial_df <- colData(cds)
    initial_df <- as.data.frame(initial_df)
    cell_ids <- which.min(initial_df$Size_Factor)

    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
    root_pr_nodes
  }


  # Making the CDS object----
  vim_genes <- vim.genes
  cell_data_set <- new_cell_data_set(cell.data, gene_metadata = cell.meta)
  cds <- preprocess_cds(cell_data_set, num_dim = 100, method = "PCA", norm_method = norm.flag)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds, random_seed = 123)
  partition_plot <- plot_cells(cds, color_cells_by = "partition")
  if (plot.cells == TRUE) {
    vim_plot <- plot_cells(cds, genes = vim.genes, cell_size = point.size)
    cds_return_list[["Cell Progression Graph"]] <- vim_plot
    
    # The gene expression graph
    ggsave(
      filename = my.monocle.graph.genes,
      plot = print(vim_plot, newpage = FALSE),
      device = my.moncole.plot.type, dpi = my.monocle.plot.dpi,
      width = my.monocle.width, height = my.monocle.height,
      units = my.monocle.unit
    )
    if(save.gene.expr==TRUE){
      return(vim_plot)
    }
    
  }
  cds <- learn_graph(cds)

  if (use.function == TRUE) {
    my_roots <- get_earliest_principal_node(cds)
    print(my_roots)
    cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
    pt_graph <- plot_cells(
      cds = cds,
      color_cells_by = "pseudotime",
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = TRUE,
      label_principal_points = TRUE,
      graph_label_size = point.size,
      show_trajectory_graph = show.traj,
      cell_size = point.size
    )
    
  } else {
    cds <- order_cells(cds, root_pr_nodes = my.root)
    pt_graph <- plot_cells(
      cds = cds,
      color_cells_by = "pseudotime",
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = TRUE,
      label_principal_points = TRUE,
      graph_label_size = point.size,
      show_trajectory_graph = show.traj,
      cell_size = point.size
    )
    
  }


  # Integrating the new pseudotime data into the count matrix.----
  # Here I extract just the pseudotime calculations from the cds
  # object to a use for the swithde function calculation later.
  # The output of this code is a data frame that has 2 columns (pseudotime
  # value and sample name[cell name]) and 375 rows
  # (cells with their associated pseudotime).
  pseudotime_data <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pseudotime_data <- as.data.frame(pseudotime_data)
  colnames(pseudotime_data) <- "Pseudotime"
  pseudotime_contents <- rownames(pseudotime_data)
  pseudotime_data$Samples <- pseudotime_contents

  # Saving the  various plots to .svg or other specified format
  # with particular parameters

  # The pseudotime graph
  ggsave(
    filename = my.monocle.graph,
    plot = print(pt_graph, newpage = FALSE),
    device = my.moncole.plot.type, dpi = my.monocle.plot.dpi,
    width = my.monocle.width, height = my.monocle.height,
    units = my.monocle.unit
  )


  # Saving the relevant R objects and data frames
  saveRDS(cds, my.cds.filename)
  write.csv(pseudotime_data, my.pt.data.filename)

  # Returning our objects----
  cds_return_list[["CDS"]] <- cds
  cds_return_list[["Pseudotime"]] <- pseudotime_data
  cds_return_list[["PT Graph"]] <- pt_graph
  return(cds_return_list)
}
