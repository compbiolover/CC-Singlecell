#' Calculate Median Absolute Deviation (MAD) metric gene lists
#'
#' This function calculates the median absolute deviation (MAD) for each gene
#' across cells in single-cell RNA sequencing data.
#'
#' @param expression_matrix A gene expression matrix with genes as rows and cells as columns
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param n_cores Number of cores to use for parallel processing (default: NULL, uses available cores - 1)
#' @param normalize Whether to normalize the scores of our ranked list of genes (default: TRUE)
#'
#' @return A named numeric vector of normalized MAD values sorted in decreasing order
#' @export
#'
#' @examples
#' \dontrun{
#' mad_genes <- calculate_mad(expression_matrix)
#' }
calculate_mad <- function(expression_matrix, parallel = FALSE, n_cores = NULL, normalize = TRUE) {
    # Input validation
    if (!is.matrix(expression_matrix)) {
        stop("expression_matrix must be a matrix")
    }

    if (length(expression_matrix) == 0) {
        stop("expression_matrix must not be empty")
    }

    # Set up parallel processing if requested
    if (parallel) {
        if (!requireNamespace("future", quietly = TRUE) ||
            !requireNamespace("furrr", quietly = TRUE)) {
            warning("Packages 'future' and 'furrr' are required for parallel processing.
              Falling back to sequential processing.")
            parallel <- FALSE
        } else {
            # Determine number of cores to use
            if (is.null(n_cores)) {
                n_cores <- max(1, parallel::detectCores() - 1)
            }

            # Set up parallel backend
            future::plan(future::multicore, workers = n_cores)
        }
    }

    # Calculate MAD values
    if (parallel) {
        # Parallel implementation using furrr
        gene_mads <- furrr::future_map_dbl(
            seq_len(nrow(expression_matrix)),
            function(i) stats::mad(expression_matrix[i, ]),
            .options = furrr::furrr_options(seed = TRUE)
        )
        names(gene_mads) <- rownames(expression_matrix)
    } else {
        # Sequential implementation
        gene_mads <- purrr::map_dbl(
            seq_len(nrow(expression_matrix)),
            function(i) stats::mad(expression_matrix[i, ])
        )
        names(gene_mads) <- rownames(expression_matrix)
    }

    # Rank and normalize
    if (normalize == TRUE) {
        ranked_mads <- gene_mads %>%
            tibble::enframe(name = "gene", value = "mad") %>%
            dplyr::arrange(dplyr::desc(mad)) %>%
            dplyr::mutate(mad_normalized = abs(mad) / sum(abs(mad)))

        # Extract the normalized MAD values as a named vector
        mad_ranking <- ranked_mads$mad_normalized
        names(mad_ranking) <- ranked_mads$gene
    } else {
        ranked_mads <- gene_mads %>%
            tibble::enframe(name = "gene", value = "mad") %>%
            dplyr::arrange(dplyr::desc(mad)) %>%
            dplyr::mutate(mad_unnormalized = mad)

        # Extract the un-normalized MAD values as a named vector
        mad_ranking <- ranked_mads$mad_unnormalized
        names(mad_ranking) <- ranked_mads$gene

        return(mad_ranking)
    }
}
