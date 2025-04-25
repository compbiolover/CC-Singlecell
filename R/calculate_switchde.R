#' Calculate Switch-like Differential Expression (SwitchDE) gene rankings
#'
#' This function builds Inference of switch-like differential expression along
#' single-cell RNA sequencing trajectories (switchde) outputs.
#' Based on the method described in Campbell & Yau (2017): switchde: inference of
#' switch-like differential expression along single-cell trajectories.
#'
#' @param denoised_sc A gene expression matrix with genes as rows and cells as columns
#' @param pseudo_time A numeric vector of pseudotime values for each cell
#' @param zero_inflated Logical indicating whether to use zero-inflated model (default: FALSE)
#' @param q_threshold Q-value threshold for filtering significant genes (default: 0.05)
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param n_cores Number of cores to use for parallel processing (default: NULL, uses available cores - 1)
#' @param normalize Logical indicating if the returned scores of genes are normalized? (default: TRUE)
#'
#' @return A named numeric vector of normalized SwitchDE k values sorted by absolute magnitude
#' @references Campbell KR, Yau C (2017). "switchde: inference of switch-like differential
#'   expression along single-cell trajectories." Bioinformatics, 33(8), 1241-1242.
#'   doi:10.1093/bioinformatics/btw798
#' @export
#'
#' @examples
#' \dontrun{
#' sde_genes <- switchde_calculator(denoised_sc = expression_matrix, pseudo_time = pseudotime_data, normalize = TRUE)
#' }
switchde_calculator <- function(
    denoised_sc,
    pseudo_time,
    zero_inflated = FALSE,
    q_threshold = 0.05,
    parallel = FALSE,
    n_cores = NULL,
    normalize = TRUE
) {
    # Input validation
    if (!is.matrix(denoised_sc)) {
        stop("denoised_sc must be a matrix")
    }

    if (length(denoised_sc) == 0) {
        stop("denoised_sc must not be empty")
    }

    if (!is.numeric(pseudo_time)) {
        # Check if pseudo_time is a dataframe with a Pseudotime column
        if (
            is.data.frame(pseudo_time) &&
                "Pseudotime" %in% colnames(pseudo_time)
        ) {
            pseudo_time <- as.numeric(pseudo_time$Pseudotime)
        } else {
            stop(
                "pseudo_time must be a numeric vector or a dataframe with a 'Pseudotime' column"
            )
        }
    }

    if (length(pseudo_time) != ncol(denoised_sc)) {
        stop(
            "Length of pseudo_time must match the number of cells (columns) in denoised_sc"
        )
    }

    # Check for required packages
    if (!requireNamespace("switchde", quietly = TRUE)) {
        stop("Package 'switchde' is required for this function")
    }

    if (!requireNamespace("tidyverse", quietly = TRUE)) {
        stop("Package 'tidyverse' is required for this function")
    }

    # Set up parallel processing if requested
    if (parallel) {
        if (
            !requireNamespace("future", quietly = TRUE) ||
                !requireNamespace("furrr", quietly = TRUE)
        ) {
            warning(
                "Packages 'future' and 'furrr' are required for parallel processing.
              Falling back to sequential processing."
            )
            parallel <- FALSE
        } else {
            # Determine number of cores to use
            if (is.null(n_cores)) {
                if (!requireNamespace("parallel", quietly = TRUE)) {
                    warning(
                        "Package 'parallel' is required for auto-detecting cores.
                  Using 2 cores."
                    )
                    n_cores <- 2
                } else {
                    n_cores <- max(1, parallel::detectCores() - 1)
                }
            }

            # Set up parallel backend
            future::plan(future::multicore, workers = n_cores)
        }
    }

    # Calculate switchde values
    if (parallel) {
        # Not implementing parallel version yet as it would require restructuring switchde internals
        warning(
            "Parallel processing for switchde is not yet implemented. Using sequential processing."
        )
    }

    # Run switchde analysis
    sde <- switchde::switchde(
        denoised_sc,
        pseudo_time,
        verbose = TRUE,
        zero_inflated = zero_inflated
    )

    # Filter by q-value threshold
    sde_filtered <- dplyr::filter(sde, qval < q_threshold)

    # If no genes pass the threshold, warn and return empty result
    if (nrow(sde_filtered) == 0) {
        warning("No genes passed the q-value threshold of ", q_threshold)
        return(numeric(0))
    }

    # Order by absolute k value
    index <- order(abs(sde_filtered$k), decreasing = TRUE)
    sde_rank <- sde_filtered[index, ]

    # Extract k values as a named vector
    sde_ranking <- sde_rank$k
    names(sde_ranking) <- sde_rank$gene

    # Normalize the values
    if (normalize) {
        sde_ranking <- abs(sde_ranking) / sum(abs(sde_ranking))
    } else {
        sde_ranking
    }

    return(sde_ranking)
}
