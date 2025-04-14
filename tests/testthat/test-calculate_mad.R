# Create example data for testing
set.seed(123)
# Small test matrix with predictable values
test_matrix_small <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    nrow = 3,
    dimnames = list(c("gene1", "gene2", "gene3"), c("cell1", "cell2", "cell3"))
)
# Larger test matrix with more random values
test_matrix_large <- matrix(
    rnorm(100),
    nrow = 10,
    dimnames = list(paste0("gene", 1:10), paste0("cell", 1:10))
)

test_that("calculate_mad produces correct output format", {
    # Test with matrix input
    result <- calculate_mad(test_matrix_small)
    # Check that result is a named numeric vector
    expect_type(result, "double")
    expect_true(is.vector(result))
    expect_true(!is.null(names(result)))
    # Check that the length matches the number of rows in input
    expect_equal(length(result), nrow(test_matrix_small))
    # Check that the result is sorted in decreasing order
    expect_true(all(diff(result) <= 0))
    # Check that all values are between 0 and 1
    expect_true(all(result >= 0 & result <= 1))
    # Check that values sum to 1
    expect_equal(sum(result), 1)
})
test_that("calculate_mad performs correct calculation", {
    # Manually calculate expected result for small test matrix
    manual_mads <- c(
        gene1 = stats::mad(c(1, 4, 7)),
        gene2 = stats::mad(c(2, 5, 8)),
        gene3 = stats::mad(c(3, 6, 9))
    )
    # Sort in decreasing order
    manual_mads <- manual_mads[order(manual_mads, decreasing = TRUE)]
    # Normalize
    manual_mads_norm <- abs(manual_mads) / sum(abs(manual_mads))
    # Calculate result using our function
    result <- calculate_mad(test_matrix_small)
    # Compare results
    expect_equal(result, manual_mads_norm)
})
test_that("calculate_mad handles edge cases correctly", {
    # Empty matrix
    empty_matrix <- matrix(0, nrow = 0, ncol = 0)
    expect_error(calculate_mad(empty_matrix))
    # Single row
    single_row <- matrix(1:5,
        nrow = 1,
        dimnames = list("gene1", paste0("cell", 1:5))
    )
    single_result <- unname(calculate_mad(single_row))
    expect_length(single_result, 1)
    expect_equal(single_result, 1) # The only MAD value normalized
    # Constant values in a row (MAD = 0)
    constant_matrix <- matrix(rep(5, 15),
        nrow = 3, ncol = 5,
        dimnames = list(paste0("gene", 1:3), paste0("cell", 1:5))
    )
    const_result <- calculate_mad(constant_matrix)
    # All MADs are 0, so normalization should give NaN or equal values
    expect_true(all(is.nan(const_result)) || all(const_result == 1 / 3))
})
test_that("calculate_mad rejects invalid input", {
    # Non-matrix/data.frame input
    expect_error(calculate_mad(c(1, 2, 3)))
    expect_error(calculate_mad(list(a = 1, b = 2)))
    # NULL input
    expect_error(calculate_mad(NULL))
})
# Only run parallel tests if future and furrr are available
test_that("calculate_mad works in parallel mode", {
    skip_if_not_installed("future")
    skip_if_not_installed("furrr")
    # Compare sequential and parallel results
    result_seq <- calculate_mad(test_matrix_large)
    result_par <- calculate_mad(test_matrix_large, parallel = TRUE, n_cores = 2)
    # Results should be identical
    expect_equal(result_seq, result_par)
})
