#' Get expected frequencies for 2x2 contingency tables under null hypothesis
#' of no difference between the rows (i.e. a test of homogeneity).

get_expected_frequencies <- function(observed) {
  # check dimensions
  if (!(all(dim(observed) == c(2,2)))) {
    stop("Contingency tables must be 2 x 2.")
  }
  row_totals <- matrix( rowSums(observed), nrow = 2, ncol = 1) 
  col_totals <- matrix( colSums(observed), nrow = 1, ncol = 2)
  
  # check for zero table 
  if (sum(col_totals) == 0) {
    expected <- matrix(0, nrow = 2, ncol = 2)
    rownames(expected) <- rownames(observed)
    colnames(expected) <- colnames(observed)
  }
  
  if (sum(col_totals) > 0) {
    expected <- row_totals %*% col_totals / sum(col_totals)
    rownames(expected) <- rownames(observed)
    colnames(expected) <- colnames(observed)
  }

  # return
  expected
}
