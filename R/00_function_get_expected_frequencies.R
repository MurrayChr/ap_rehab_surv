#' Get expected frequencies for contingency tables under null hypothesis
#' of no difference between the rows (i.e. a test of homogeneity).

get_expected_frequencies <- function(observed) {
  # size of table
  m <- nrow(observed)
  n <- ncol(observed)
  
  row_totals <- matrix( rowSums(observed), nrow = m, ncol = 1) 
  col_totals <- matrix( colSums(observed), nrow = 1, ncol = n)
  
  # if observed is all zeros, set expected the same 
  if (sum(col_totals) == 0) {
    expected <- matrix(0, nrow = m, ncol = n)
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
