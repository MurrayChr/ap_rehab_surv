#' Functions to get expected frequencies for contingency tables, either
#' (i) under null hypothesis of no difference between the rows i.e. for tests 
#'     of homogeneity (3G.SR and WBWA), or
#' (ii) under the null hypothesis that each multinomial in one set is a mixture
#'      of the multinomials in a second set (M.ITEC)     

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Expected frequencies for tests of homogeneity ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Expected frequencies for tests of mixtures ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_expected_frequencies_mixtures <- function(obs_comp, obs_mix, cell_probs_comp, cell_probs_mix) {
  # test that everything is the size it should to be
  if (ncol(obs_comp) != ncol(obs_mix)) {
    stop("Matrices for component and mixture data should have the same number of columns.")
  }
  if (any(dim(obs_comp) != dim(cell_probs_comp))) {
    stop("Component data matrix and cell probability matrix should have the same dimensions.")
  }
  if (any(dim(obs_mix) != dim(cell_probs_mix))) {
    stop("Mixture data matrix and cell probability matrix should have the same dimensions.")
  }
  M <- ncol(obs_comp)
  N_comp <- nrow(obs_comp)
  N_mix <- nrow(obs_mix)
  row_totals_comp <- rowSums(obs_comp)
  col_totals_comp <- colSums(obs_comp)
  row_totals_mix <- rowSums(obs_mix)
  col_totals_mix <- colSums(obs_mix)

  # if obs_comp is all zeros return NA (is there a better option?)
  if (sum(col_totals_comp) == 0) {
    return(NA)
  }
  
  # if obs_mix is all zeros, return zero matrix
  if (sum(col_totals_mix) == 0) {
    exp_mix <- matrix(0, nrow = N_mix, ncol = M)
    rownames(exp_mix) <- rownames(obs_mix)
    colnames(exp_mix) <- colnames(obs_mix)
  }
  
  # if neither of above, calculated expected values
  if ((sum(col_totals_mix) > 0) & (sum(col_totals_comp) > 0)) { 
    # expected values for component data
    exp_comp <- row_totals_comp * cell_probs_comp
    colnames(exp_comp) <- colnames(obs_comp)
    rownames(exp_comp) <- rownames(obs_comp)
    # expected values for mixture data
    exp_mix <- row_totals_mix * cell_probs_mix
    colnames(exp_mix) <- colnames(obs_mix)
    rownames(exp_mix) <- rownames(obs_mix)
  }
  
  # return
  list(expected_comp = exp_comp, expected_mix = exp_mix)
}

