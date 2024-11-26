# Function to get point estimates for cell probabilities and mixture weights
# in tests of mixtures (M.ITEC)
# see Choquet et al U-CARE User Guide, section 5.5 for details
library(cmdstanr)
library(tidyverse)

get_cell_probs_and_mixture_weights <- function(M, N_mix, N_comp, Y_comp, Y_mix) {
  # check dimensions
  if (any(dim(Y_comp) != c(N_comp, M))) {
    stop("Dimensions of Y_comp declared do not equal the dimensions provided.")
  }
  if (any(dim(Y_mix) != c(N_mix, M))) {
    stop("Dimensions of Y_mix declared do not equal the dimensions provided.")
  }
  
  # get maximum-likelihood estimates 
  file <- "stan/00_multinomial_mixture.stan" 
  mod <- cmdstan_model(file)
  stan_data <- list(M = M, N_comp = N_comp, N_mix = N_mix, Y_comp = Y_comp, Y_mix = Y_mix)
  ml <- mod$optimize(stan_data)
  
  # format estimates in matrices 
  cell_probs <- ml$summary("p_comp") %>%
    mutate(
      row = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
      col = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") )
    ) %>%
    arrange(row,col) %>%
    pull(estimate) %>%
    matrix(nrow = N_comp, ncol = M, byrow = TRUE)
  
   mix_wts <- ml$summary("mix_wts") %>%
    mutate(
      row = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
      col = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") )
    ) %>%
    arrange(row,col) %>%
    pull(estimate) %>%
    matrix(nrow = N_mix, ncol = N_comp, byrow = TRUE)
  
  # output
  list(cell_probs = cell_probs, mix_wts = mix_wts)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              ---- Test the function on simulated data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::install_github("dkahle/dirichlet")
# library(dirichlet)
# 
# #### Generate data ----
# M <- 5
# N_comp <- 3
# N_mix <- 2
# 
# # cell probabilities component multinomials
# cell_probs <- rdirichlet(N_comp, rep(1,M))
# # rowSums(cell_probs)
# 
# # mixture weights
# mix_wts <- rdirichlet(N_mix, rep(1,N_comp))
# # rowSums(mix_wts)
# 
# # cell probabilities for previously released states as mixtures of
# # currently released states
# prob <- mix_wts %*% cell_probs
# # rowSums(mix_wts)
# 
# # totals for each row
# N_trials <- 1e5
# 
# # generate the table Y
# Y_comp <- matrix(NA, N_comp, M)
# Y_mix <- matrix(NA, N_mix, M)
# for (i in 1:N_mix) {
#   Y_mix[i,] <- rmultinom(1, N_trials, prob[i,]) 
# }
# for (i in 1:N_comp) {
#   Y_comp[i,] <- rmultinom(1, N_trials, cell_probs[i,])
# }
# 
# #### Estimate cell probs and mixture weights
# ml_estimates <- get_cell_probs_and_mixture_weights(M, N_mix, N_comp, Y_comp, Y_mix)
# 
# # compare cell probabilities for lower half of table
# round(ml_estimates$cell_probs, 3)
# round(cell_probs, 3)
# 
# # compare mixture weights
# round(ml_estimates$mix_wts, 3)
# round(mix_wts, 3)
