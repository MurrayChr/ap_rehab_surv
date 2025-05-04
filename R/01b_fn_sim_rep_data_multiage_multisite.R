#' Simulate cmr data under the multi-age, multi-site model '01_multiage_multisite.stan'
#' (that does not incorporate trap-dependence or transience).
#' The code is factored so that data can be simulated under the model either
#'  (i) using parameter values from a posterior draw of the fit model and the 
#'  structure of the dataset to which it was fit (for posterior predictive checks),
#'  or
#' (ii) user-specified parameter values and compatible dataset structure (useful
#' when simulation-testing the model) 
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    ---- Function to sim data under '01_multiage_multisite.stan' ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' The function arguments are:
#' 'pars' - a list of matrices of parameter values 
#' 'data_str' - a list of dataset structure attributes, including study length,
#'   total number of individuals, a vector of first capture occasions, and a 
#'   vector of codes at first capture (which encode age and site information)

sim_cmr_data <- function(pars, data_str) {
  
  # check that all parameters required are present
  pars_required <- c("phi_ad", "phi_jv", "p_im", "p_ad", "m_jv", "m_ad")
  if (!setequal(names(pars), pars_required)) {
    stop(str_c("Error! Required parameters missing or unknown parameters given: ",
               str_flatten_comma(pars_required), " are required."))
  }
  
  # decant dataset structure
  N <- data_str$N
  T <- data_str$T
  fc <- data_str$fc
  fc_code <- data_str$fc_code
  
  # remove individuals whose first capture is the last capture occasion
  if (any(fc==T)) {
    fc_code <- fc_code[fc!=T] # NB: must be done before subsetting fc!
    fc <- fc[fc!=T]
    N <- length(fc)
    print("Warning: removed individuals first captured on last capture occasion!")
  }
  
  # check that parameters have the correct dimensions
  if (any(dim(pars$phi_ad)!=c(3,T-1))) {
    stop(str_c("Error! Incorrect dimensions for 'phi_ad'; should be 3 by ",T-1,"."))
  }
  if (any(dim(pars$phi_jv)!=c(3,T-1))) {
    stop(str_c("Error! Incorrect dimensions for 'phi_jv'; should be 3 by ",T-1,"."))
  }
  if (any(dim(pars$p_im)!=c(3,T))) {
    stop(str_c("Error! Incorrect dimensions for 'p_im'; should be 3 by ",T,"."))
  }
  if (any(dim(pars$p_ad)!=c(3,T))) {
    stop(str_c("Error! Incorrect dimensions for 'p_ad'; should be 3 by ",T,"."))
  }
  if (any(dim(pars$m_jv)!=c(3,3))) {
    stop("Error! Incorrect dimensions for 'm_jv'; should be 3 by 3.")
  }
  if (any(dim(pars$m_ad)!=c(3,3))) {
    stop("Error! Incorrect dimensions for 'm_ad'; should be 3 by 3.")
  }
  
  # decant parameters
  phi_jv <- pars$phi_jv
  phi_ad <- pars$phi_ad
  p_im <- pars$p_im
  p_ad <- pars$p_ad
  m_jv <- pars$m_jv
  m_ad <- pars$m_ad
  
  # create transition matrices 
  trans <- list()
  for (t in 1:(T-1)) {
    trans[[t]] <- matrix(NA, nrow=10, ncol=10)
    for (I in 1:3) {
      first_i <- 3*(I-1)+1 
      last_i <- 3*I  
      for (J in 1:3) {
        first_j <- 3*(J-1)+1
        last_j <- 3*J
        trans[[t]][first_i:last_i, first_j:last_j] <- matrix(
          c( 
            0, phi_jv[I,t]*m_jv[I,J],                     0,
            0,                     0, phi_ad[I,t]*m_jv[I,J], 
            0,                     0, phi_ad[I,t]*m_ad[I,J]
          ),
          nrow=3, ncol=3, byrow=TRUE
        )
      }
    }
    # bottom row for dead state
    trans[[t]][10, 1:9] <- 0
    # enforce row-sum-to-one constraint
    trans[[t]][,10] = 1 - rowSums(trans[[t]][,1:9])
  }
  
  # create observation matrices (these work a bit differently than usual)
  obs <- list()
  for (t in 1:T) {
    obs[[t]] <- matrix(0, nrow=10, ncol=10)
    # define diagonal blocks (off-diagonals in top-left 9x9 submatrix are zero)
    for (I in 1:3) {
      first_i <- 3*(I-1)+1
      last_i <- 3*I
      first_j <- first_i
      last_j <- last_i
      obs[[t]][first_i:last_i, first_j:last_j] <- matrix(
        c(
          1,         0,         0,
          0, p_im[I,t],         0,
          0,         0, p_ad[I,t]
        ),
        nrow=3, ncol=3, byrow=TRUE
      )
    }
    # bottom row for dead state
    obs[[t]][10,1:9] <- 0
    # enforce row-sum-to-one constraint
    obs[[t]][,10] = 1 - rowSums(obs[[t]][,1:9])
  }
  
  # initialise z matrix of states
  z <- matrix(0, nrow = N, ncol = T)
  for (i in 1:N) {
    z[i, fc[i]] <- fc_code[i]
  }
  
  # simulate states
  for (i in 1:nrow(z)) {
    for (t in fc[i]:(T-1)) {
      z[i, t+1] <- sample(1:10, 1, prob=trans[[t]][z[i,t], ])
    }
  }
  
  # initialise y matrix of observations
  y <- matrix(0, nrow = N, ncol = T)
  for (i in 1:nrow(y)) {
    y[i, fc[i]] <- z[i,fc[i]]
    for (t in (fc[i]+1):T) {
      y[i,t] <- sample(1:10, 1, prob=obs[[t]][z[i,t], ]) # ERROR!!! WE SHOULDN'T SAMPLE y on fc...
    }
  }
  
  # return
  list(fc=fc, fc_code=fc_code, z=z, y=y)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         --- Function to format a posterior draw for sim_cmr_data() ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function arguments:
# 'fit' - a cmdstanr fit of model '01_multiage_multisite.stan'
# 'draw' - a positive integer (usually between 1 and 4000) 

get_pars_from_posterior <- function( fit, draw ) {
  # make sure posterior draw is available
  if (draw > nrow(fit$draws(format="df"))) {
    stop(str_c("Draw ", draw ," requested but only ",
               nrow(fit$draws(format="df"))," draws available!"))
  }
  
  # format posterior sample
  vars <- c("phi_jv", "phi_ad", "p_im", "p_ad", "m_jv", "m_ad")
  post_sample <- fit$draws(vars, format="df") %>%
    filter(.draw==draw) %>%
    select(-starts_with("."))
  
  # decant values to matrices 
  extract_matrix <- function(post_sample, par) {
    temp <- post_sample %>%
      select(starts_with(par)) %>%
      pivot_longer(everything()) %>%
      mutate(
        site = as.integer(str_extract(name, "(?<=\\[)[1-9](?=,)"))
      )
    mat <- lapply(1:3, function(i) filter(temp, site==i)$value) %>%
      do.call(rbind,.)
    mat
  }
  
  # package as list and return
  pars_list <- lapply( vars, function(par) extract_matrix(post_sample, par) )
  names(pars_list) <- vars
  pars_list
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Wrapper for convenience ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim_rep_cmr_data <- function(fit, draw, data_str) {
  pars <- get_pars_from_posterior(fit, draw)
  sim_cmr_data(pars, data_str)
}

