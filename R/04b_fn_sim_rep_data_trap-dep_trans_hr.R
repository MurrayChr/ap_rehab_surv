#' Simulate cmr data under the multi-age, multi-site model 
#' '04_td_tr_hr.stan' that incorporates trap-dependence and transience and a 
#' hand-rearing covariate on adult and juvenile survival.
#' The code is factored so that data can be simulated under the model either
#'  (i) using parameter values from a posterior draw of the fit model and the 
#'  structure of the dataset to which it was fit (for posterior predictive checks),
#'  or
#' (ii) user-specified parameter values and compatible dataset structure (useful
#' when simulation-testing the model) 
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ---- Function to sim data under '04_td_tr_hr.stan' ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' The function arguments are:
#' 'pars' - a list of matrices/scalars of parameter values 
#' 'data_str' - a list of dataset structure attributes, including study length,
#'   total number of individuals, a vector of first capture occasions, a 
#'   vector of codes at first capture (which encode age and site information),
#'   and a vector of hand-rearing status

sim_cmr_data <- function(pars, data_str) {
  
  # check that all parameters required are present
  pars_required <- c("phi_ad_wr", "phi_jv_wr", "hr_jv", "hr_ad", "pi_r", 
                     "p_im", "p_ad_A", "p_ad_U", "p_ad_N",
                     "m_jv", "m_ad")
  if (!setequal(names(pars), pars_required)) {
    stop(str_c("Error! Required parameters missing or unknown parameters given: ",
               str_flatten_comma(pars_required), " are required."))
  }
  
  # decant dataset structure
  N <- data_str$N
  T <- data_str$T
  fc <- data_str$fc
  fc_code <- data_str$fc_code
  hr <- data_str$hr
  
  # remove individuals whose first capture is the last capture occasion
  if (any(fc==T)) {
    fc_code <- fc_code[fc!=T] # NB: must be done before subsetting fc!
    fc <- fc[fc!=T]
    N <- length(fc)
    print("Warning: removed individuals first captured on last capture occasion!")
  }
  
  # check that parameters have the correct dimensions
  if (any(dim(pars$phi_ad_wr)!=c(3,T-1))) {
    stop(str_c("Error! Incorrect dimensions for 'phi_ad_wr'; should be 3 by ",T-1,"."))
  }
  if (any(dim(pars$phi_jv_wr)!=c(3,T-1))) {
    stop(str_c("Error! Incorrect dimensions for 'phi_jv_wr'; should be 3 by ",T-1,"."))
  }
  if (any(dim(pars$pi_r)!=c(3,T-1))) {
    stop(str_c("Error! Incorrect dimensions for 'pi_r'; should be 3 by ",T-1,"."))
  }
  if (any(dim(pars$p_im)!=c(3,T))) {
    stop(str_c("Error! Incorrect dimensions for 'p_im'; should be 3 by ",T,"."))
  }
  if (any(dim(pars$p_ad_A)!=c(3,T))) {
    stop(str_c("Error! Incorrect dimensions for 'p_ad_A'; should be 3 by ",T,"."))
  }
  if (any(dim(pars$p_ad_U)!=c(3,T))) {
    stop(str_c("Error! Incorrect dimensions for 'p_ad_U'; should be 3 by ",T,"."))
  }
  if (any(dim(pars$p_ad_N)!=c(3,T))) {
    stop(str_c("Error! Incorrect dimensions for 'p_ad_N'; should be 3 by ",T,"."))
  }
  if (any(dim(pars$m_jv)!=c(3,3))) {
    stop("Error! Incorrect dimensions for 'm_jv'; should be 3 by 3.")
  }
  if (any(dim(pars$m_ad)!=c(3,3))) {
    stop("Error! Incorrect dimensions for 'm_ad'; should be 3 by 3.")
  }
  
  # decant parameters
  phi_jv_wr <- pars$phi_jv_wr
  phi_ad_wr <- pars$phi_ad_wr
  hr_ad <- pars$hr_ad
  hr_jv <- pars$hr_jv
  pi_r <- pars$pi_r
  p_im <- pars$p_im
  p_ad_A <- pars$p_ad_A
  p_ad_U <- pars$p_ad_U
  p_ad_N <- pars$p_ad_N
  m_jv <- pars$m_jv
  m_ad <- pars$m_ad
  
  # calculate juvenile and adult survival for hand-reared birds
  logit <- function(p) {log(p / (1 - p))}  # check for 0 < p < 1 first ?
  inv_logit <- function(x) {1 / (1 + exp(-x))}
  phi_ad_hr <- inv_logit(logit(phi_ad_wr) + hr_ad)
  phi_jv_hr <- inv_logit(logit(phi_jv_wr) + hr_jv)
  
  # create transition matrices for wild-raised birds
  trans_wr <- list()
  for (t in 1:(T-1)) {
    trans_wr[[t]] <- matrix(NA, nrow=13, ncol=13, byrow=TRUE)
    for (I in 1:3) {
      first_i <- 4*(I-1)+1 
      last_i <- 4*I  
      for (J in 1:3) {
        first_j <- 4*(J-1)+1
        last_j <- 4*J
        trans_wr[[t]][first_i:last_i, first_j:last_j] <- matrix(
          c( 
            0, phi_jv_wr[I,t]*m_jv[I,J],                                   0,                                       0,
            0,                     0, phi_ad_wr[I,t]*m_jv[I,J]*p_ad_N[J,t+1], phi_ad_wr[I,t]*m_jv[I,J]*(1-p_ad_N[J,t+1]),
            0,                     0, phi_ad_wr[I,t]*m_ad[I,J]*p_ad_A[J,t+1], phi_ad_wr[I,t]*m_ad[I,J]*(1-p_ad_A[J,t+1]),
            0,                     0, phi_ad_wr[I,t]*m_ad[I,J]*p_ad_U[J,t+1], phi_ad_wr[I,t]*m_ad[I,J]*(1-p_ad_U[J,t+1])
          ),
          nrow=4, ncol=4, byrow=TRUE
        )
      }
    }
    # bottom row for dead state
    trans_wr[[t]][13, 1:12] <- 0
    # enforce row-sum-to-one constraint
    trans_wr[[t]][,13] = 1 - rowSums(trans_wr[[t]][,1:12])
  }
  
  # create transition matrices for hand-reared birds
  trans_hr <- list()
  for (t in 1:(T-1)) {
    trans_hr[[t]] <- matrix(NA, nrow=13, ncol=13, byrow=TRUE)
    for (I in 1:3) {
      first_i <- 4*(I-1)+1 
      last_i <- 4*I  
      for (J in 1:3) {
        first_j <- 4*(J-1)+1
        last_j <- 4*J
        trans_hr[[t]][first_i:last_i, first_j:last_j] <- matrix(
          c( 
            0, phi_jv_hr[I,t]*m_jv[I,J],                                   0,                                       0,
            0,                     0, phi_ad_hr[I,t]*m_jv[I,J]*p_ad_N[J,t+1], phi_ad_hr[I,t]*m_jv[I,J]*(1-p_ad_N[J,t+1]),
            0,                     0, phi_ad_hr[I,t]*m_ad[I,J]*p_ad_A[J,t+1], phi_ad_hr[I,t]*m_ad[I,J]*(1-p_ad_A[J,t+1]),
            0,                     0, phi_ad_hr[I,t]*m_ad[I,J]*p_ad_U[J,t+1], phi_ad_hr[I,t]*m_ad[I,J]*(1-p_ad_U[J,t+1])
          ),
          nrow=4, ncol=4, byrow=TRUE
        )
      }
    }
    # bottom row for dead state
    trans_hr[[t]][13, 1:12] <- 0
    # enforce row-sum-to-one constraint
    trans_hr[[t]][,13] = 1 - rowSums(trans_hr[[t]][,1:12])
  }
  
  # create observation matrices (these work a bit differently than usual)
  obs <- list()
  for (t in 1:T) {
    obs[[t]] <- matrix(0, nrow=13, ncol=13)
    # define diagonal blocks (off-diagonals in top-left 12x12 submatrix are zero)
    for (I in 1:3) {
      first_i <- 4*(I-1)+1
      last_i <- 4*I
      first_j <- first_i
      last_j <- last_i
      obs[[t]][first_i:last_i, first_j:last_j] <- matrix(
        c(
          1,         0, 0, 0,
          0, p_im[I,t], 0, 0,
          0,         0, 1, 0,
          0,         0, 0, 0
        ),
        nrow=4, ncol=4, byrow=TRUE
      )
    }
    # bottom row for dead state
    obs[[t]][13,1:12] <- 0
    # enforce row-sum-to-one constraint
    obs[[t]][,13] = 1 - rowSums(obs[[t]][,1:12])
  }
  
  # initialise z matrix of states
  z <- matrix(0, nrow = N, ncol = T)
  for (t in 1:(T-1)) {
    z[fc==t, t] <- fc_code[fc==t]
  }
  
  # simulate states
  for (i in 1:nrow(z)) {
    if (fc_code[i] %in% c(3, 7, 11)) {  # individuals marked as adults
      # sample if resident or not
      site_ind <- as.integer((fc_code[i] - 3)/4 + 1) # maps 3,7,11 to 1,2,3
      res <- runif(1) < pi_r[site_ind, fc[i]]
      # if transient then 'dies'
      if (!res) {
        z[i, (fc[i]+1):T] <- 13
      }
      # if resident then fate is sampled from transition matrices with 
      # survival prob corresponding to hand-reared/wild-raised status
      if (res) {
        if (hr[i] == 1) {
          for (t in fc[i]:(T-1)) {
            z[i, t+1] <- sample(1:13, 1, prob=trans_hr[[t]][z[i,t], ])
          }
        } else {
          for (t in fc[i]:(T-1)) {
            z[i, t+1] <- sample(1:13, 1, prob=trans_wr[[t]][z[i,t], ])
          }
        }
      }
    } else {  # individuals not marked as adults (i.e. marked as juveniles)  
      if (hr[i] == 1) {
        for (t in fc[i]:(T-1)) {
          z[i, t+1] <- sample(1:13, 1, prob=trans_hr[[t]][z[i,t], ])
        }
      } else {
        for (t in fc[i]:(T-1)) {
          z[i, t+1] <- sample(1:13, 1, prob=trans_wr[[t]][z[i,t], ])
        }
      }
    } 
  }
  
  # initialise y matrix of observations
  y <- matrix(0, nrow = N, ncol = T)
  for (i in 1:nrow(y)) {
    y[i, fc[i]] <- z[i, fc[i]]   
    for (t in (fc[i]+1):T) {
      y[i,t] <- sample(1:13, 1, prob=obs[[t]][z[i,t], ])
    }
  }
  
  # return
  truth <- list(
    "phi_jv_wr" = phi_jv_wr,
    "phi_ad_wr" = phi_ad_wr,
    "phi_jv_hr" = phi_jv_hr,
    "phi_ad_hr" = phi_ad_hr,
    "hr_ad" = hr_ad,
    "hr_jv" = hr_jv,
    "pi_r" = pi_r,
    "p_im" = p_im,
    "p_ad_A" = p_ad_A,
    "p_ad_U" = p_ad_U,
    "p_ad_N" = p_ad_N,
    "m_jv" = m_jv,
    "m_ad" = m_ad
  )
  return(list(truth=truth, fc=fc, fc_code=fc_code, hr=hr, z=z, y=y))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         --- Function to format a posterior draw for sim_cmr_data() ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function arguments:
# 'fit' - a cmdstanr fit of model '04_td_tr_hr.stan'
# 'draw' - a positive integer (usually between 1 and 4000) 

get_pars_from_posterior <- function( fit, draw ) {
  # make sure posterior draw is available
  if (draw > nrow(fit$draws(format="df"))) {
    stop(str_c("Draw ", draw ," requested but only ",
               nrow(fit$draws(format="df"))," draws available!"))
  }
  
  vector_vars <- c("phi_jv_wr", "phi_ad_wr", "pi_r", 
                   "p_im", "p_ad_A", "p_ad_U", "p_ad_N", 
                   "m_jv", "m_ad")
  scalar_vars <- c("hr_jv", "hr_ad")
  vars <- c(scalar_vars, vector_vars)
  
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
  vector_pars_list <- lapply(vector_vars, function(par) extract_matrix(post_sample, par) )
  names(vector_pars_list) <- vector_vars
  scalar_pars_list <- lapply(scalar_vars, function(par) post_sample[[par]])
  names(scalar_pars_list) <- scalar_vars
  pars_list <- c(vector_pars_list, scalar_pars_list)
  
  return(pars_list)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Wrapper for convenience ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim_rep_cmr_data <- function(fit, draw, data_str) {
  pars <- get_pars_from_posterior(fit, draw)
  sim_cmr_data(pars, data_str)
}
