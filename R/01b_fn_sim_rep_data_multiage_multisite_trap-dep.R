# Function to simulate a replicate dataset under model under the multi-age,
# multi-site model '01_multiage_multisite_trap-dep.stan' that incorporates trap-
# dependence
library(tidyverse)

#' The function arguments are:
#' 'fit' - a cmdstanr model fit of '01_multiage_multisite_trap-dep.stan' to
#'   the data 'data/00b_cmr_data_multisite_multiage_trapdep.RDS'
#' 'draw' - a posterior sample to use for the data-generating parameter values
#' 'data_str' - a list of dataset structure attributes, including study length,
#'   total number of individuals, and a year x age x site array whose (i,j,k)
#'   entry is the number of individuals first captured in year i, of age-class j
#'   at site k

sim_rep_cmr_data <- function( fit, draw, data_str ) {
  
  # make sure posterior draw is available
  if (draw > nrow(fit$draws(format="df"))) {
    stop(str_c("Draw ", draw ," requested but only ",
               nrow(fit$draws(format="df"))," draws available!"))
  }
  
  # unpack dataset structure
  T <-data_str$T
  N <- data_str$N
  fc <- data_str$fc
  fc_code <- data_str$fc_code
  
  # remove individuals whose first capture is the last capture occasion
  if (any(fc==T)) {
    fc_code <- fc_code[fc!=T] # NB: must be done before subsetting fc!
    fc <- fc[fc!=T]
    N <- length(fc)
    print("Warning: removed individuals first captured on last capture occasion!")
  }
  
  # format posterior sample
  vars <- c("phi_jv", "phi_ad", "p_im", "p_ad_A", "p_ad_U", "p_ad_N", "m_jv", "m_ad")
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
  
  # parking this simpler but possibly less safe version of extract_matrix() here
  # extract_matrix <- function(post_sample, par) {
  #   post_sample %>%
  #     select(starts_with(par)) %>%
  #     pivot_longer( everything() ) %>%
  #     pull(value) %>%
  #     matrix(nrow=3)
  # }
  
  # site x year matrices for survival and detection
  phi_jv <- extract_matrix(post_sample, "phi_jv")
  phi_ad <- extract_matrix(post_sample, "phi_ad")
  p_im <- extract_matrix(post_sample, "p_im")
  p_ad_A <- extract_matrix(post_sample, "p_ad_A")
  p_ad_U <- extract_matrix(post_sample, "p_ad_U")
  p_ad_N <- extract_matrix(post_sample, "p_ad_N")
  #  source site x destinatino site matrices for movement
  m_jv <- extract_matrix(post_sample, "m_jv")
  m_ad <- extract_matrix(post_sample, "m_ad")
  
  # create transition matrices 
  trans <- list()
  for (t in 1:(T-1)) {
    trans[[t]] <- matrix(NA, nrow=13, ncol=13, byrow=TRUE)
    for (I in 1:3) {
      first_i <- 4*(I-1)+1 
      last_i <- 4*I  
      for (J in 1:3) {
        first_j <- 4*(J-1)+1
        last_j <- 4*J
        trans[[t]][first_i:last_i, first_j:last_j] <- matrix(
          c( 
            0, phi_jv[I,t]*m_jv[I,J],                                   0,                                       0,
            0,                     0, phi_ad[I,t]*m_jv[I,J]*p_ad_N[J,t+1], phi_ad[I,t]*m_jv[I,J]*(1-p_ad_N[J,t+1]),
            0,                     0, phi_ad[I,t]*m_ad[I,J]*p_ad_A[J,t+1], phi_ad[I,t]*m_ad[I,J]*(1-p_ad_A[J,t+1]),
            0,                     0, phi_ad[I,t]*m_ad[I,J]*p_ad_U[J,t+1], phi_ad[I,t]*m_ad[I,J]*(1-p_ad_U[J,t+1])
          ),
          nrow=4, ncol=4, byrow=TRUE
        )
      }
    }
    # bottom row for dead state
    trans[[t]][13, 1:12] <- 0
    # enforce row-sum-to-one constraint
    trans[[t]][,13] = 1 - rowSums(trans[[t]][,1:12])
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
    for (t in fc[i]:(T-1)) {
      z[i, t+1] <- sample(1:13, 1, prob=trans[[t]][z[i,t], ])
    }
  }
  
  # initialise y matrix of observations
  y <- matrix(0, nrow = N, ncol = T)
  for (i in 1:nrow(y)) {
    for (t in fc[i]:T) {
      y[i,t] <- sample(1:13, 1, prob=obs[[t]][z[i,t], ])
    }
  }
  
  # return
  list(fc=fc, fc_code=fc_code, z=z, y=y)
}
  
