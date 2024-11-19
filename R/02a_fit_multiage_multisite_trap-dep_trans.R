# Fit a multi-age, multi-site model that incorporates trap-dependence and
# transience (but without a hand-reading covariate).
library(tidyverse)
library(cmdstanr)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            ---- Load cmr data and format into m-array ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# add last capture occasion to find individuals not recaptured after first capture
cmr_data <- cmr_data %>%
  rowwise() %>%
  mutate( 
    lc = max( which( c_across( starts_with("yr") )!=0 ) ),
    .after = fc
  ) %>%
  ungroup()

# split data 
cmr_data_not_for_marr <- cmr_data %>%
  filter(
    (marking_age == 3) & (fc == lc)
  )
cmr_data_for_marr <- cmr_data %>%
  filter(
    !((marking_age == 3) & (fc == lc))
  )

# construct m-array
y <- cmr_data_for_marr %>%
  select( starts_with("yr") ) %>%
  as.matrix()
fc <- cmr_data_for_marr$fc

# change '0's after first capture to nStates + 1 = 13 for get_marray()
T <- ncol(y)
for (i in 1:nrow(y)) {
  y[i,fc[i]:T][y[i, fc[i]:T]==0] <- 13 
}

# m-array for model '02_multiage_multisite_trap-dep_trans.stan' which has 12 states
# states 4,8,12 are unobservable so the warning about their absence can be safely ignored
marr <- get_marray(y, nStates=12)

# count number of birds marked as adults that are ever/never recaptured
# per site and cohort
N_1 <- matrix(NA, 3, T-1)
N_0 <- matrix(NA, 3, T-1)
for (i in 1:3) {
  for (t in 1:(T-1)) {
    N_1[i,t] <- cmr_data_for_marr %>%
      filter( fc_code==(3+4*(i-1)) & fc==t ) %>%
      nrow()
    N_0[i,t] <- cmr_data_not_for_marr %>%
      filter( fc_code==(3+4*(i-1)) & fc==t ) %>%
      nrow()
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ---- Fit the model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compile and fit model
file <- "stan/02_multiage_multisite_trap-dep_trans.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, marr=marr, N_1=N_1, N_0=N_0)
fit <- mod$sample(stan_data, parallel_chains = 4) 
# fit$save_object("outputs/02a_multiage_multisite_trap-dep_trans_fit.RDS")

# sampler diagnostics
fit$diagnostic_summary()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        ---- Plot estimates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~