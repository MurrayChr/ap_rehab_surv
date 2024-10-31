# Fit a multi-age, multi-site model that incorporates trap-dependence but not
# transience (and without a hand-reading covariate).

library(tidyverse)
library(cmdstanr)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            ---- Load cmr data and format into m-array ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# extract capture history matrix
y <- cmr_data %>%
  select( starts_with("yr") ) %>%
  as.matrix()
fc <- cmr_data$fc

# change '0's after first capture to nStates + 1 = 13 for get_marray()
T <- ncol(y)
for (i in 1:nrow(y)) {
  y[i,fc[i]:T][y[i, fc[i]:T]==0] <- 13 
}

# m-array for model '01_multiage_multisite_trap-dep.stan' which has 12 states
# states 4,8,12 are unobservable so the warning of their absence can be ignored
marr <- get_marray(y, nStates=12) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ---- Fit the model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compile and fit model
file <- "stan/01_multiage_multisite_trap-dep.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, marr=marr)
fit <- mod$sample(stan_data, parallel_chains = 4)

# sampler diagnostics
fit$diagnostic_summary()

# plot adult survival estimates
fit$summary("phi_ad") %>%
  mutate(
    site_index = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
    t = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") ),
    year = 2012 + t,
    site = case_when(
      site_index == 1 ~ "Robben",
      site_index == 2 ~ "Boulders",
      site_index == 3 ~ "Stony"
    ),
    xshift = 0.1*case_when(
      site_index == 1 ~ -1,
      site_index == 2 ~ 0,
      site_index == 3 ~ 1
    )
  ) %>%
  ggplot( aes(x=year+xshift, colour=site) ) +
  geom_pointrange( aes(y=median, ymin=q5, ymax=q95), size=0.8, linewidth=0.8 ) +
  scale_x_continuous( breaks=2013:2023 ) +
  coord_cartesian( ylim=c(0,1) ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line()
  ) + 
  labs( x= "year", y="estimate", title = "Adult survival")

# plot juvenile survival estimates
fit$summary("phi_jv") %>%
  mutate(
    site_index = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
    t = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") ),
    year = 2012 + t,
    site = case_when(
      site_index == 1 ~ "Robben",
      site_index == 2 ~ "Boulders",
      site_index == 3 ~ "Stony"
    ),
    xshift = 0.1*case_when(
      site_index == 1 ~ -1,
      site_index == 2 ~ 0,
      site_index == 3 ~ 1
    )
  ) %>%
  ggplot( aes(x=year+xshift, colour=site) ) +
  geom_pointrange( aes(y=median, ymin=q5, ymax=q95), size=0.8, linewidth=0.8 ) +
  scale_x_continuous( breaks=2013:2023 ) +
  coord_cartesian( ylim=c(0,1) ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line()
  ) + 
  labs( x= "year", y="estimate", title = "Juvenile survival")
