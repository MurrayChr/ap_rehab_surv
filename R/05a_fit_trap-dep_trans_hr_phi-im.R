# Fit model with trap-dependence, transience and additive hand-rearing effect 
# on survival which additionally includes separate survival for immatures 
# (as opposed to them sharing survival with adults)
library(tidyverse)
library(cmdstanr)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           ---- Load cmr data and format into m-arrays ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# set number of years
T <- sum(str_starts(colnames(cmr_data),"yr"))

# add last capture occasion to find individuals not recaptured after first capture
cmr_data <- cmr_data %>%
  rowwise() %>%
  mutate( 
    lc = max( which( c_across( starts_with("yr") )!=0 ) ),
    .after = fc
  ) %>%
  ungroup()

# split data into three parts
cmr_data_not_for_marr <- cmr_data %>%
  filter(
    (marking_age == 3) & (fc == lc)
  )
cmr_data_for_marr_wr <- cmr_data %>%
  filter(
    !((marking_age == 3) & (fc == lc)) & (hr == 0)
  )
cmr_data_for_marr_hr <- cmr_data %>%
  filter(
    !((marking_age == 3) & (fc == lc)) & (hr == 1)
  )

# check that these add up
# nrow(cmr_data_for_marr_hr) + nrow(cmr_data_for_marr_wr) +
#   nrow(cmr_data_not_for_marr) == nrow(cmr_data)

# check that there are no marking_age == 3 in cmr_data_for_marr_hr; all
# hand-reared birds should be marked as juveniles, not adults
# all(cmr_data_for_marr_hr$marking_age != 3)

# prepare data for get_marray()
prepare_cmr_data <- function(cmr_data_for_marr) {
  # split into capture history matrix and first-capture vector
  y <- cmr_data_for_marr %>%
    select( starts_with("yr") ) %>%
    as.matrix()
  fc <- cmr_data_for_marr$fc
  
  # change '0's after first capture to nStates + 1 = 13 for get_marray()
  for (i in 1:nrow(y)) {
    y[i,fc[i]:T][y[i, fc[i]:T]==0] <- 13 
  }
  
  # return
  y
}
y_wr <- prepare_cmr_data(cmr_data_for_marr_wr)
y_hr <- prepare_cmr_data(cmr_data_for_marr_hr)

# m-array for model '02_multiage_multisite_trap-dep_trans.stan' which has 12 states
# states 4,8,12 are unobservable so the warning about their absence can be safely ignored
marr_wr <- get_marray(y_wr, nStates=12)
marr_hr <- get_marray(y_hr, nStates=12)

# count number of birds marked as adults that are ever/never recaptured
# per site and cohort
N_1 <- matrix(NA, 3, T-1)
N_0 <- matrix(NA, 3, T-1)
for (i in 1:3) {
  for (t in 1:(T-1)) {
    N_1[i,t] <- cmr_data_for_marr_wr %>%
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

# compile and fit (~ 17 minutes)
file <- "stan/05_td_tr_hr_phi-im.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, marr_wr=marr_wr, marr_hr = marr_hr, N_1=N_1, N_0=N_0)
fit <- mod$sample(stan_data, parallel_chains = 4)
fit$save_object("outputs/05a_td_tr_hr_phi-im_fit.RDS")

# diagnostic summary
fit$diagnostic_summary()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Plot estimates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# load fitted model object 
fit <- readRDS("outputs/05a_td_tr_hr_phi-im_fit.RDS")

# posteriors for hand-rearing coefficient (logit scale)
fit$draws(c("hr_jv", "hr_ad"), format = "df") %>%
  select(-starts_with(".")) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(-1,1)) +
  geom_vline(xintercept = 0, colour="grey", linetype = "dotted") +
  scale_fill_manual(values=c("darkorange", "navyblue"), labels = c("adult", "juvenile"), name = "") +
  theme_classic() +
  theme(legend.position = "inside", legend.position.inside = c(0.7,0.85)) +
  labs(
    # title = "Age-specific hand-rearing effects",
    x = "Hand-rearing coefficient (logit scale)", 
    y = "Posterior density" 
  )


# immature survival
fit$summary("phi_im") %>%
  mutate(
    par = str_extract(variable, "[a-z_]+(?=\\[)"),
    t = as.integer(str_extract(variable,"(?<=\\[)[0-9]+(?=\\])")),
  ) %>%
  ggplot(aes(x = t)) +
  geom_pointrange(aes(y=median, ymin = q5, ymax = q95)) +
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(colour = "grey95")
  ) +
  labs(
    x = "Year", y = "Posterior median and 90% CrI", 
    title = "Immature survival",
    subtitle = "Shared between sites and hand-rearing categories, with year random effect"
  )





