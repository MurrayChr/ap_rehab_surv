#' Fit model with additive hand-rearing effect on survival (the same model as in
#' 04a_fit_trap-dep_trans_hr.R), but using a hierarchical hidden Markov model likelihood.
library(tidyverse)
library(cmdstanr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           ---- Load cmr data and create reduced representation ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# set number of years
T <- sum(str_starts(colnames(cmr_data),"yr"))

# recode hr indicator, which is a factor, as an index
cmr_data <- cmr_data %>%
  mutate(hr = ifelse(hr==1, 1, 2)) # so hr = 1, wr = 2

# remove individuals first captured in the last year
cmr_data <- filter(cmr_data, fc < T)

# get reduced representation of the data (ADD Notes about same capt hist distinct if hr differs)
reduced_representation <- cmr_data %>%
  ungroup() %>%
  select(starts_with("yr"), fc, hr) %>%
  group_by_all() %>%
  summarise(n=n()) %>%
  ungroup()
y_rr <- reduced_representation %>%
  select(starts_with("yr")) %>%
  as.matrix()
fc_rr <- reduced_representation$fc
mult_rr <- reduced_representation$n
hr_rr <- reduced_representation$hr

# state at first capture
N_rr <- nrow(y_rr)
fc_state_rr <- rep(NA, N_rr)
for (n in 1:N_rr) {
  fc_state_rr[n] <- y_rr[n, fc_rr[n]]
}

# change 0s after first capture to not captured code '13'
for (n in 1:N_rr) {
  y_rr[n, (fc_rr[n]+1):T][y_rr[n, (fc_rr[n]+1):T]==0] <- 13
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ---- Fit the model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compile and fit (~ 2 hours)
file <- "stan/04_td_tr_hr_hhmm.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, N=N_rr, y=y_rr, fc=fc_rr, fc_state=fc_state_rr, 
                  mult=mult_rr, hr=hr_rr)
fit <- mod$sample(stan_data, parallel_chains = 4)
# fit$save_object("outputs/04a_td_tr_hr_hhmm_fit.RDS")

# diagnostic summary
fit$diagnostic_summary()

# rhats and ess's
fit_summary <- fit$summary()
max(fit_summary$rhat)
min(fit_summary$ess_bulk)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   ---- Compare posteriors ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_hhmm <- readRDS("outputs/04a_td_tr_hr_hhmm_fit.RDS")
fit_epm <- readRDS("outputs/04a_td_tr_hr_fit.RDS")

# get all parameters common to two models (all parameters in this case)
get_common_params <- function(fit_1, fit_2) {
  params_1 <- fit_1$metadata()$model_params
  params_2 <- fit_2$metadata()$model_params
  intersect(params_1, params_2)
}
common_params <- get_common_params(fit_epm, fit_hhmm)
common_params <- common_params[common_params != "lp__"]

# format posterior draws and add model name
prepare_draws <- function(fit, params, name) {
  fit$draws(params, format = "df") %>%
    add_column(lk = name) %>%
    select(-starts_with(".")) %>%
    pivot_longer(-"lk")
}
post_epm <- prepare_draws(fit_epm, common_params, "epm")
post_hhmm <- prepare_draws(fit_hhmm, common_params, "hhmm")
post <- rbind(post_epm, post_hhmm) 

# plots comparing marginal posterior densities
post %>%
  filter(str_starts(name, "phi_ad")) %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density(bounds=c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free") +
  labs(title = "Adult survival probabilities")

post %>%
  filter(str_starts(name, "phi_jv")) %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density(bounds=c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free") +
  labs(title = "Juvenile survival probabilities")

post %>%
  filter(str_starts(name, "pi_r")) %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density(bounds=c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free") +
  labs(title = "Residency probabilities")

post %>%
  filter(str_starts(name, "m_")) %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density(bounds=c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free") +
  labs(title = "Movement probabilities")

post %>%
  filter(str_starts(name, "hr")) %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density() +
  # coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free") +
  labs(title = "Hand-rearing effects")

post %>%
  filter(str_starts(name, "p_")) %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density(bounds=c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free", ncol = 11) +
  labs(title = "Detection probabilities")


