#' Single-site, multi-age, model comparing survival between hand-reared and 
#' wild-raised birds. The model is fit to a subset of the data that is more 
#' balanced than the entire dataset, namely, all birds 'belonging' to Stony that
#' were marked as juveniles. 
#' 
#' The model
#'   (i) incorporates trap-dependence in the adult states
#'   (ii) does not incorporate transience in adult states, because no birds marked
#'        as adults are included in this dataset,
#'   (iii) shares immature survival with adult survival; although this could bias 
#'         adult and juvenile survival as discussed in the manuscript, it should
#'         not bias the hand-rearing vs wild-raised comparison because each group
#'         consists purely of birds marked as juveniles in this dataset.
#'   (iv) includes age-specific movement probabilities via informative priors set 
#'        using the fitted model 05a_td_tr_hr_phi-im; these could not be estimated
#'        from this data along (because of the usual confounding of survival and
#'        permanent emigration), but can be used to reduce bias in survival from
#'        emigration from Stony to Robben or Boulders.
#'   
library(tidyverse)
library(cmdstanr)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Load, subset and prepare data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data with states coded with trap-dependence
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# filter for Stony birds marked as juveniles
cmr_data <- cmr_data %>%
  filter(
    marking_site == "Stony",
    marking_age == 1
  ) 

# set number of years
T <- sum(str_starts(colnames(cmr_data),"yr"))

# extract capture history matrix
y <- cmr_data %>%
  select( starts_with("yr") ) %>%
  as.matrix()
fc <- cmr_data$fc
hr <- cmr_data$hr

# recode entries for get_marray()
# the codes in y are:
table(c(y))
# first, all codes relating to Robben and Boulders changed to zero
y[y %in% 1:8] <- 0
# second, codes for Stony mapped to 1:4
y[y==9] <- 1    # juvenile
y[y==10] <- 2   # immature
y[y==11] <- 3   # adult, trap-aware
y[y==12] <- 4   # adult, trap-unaware

# change zeros after first capture to nStates + 1 = 5
for (i in 1:nrow(y)) {
  y[i,fc[i]:T][y[i, fc[i]:T]==0] <- 5
}

# separate hand-reared and wild-raised birds
y_wr <- y[hr==0,]
y_hr <- y[hr==1,]

# create m-arrays for model 07_td_hr.stan
# the trap-unaware state '4' is unobservable, so the warnings can be safely ignored
marr_wr <- get_marray(y_wr, nStates=4)
marr_hr <- get_marray(y_hr, nStates=4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              ---- Create informative prior for movement ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_ <- readRDS("outputs/05a_td_tr_hr_phi-im_fit.RDS")
post_ <- fit_$draws(c("m_jv", "m_ad"), format = "df")
post_ <- post_ %>%
  select("m_jv[3,3]", "m_ad[3,3]") %>% # probabilities for Stony -> Stony movement
  rename("m_jv" = "m_jv[3,3]", "m_ad" = "m_ad[3,3]")

get_beta_dist_pars <- function(mu, sig2) {
  alpha_plus_beta <- ((mu * (1-mu))/sig2) - 1
  alpha <- mu * alpha_plus_beta
  beta <- (1 - mu) * alpha_plus_beta 
  list("alpha" = alpha, "beta" = beta)
}
m_jv_pars <- get_beta_dist_pars(mean(post_$m_jv), var(post_$m_jv))
m_ad_pars <- get_beta_dist_pars(mean(post_$m_ad), var(post_$m_ad))

post_models <- tibble(
  m_ad = rbeta(4000, m_ad_pars$alpha, m_ad_pars$beta),
  m_jv = rbeta(4000, m_jv_pars$alpha, m_jv_pars$beta)
)

# compare posteriors to their moment-matched beta distributions
rbind(
  post_ %>% add_column(draws = "post"),
  post_models %>% add_column(draws = "matched")
) %>%
  pivot_longer(starts_with("m_")) %>%
  ggplot(aes(x=value, colour = draws)) +
  geom_density() +
  # coord_cartesian(xlim=c(0,1)) +
  theme_classic() +
  theme(
    legend.position.inside = TRUE,
    legend.position = c(0.8, 0.75)
  ) +
  facet_wrap(vars(name), scales = "fixed")
# looks good enough

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        ---- Fit model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compile and fit (~ 2 minutes)
file <- "stan/07_stony_td_hr.stan"
mod <- cmdstan_model(file)
stan_data <- list(
  T=T, 
  marr_wr=marr_wr, marr_hr = marr_hr, 
  m_jv_alpha = m_jv_pars$alpha, m_jv_beta = m_jv_pars$beta,
  m_ad_alpha = m_ad_pars$alpha, m_ad_beta = m_ad_pars$beta
)
fit <- mod$sample(stan_data, parallel_chains = 4)
# fit$save_object("outputs/07a_stony_td_hr_fit.RDS")

# diagnostic summary
fit$diagnostic_summary()

# rhats and ess's
fit_summary <- fit$summary()
max(fit_summary$rhat)
min(fit_summary$ess_bulk)
min(fit_summary$ess_tail)
# great

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         ---- Plot estimates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit <- readRDS("outputs/07a_stony_td_hr_fit.RDS")

# look at the hand-rearing estimates
fit$draws(c("hr_jv", "hr_ad"), format = "df") %>%
  select(-starts_with(".")) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(-1,1)) +
  geom_vline(xintercept = 0, colour="grey", linetype = "dotted") +
  scale_fill_manual(
    values=c("darkorange", "navyblue"),
    labels = c("adult", "juvenile"), 
    name = ""
  ) +
  theme_classic() +
  theme(legend.position = "inside", legend.position.inside = c(0.7,0.85)) +
  labs(
    # title = "Age-specific hand-rearing effects",
    x = "Hand-rearing coefficient (logit scale)", 
    y = "Posterior density" 
  )
