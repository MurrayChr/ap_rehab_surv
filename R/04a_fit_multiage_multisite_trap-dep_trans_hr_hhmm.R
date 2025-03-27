#' Fit model with additive hand-rearing effect on survival (the same model as in
#' 03a_...), but using a hierarchical hidden Markov model likelihood.
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

# compile and fit
file <- "stan/04_multiage_multisite_trap-dep_trans_hr_hhmm.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, N=N_rr, y=y_rr, fc=fc_rr, fc_state=fc_state_rr, 
                  mult=mult_rr, hr=hr_rr)
fit <- mod$sample(stan_data, parallel_chains = 4, 
                  iter_warmup = 250, iter_sampling = 250, refresh = 10)
# fit$save_object("outputs/04a_multiage_multisite_trap-dep_trans_hr_hhmm_fit.RDS")

# diagnostic summary
fit$diagnostic_summary()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ---- Plot estimates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# load fitted model object 
# fit <- readRDS("outputs/04a_multiage_multisite_trap-dep_trans_hr_hhmm_fit.RDS")

# posteriors for hand-rearing effects (logit scale)
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
    x = "Hand-rearing effect (logit scale)", 
    y = "Posterior density" 
  )
# ggsave("figs/03a_logit_hr_effects.png", height = 3, width = 5)

# what's 'the right' way to get posterior on the probability scale
fit$draws(c("phi_ad_hr", "phi_ad_wr"), format = "df") %>%
  mutate(
    draw = .draw
  ) %>%
  select(!starts_with(".")) %>%
  pivot_longer(-"draw") %>%
  mutate(
    par = str_extract(name, "[a-zA-Z_]+(?=\\[)"),
    site_index = as.integer( str_extract(name, "(?<=\\[)[1-9](?=,)")  ),
    t = as.integer( str_extract(name,"(?<=,)[0-9]+(?=\\])") )
  ) %>% 
  # maybe I need a reframe in here?
  group_by(site_index, t, draw) %>%
  pivot_wider(names_from = par, values_from = value)
