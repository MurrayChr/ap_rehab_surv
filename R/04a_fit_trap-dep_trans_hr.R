# Fit model with trap-dependence, transience and additive hand-rearing effect 
# on survival
library(tidyverse)
library(cowplot)
library(cmdstanr)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           ---- Load cmr data and format into m-arrays ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# set number of years
T <- sum(str_starts(colnames(cmr_data),"yr"))

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

# m-array for model '04_td_tr_hr.stan' which has 12 states for each hand-rearing category
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

# compile and fit (~ 18 minutes)
file <- "stan/04_td_tr_hr.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, marr_wr=marr_wr, marr_hr = marr_hr, N_1=N_1, N_0=N_0)
fit <- mod$sample(stan_data, parallel_chains = 4)
# fit$save_object("outputs/04a_td_tr_hr_fit.RDS")

# sampler diagnostics
fit$diagnostic_summary()

# rhats and ess's
fit_summary <- fit$summary() %>% filter(variable != "lp__") 
max(fit_summary$rhat)
min(fit_summary$ess_bulk)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ---- Plot estimates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# load fitted model object 
fit <- readRDS("outputs/04a_td_tr_hr_fit.RDS")

# plot adult survival estimates for wild-raised birds
plt_phi_ad_wr <- fit$summary("phi_ad_wr") %>%
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
    panel.grid.major = element_line(),
    legend.position = "inside",
    legend.position.inside = c(0.75,0.25)
  ) + 
  labs( x= "year", y="estimate", title = "Adult survival for wild-raised birds")
plt_phi_ad_wr

# plot adult survival estimates for hand-reared birds
plt_phi_ad_hr <- fit$summary("phi_ad_hr") %>%
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
    panel.grid.major = element_line(),
    legend.position = "none"
  ) + 
  labs( x= "year", y="estimate", title = "Adult survival for hand-reared birds")
plt_phi_ad_hr

# plot juvenile survival estimates for wild-raised birds
plt_phi_jv_wr <- fit$summary("phi_jv_wr") %>%
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
    panel.grid.major = element_line(),
    legend.position = "inside",
    legend.position.inside = c(0.225,0.8)
  ) + 
  labs( x= "year", y="estimate", title = "Juvenile survival for wild-raised birds")
plt_phi_jv_wr

# plot juvenile survival estimates for hand-reared birds
plt_phi_jv_hr <- fit$summary("phi_jv_hr") %>%
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
    panel.grid.major = element_line(),
    legend.position = "none"
  ) + 
  labs( x= "year", y="estimate", title = "Juvenile survival for hand-reared birds")
plt_phi_jv_hr

# plot 'residency' probabilities for newly marked adults
plt_pi_r <- fit$summary("pi_r") %>%
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
    panel.grid.major = element_line(),
    legend.position = "inside",
    legend.position.inside = c(0.775,0.225)
  ) + 
  labs( x= "year", y="estimate", title = "Residency probability")
plt_pi_r

# plot the detection parameters for 'trap-aware' adults
plt_p_ad_A <- fit$summary("p_ad_A") %>%
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
  scale_x_continuous( breaks=2013:2024 ) +
  coord_cartesian( ylim=c(0,1) ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(),
    legend.position = "inside",
    legend.position.inside = c(0.775,0.25)
  ) + 
  labs( x= "year", y="estimate", title = "Detection of trap-aware adults")
plt_p_ad_A

# plot the detection parameters for 'trap-unaware' adults
plt_p_ad_U <- fit$summary("p_ad_U") %>%
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
  scale_x_continuous( breaks=2013:2024 ) +
  coord_cartesian( ylim=c(0,1) ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(),
    legend.position = "none"
  ) + 
  labs( x= "year", y="estimate", title = "Detection of trap-unaware adults")
plt_p_ad_U

# plot the detection parameters for 'newly matured' adults 
plt_p_ad_N <- fit$summary("p_ad_N") %>%
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
  scale_x_continuous( breaks=2013:2024 ) +
  coord_cartesian( ylim=c(0,1) ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(),
    legend.position = "none"
  ) + 
  labs( x= "year", y="estimate", title = "Detection of newly-matured adults")
plt_p_ad_N

# immature detection
plt_p_im <- fit$summary("p_im") %>%
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
  scale_x_continuous( breaks=2013:2024 ) +
  coord_cartesian( ylim=c(0,1) ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(),
    legend.position = "inside",
    legend.position.inside = c(0.375,0.7)
  ) + 
  labs( x= "year", y="estimate", title = "Detection of immatures")
plt_p_im

# movement parameters
plt_mv <- fit$draws(format = "df") %>%
  select( starts_with("m_") ) %>%
  pivot_longer( everything(), names_to = "variable", values_to = "draw" ) %>%
  mutate(
    site_from_index = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
    site_to_index = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") ),
    site_from = case_when(
      site_from_index == 1 ~ "Robben",
      site_from_index == 2 ~ "Boulders",
      site_from_index == 3 ~ "Stony"
    ),
    site_to = case_when(
      site_to_index == 1 ~ "Robben",
      site_to_index == 2 ~ "Boulders",
      site_to_index == 3 ~ "Stony"
    ),
    age_class = str_extract(variable, "m_[a-z]+")
  ) %>% 
  ggplot( aes(x = draw, fill = age_class ) ) +
  geom_density( bounds = c(0,1), alpha = 0.5 ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank() 
  ) +
  scale_fill_manual(
    values=c("m_ad"="red", "m_jv"="blue"),
    labels=c("m_ad"="adult", "m_jv"="jv. & imm.")
  ) +
  coord_cartesian( xlim = c(0,1), ylim = c(0,30) ) +
  facet_grid(rows = vars(site_from), cols = vars(site_to), switch = "y" ) +
  labs(
    title = "Posterior distributions for movement probabilities",
    y = "Density",
    x = "Probability",
    fill = "Age class"
  )
plt_mv

# posteriors for hand-rearing coefficient (logit scale)
plt_hr_diffs <- fit$draws(c("hr_jv", "hr_ad"), format = "df") %>%
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
plt_hr_diffs

# arrange plots in grid and save
plt_estimates <- plot_grid(
  plt_phi_ad_wr, plt_phi_ad_hr,
  plt_phi_jv_wr, plt_phi_jv_hr,
  plt_p_ad_A, plt_p_im, 
  plt_p_ad_U, plt_p_ad_N,
  plt_mv, plt_pi_r, 
  plt_hr_diffs,
  ncol = 2
)
# ggsave(plot = plt_estimates, "figs/04a_estimates.pdf", height = 20, width = 10)
