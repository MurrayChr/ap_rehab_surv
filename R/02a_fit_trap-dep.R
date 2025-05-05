# Fit a multi-age, multi-site model that incorporates trap-dependence but not
# transience (and without a hand-reading covariate).
library(tidyverse)
library(cmdstanr)
library(cowplot)
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

# m-array for model '02_td.stan' which has 12 states
# states 4,8,12 are unobservable so the warning about their absence can be safely ignored
marr <- get_marray(y, nStates=12) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ---- Fit the model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compile and fit model (~ 8 minutes)
file <- "stan/02_td.stan"
mod <- cmdstan_model(file)
stan_data <- list(T=T, marr=marr)
fit <- mod$sample(stan_data, parallel_chains = 4) 
# fit$save_object("outputs/02a_td_fit.RDS")

# sampler diagnostics
fit$diagnostic_summary()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        ---- Plot estimates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit <- readRDS("outputs/02a_td_fit.RDS")

# plot adult survival estimates
plt_phi_ad <-fit$summary("phi_ad") %>%
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
    legend.position.inside = TRUE,
    legend.position = c(0.75,0.25)
  ) + 
  labs( x= "year", y="estimate", title = "Adult survival")
plt_phi_ad

# plot juvenile survival estimates
plt_phi_jv <- fit$summary("phi_jv") %>%
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
    legend.position.inside = TRUE,
    legend.position = c(0.225,0.8)
  ) + 
  labs( x= "year", y="estimate", title = "Juvenile survival")
plt_phi_jv

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
    legend.position.inside = TRUE,
    legend.position = c(0.775,0.25)
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
    legend.position.inside = TRUE,
    legend.position = c(0.375,0.7)
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

# arrange plots in grid and save
plt_estimates <- plot_grid(
  plt_phi_ad, plt_phi_jv, plt_p_ad_A, plt_p_im, plt_p_ad_U, plt_p_ad_N,
  plt_mv,
  ncol = 2
)
# ggsave(plot = plt_estimates, "figs/02a_estimates.pdf", height = 13.3, width = 10)

