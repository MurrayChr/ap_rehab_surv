#' Simulation test of model 05_td_tr_hr_phi-im.stan
#' We do not implement full simulation-based calibration (Talts et al. 2020),
#' taking the following approach instead:
#' (i) Choose data-generating 'true' parameter values in one of three ways (see 
#'     the function generate_parameter_values() below for precise details):
#'   - For (logit-scale) hand-rearing 'effects' we sample randomly from Normal(0, 0.5)
#'   - For wild-raised survival probabilities, immature survival probabilities,
#'     movement probabilities and residency probabilities we sample from distributions 
#'     modeled on the posterior of the fitted model 05_td_tr_hr_phi-im.stan
#'   - For detection probabilities we sample randomly from marginal posteriors
#'     of fitted model 05_td_tr_hr_phi-im.stan
#' (ii) Simulate a replicate dataset under each of these.
#' (iii) Fit the model to each of these replicate datasets and compare point
#'       estimates and coverage.
library(tidyverse)
library(cmdstanr)
source("R/00_function_get_marray.R")
source("R/05b_fn_sim_rep_data_trap-dep_trans_hr_phi-im.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Function to draw data-generating parameter values ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_parameter_values <- function(fit, real_cmr_data) {
  # number of years
  T <- sum(str_starts(colnames(real_cmr_data), "yr"))
  
  # survival probabilities for wild-raised birds
  phi_ad_wr <- matrix(rbeta(3*(T-1), 16, 4), nrow = 3, ncol = T-1)
  phi_jv_wr <- matrix(rbeta(3*(T-1), 12, 15), nrow = 3, ncol = T-1)
  
  # hand-rearing 'effects'
  hr_jv <- rnorm(1, 0, 0.5)
  hr_ad <- rnorm(1, 0, 0.5)
  
  # immature survival
  phi_im <- rbeta(T-1, 15, 5)
  
  # movement probabilities
  m_jv <- matrix(NA, nrow = 3, ncol = 3)
  m_ad <- matrix(NA, nrow = 3, ncol = 3)
  alpha_jv <- matrix(c(
    64, 8, 8,
    8, 64, 8,
    8, 8, 64
  ), nrow = 3, ncol = 3, byrow = TRUE)
  alpha_ad <- matrix(c(
    95, 2.5, 2.5,
    2.5, 95, 2.5,
    2.5, 2.5, 95
  ), nrow = 3, ncol = 3, byrow = TRUE)
  for (i in 1:3) {
    # draw from Dirichlet by drawing indep gammas and normalising
    m_jv[i,] <- rgamma(3, shape = alpha_jv[i,])
    m_jv[i,] <- m_jv[i,]/sum(m_jv[i,])
    m_ad[i,] <- rgamma(3, shape = alpha_ad[i,])
    m_ad[i,] <- m_ad[i,]/sum(m_ad[i,])
  }
  
  # residency probabilities drawn from distributions roughly modeled on
  # posterior distributions from the fit to real data: for most site-years
  # we use a Beta(9,1), but for a few site-years where model 01's goodness-of-fit 
  # flagged transients and the residence probabilities were lower, we used 
  # different distributions
  pi_r <- matrix(rbeta(3*(T-1), 9, 1), nrow = 3, ncol = T-1)
  pi_r[2,5] <- rbeta(1, 37, 63) 
  pi_r[2,7] <- rbeta(1, 60, 40)  
  
  # detection probabilities drawn randomly from the posterior samples
  post <- fit$draws(
    variables = c("p_im", "p_ad_A", "p_ad_U", "p_ad_N"), 
    format = "df"
  )
  p_im <- matrix(NA, nrow = 3, ncol = T)
  p_ad_A <- matrix(NA, nrow = 3, ncol = T)
  p_ad_U <- matrix(NA, nrow = 3, ncol = T)
  p_ad_N <- matrix(NA, nrow = 3, ncol = T)
  for (i in 1:3) {
    for (t in 1:T) {
      p_im_var <- str_c("p_im[",i,",",t,"]")
      p_im[i,t] <- sample(post[[p_im_var]], 1)
      p_ad_A_var <- str_c("p_ad_A[",i,",",t,"]")
      p_ad_A[i,t] <- sample(post[[p_ad_A_var]], 1)
      p_ad_U_var <- str_c("p_ad_U[",i,",",t,"]")
      p_ad_U[i,t] <- sample(post[[p_ad_U_var]], 1)
      p_ad_N_var <- str_c("p_ad_N[",i,",",t,"]")
      p_ad_N[i,t] <- sample(post[[p_ad_N_var]], 1)
    }
  }
  
  # return list of parameter values
  list(
    phi_ad_wr = phi_ad_wr,
    phi_jv_wr = phi_jv_wr,
    hr_jv = hr_jv,
    hr_ad = hr_ad,
    phi_im = phi_im,
    m_jv = m_jv,
    m_ad = m_ad,
    pi_r = pi_r,
    p_im = p_im,
    p_ad_A = p_ad_A,
    p_ad_U = p_ad_U,
    p_ad_N = p_ad_N
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        ---- Load fitted model and simulate replicate datasets ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data and model fit to data
real_cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")
fit <- readRDS("outputs/05a_td_tr_hr_phi-im_fit.RDS")

# remove individuals first captured on last occasion
T <- sum(str_starts(names(real_cmr_data), "yr"))
real_cmr_data <- filter( real_cmr_data, fc!=T )

# define dataset structure; replicate data will have the same structure
data_str <- list(
  T = sum(str_starts(names(real_cmr_data), "yr")),
  N = nrow(real_cmr_data),
  fc = real_cmr_data$fc,
  fc_code = real_cmr_data$fc_code,
  hr = real_cmr_data$hr
)

# simulate replicate data and save to outputs/sim_05/rep_data
n_reps <- 100
for (i in 1:n_reps) {
  truth <- generate_parameter_values(fit, real_cmr_data)
  rep_data <- sim_cmr_data(pars = truth, data_str = data_str)
  file <- str_c("outputs/sim_05/rep_data/05c_", str_pad(i, 3, "left", "0"), ".RDS")
  saveRDS(rep_data, file)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             ---- Function to prepare data for Stan model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_stan_data <- function(rep_cmr_data) {
  # format rep data like real data
  y <- rep_cmr_data$y
  colnames(y) <- str_c("yr", 2012 + 1:ncol(y))
  cmr_data <- as_tibble(y) %>%
    add_column(
      fc = rep_cmr_data$fc,
      fc_code = rep_cmr_data$fc_code,
      hr = rep_cmr_data$hr
    ) %>%
    mutate(
      marking_age = case_when(
        fc_code %in% c(1, 5, 9) ~ 1,    # juveniles
        fc_code %in% c(3, 7, 11) ~ 3    # adults
      )
    )
  
  # set number of years
  T <- sum(str_starts(colnames(cmr_data),"yr"))
  
  # add last capture occasion to find individuals not recaptured after first capture
  cmr_data <- cmr_data %>%
    rowwise() %>%
    mutate( 
      lc = max( which( c_across( starts_with("yr") ) != 13 ) ), # THIS WAS '0' AT ONE STAGE AND BUGGED OUT!!!
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
  
  # m-array for model '05_td_tr_hr_phi-im.stan' which has 12 states
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
  
  return(list(T = T, marr_wr = marr_wr, marr_hr = marr_hr, N_1 = N_1, N_0 = N_0))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     ---- Fit the models ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fit each replicate dataset and save to outputs/sim_05/fitted_models
mod_file <- "stan/05_td_tr_hr_phi-im.stan"
mod <- cmdstan_model(mod_file)
for (i in 1:n_reps) {
  file_num <-  str_pad(i, 3, "left", "0")
  cmr_data_file <- str_c("outputs/sim_05/rep_data/05c_", file_num, ".RDS")
  rep_cmr_data <- readRDS(cmr_data_file)
  stan_data <- get_stan_data(rep_cmr_data)
  fit <- mod$sample(stan_data, parallel_chains = 4, 
                    iter_warmup = 500, iter_sampling = 500)
  fit_file <- str_c("outputs/sim_05/fitted_models/05c_", file_num, "_fit.RDS")
  fit$save_object(fit_file)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Check model diagnostiscs across simulations ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n_reps <- 100
n_div <- rep(NA, n_reps)
for (i in 1:n_reps) {
  # read in sim and fit
  file_num <-  str_pad(i, 3, "left", "0")
  fit <- readRDS(str_c("outputs/sim_05/fitted_models/05c_", file_num, "_fit.RDS"))
  n_div[i] <- fit$diagnostic_summary()$num_divergent %>% sum()
}
table(n_div) # majority have some divergences - we might need to rerun the most egregious of these

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Compare estimates to truth in a single replicate ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_par_vs_truth <- function(rep_cmr_data, fit, par) {
  par_truth <- rep_cmr_data$truth[[par]] 
  years <- 2012 + 1:ncol(par_truth)
  colnames(par_truth) <- years
  par_truth <- as_tibble(par_truth) %>%
    add_column(site = 1:3) %>%
    pivot_longer(-"site") %>%
    mutate(
      year = as.integer(name)
    )
  
  fit$summary(par) %>%
    mutate(
      site = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
      year = 2012L + as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") )
    ) %>%
    ggplot(aes(x = year)) +
    geom_pointrange(aes(y= median, ymin = q5, ymax = q95)) +
    geom_point(data = par_truth, mapping = aes(y = value), colour = "red") +
    scale_x_continuous(breaks = years) +
    facet_grid(site ~ .) +
    labs(title = par, y = "")
}

# select replicate and load data and fit
i <- 1
file_num <-  str_pad(i, 3, "left", "0")
rep_cmr_data <- readRDS(str_c("outputs/sim_05/rep_data/05c_", file_num, ".RDS"))
fit <- readRDS(str_c("outputs/sim_05/fitted_models/05c_", file_num, "_fit.RDS"))

# plots
plot_par_vs_truth(rep_cmr_data, fit, "phi_jv_wr")
plot_par_vs_truth(rep_cmr_data, fit, "phi_ad_wr")
plot_par_vs_truth(rep_cmr_data, fit, "pi_r")  
plot_par_vs_truth(rep_cmr_data, fit, "p_im")
plot_par_vs_truth(rep_cmr_data, fit, "p_ad_A")
plot_par_vs_truth(rep_cmr_data, fit, "p_ad_U")
plot_par_vs_truth(rep_cmr_data, fit, "p_ad_N")

fit$summary(c("hr_jv", "hr_ad"))
rep_cmr_data$truth[c("hr_jv", "hr_ad")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            ---- Estimates vs truth across replicates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# uncomment to select a paramter
# par <- "phi_ad_wr"
par <- "phi_jv_wr"

# collect truth and estimates from each simulation
n_reps <- 100
est <- tibble()
for (i in 1:n_reps) {
  # read in sim and fit
  file_num <-  str_pad(i, 3, "left", "0")
  rep_cmr_data <- readRDS(str_c("outputs/sim_05/rep_data/05c_", file_num, ".RDS"))
  fit <- readRDS(str_c("outputs/sim_05/fitted_models/05c_", file_num, "_fit.RDS"))
  
  # get posterior summaries for parameters
  sim_est <- fit$summary(par)[,c("variable", "median", "q5", "q95")] %>%
    mutate(
      site = as.integer(str_extract(variable, "(?<=\\[)[1-9](?=,)")),
      .after = variable
    ) %>%
    mutate(
      t = as.integer(str_extract(variable,"(?<=,)[0-9]+(?=\\])")),
      .after = site
    )
  
  # get true data-generating values and format as sim_est
  sim_truth <- rep_cmr_data$truth[[par]] 
  colnames(sim_truth) <- 1:ncol(sim_truth)   
  sim_truth <- as_tibble(sim_truth) %>%
    add_column(site = 1:3) %>%
    pivot_longer(-"site", names_to = "t", values_to = "truth") %>%
    mutate(t = as.integer(t))
  
  # join estimates and truth
  sim_tib <- left_join(sim_est, sim_truth, by = c("site", "t")) %>%
    add_column(sim_id = i, .before = "variable")
  sim_tib
  
  # add values from this sim to others
  est <- rbind(est, sim_tib)
}

# add year and site name (used to label facets in plot)
est <- est %>%
  mutate(
    year = t + 2012,
    site_name = case_when(
      site == 1 ~ "Robben",
      site == 2 ~ "Boulders",
      site == 3 ~ "Stony"
    )
  )

# plot
plt <- est %>%
  ggplot(aes(x = truth)) +
  geom_point(aes(y = median), size = 0.5) +
  geom_linerange(aes(ymin = q5, ymax = q95), alpha = 0.2 ) +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.1) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1,length.out = 3)) +
  scale_y_continuous(breaks = seq(0,1,length.out = 3)) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3)
  ) +
  facet_grid(factor(site_name, levels = c("Robben", "Boulders", "Stony")) ~ year) + # convert site_name to factor so that levels can be ordered manually
  labs(
    x = "Truth",
    y = "Estimate",
    title = case_when(
      par == "phi_ad_wr" ~ "Adult survival, wild-raised birds",
      par == "phi_jv_wr" ~ "Juvenile survival, wild-raised birds"
    )
  )
plt

# calculate coverage, bias and root mean squared error
est_bc <- est %>%
  mutate(
    cover = (q5 < truth) & (truth < q95),
    bias = median - truth,
    bias_squared = bias^2
  ) %>%
  group_by(year, site_name) %>%
  summarise(
    coverage = mean(cover),
    bias = mean(bias),
    rmse = sqrt(mean(bias_squared)),
  ) %>%
  ungroup() 

bc_text <- est_bc %>%
  mutate(
    label = str_c(
      "cover = ", coverage, "\n",
      "bias = ", round(bias,2), "\n",
      "rmse = ", signif(rmse, 1)
    )
  ) %>%
  add_column(x = 0.75, y = 0.13)

plt +
  geom_text(
    data = bc_text,
    mapping = aes(x = x, y = y, label = label),
    colour = "grey30",
    size = 2
  )
# ggsave(str_c("figs/05c_",par,"_sim_results.png"), height = 6, width = 18)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            ---- Hand-rearing parameters across replicates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hr_est <- tibble()
n_reps <- 100
for (i in 1:n_reps) {
  # read in sim and fit
  file_num <-  str_pad(i, 3, "left", "0")
  rep_cmr_data <- readRDS(str_c("outputs/sim_05/rep_data/05c_", file_num, ".RDS"))
  fit <- readRDS(str_c("outputs/sim_05/fitted_models/05c_", file_num, "_fit.RDS"))
  # get posterior summaries for hr_ parameters
  sim_est <- fit$summary(c("hr_jv", "hr_ad"))[,c("variable", "median", "q5", "q95")]
  # get true data-generating values and format as sim_est
  sim_truth <- rep_cmr_data$truth[c("hr_ad", "hr_jv")] %>% 
    as_tibble() %>% 
    pivot_longer(everything(), names_to = "variable", values_to = "truth")
  # join estimates and truth
  sim_tib <- left_join(sim_est, sim_truth, by = "variable") %>%
    add_column(sim_id = i, .before = "variable")
  # add values from this sim to others
  hr_est <- rbind(hr_est, sim_tib)
}

# plot
facet_labs <- c("hr_ad" = "Adults", "hr_jv" = "Juveniles")
hr_plot <- hr_est %>%
  ggplot(aes(x = truth)) +
  geom_point(aes(y = median), size = 0.5) +
  geom_linerange(aes(ymin = q5, ymax = q95), alpha = 0.2 ) +
  scale_x_continuous(breaks = seq(-1,1, length.out = 5)) +
  scale_y_continuous(breaks = seq(-1,1, length.out = 5)) +
  coord_fixed(xlim=c(-1.5,1.5), ylim = c(-1.5,1.5)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(linewidth = 0.15)
  ) +
  facet_wrap(vars(variable), labeller = as_labeller(facet_labs)) +
  labs(x = "Truth", y = "Estimate", 
       title = "Age-dependent hand-rearing differences (logit-scale)" )
hr_plot

# compute bias, coverage and rmse
hr_bc <- hr_est %>%
  mutate(
    bias = median - truth,
    cover = (truth > q5) & (truth < q95),
  ) %>%
  group_by(variable) %>%
  summarise(
    rmse = sqrt(mean((bias)^2)),
    bias = mean(bias),
    coverage = mean(cover),
    n_covered = sum(cover)
  ) %>%
  ungroup()

# code to add annotations from the R Language-recommended answer at
# https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2

bc_text <- hr_bc %>%
  mutate(
    label = str_c(
      "cover = ", coverage, "\n",
      "bias = ", signif(bias,1), "\n",
      "rmse = ", signif(rmse, 1)
    )
  ) %>%
  add_column(x = 0.75, y = -1)

hr_plot +
  geom_text(
    data = bc_text,
    mapping = aes(x = x, y = y, label = label),
    colour = "grey20",
    size = 2.5
  ) 
# ggsave("figs/05c_hr_sim_results.png", scale = 1.5)
