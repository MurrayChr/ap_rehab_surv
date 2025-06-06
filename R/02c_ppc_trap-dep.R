#' Posterior predictive checks of model fit for '02_td.stan' 
#' to the data 'data/00b_cmr_data_multisite_multiage_trapdep.RDS'
library(tidyverse)
library(cmdstanr)
source("R/02b_fn_sim_rep_data_trap-dep.R")
source("R/00_functions_get_gof_tables.R")
source("R/00_function_get_expected_frequencies.R")
source("R/00_function_get_ft.R")
source("R/00_function_fit_mixtures.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Load real data and simulate replicate datasets ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data and model fit to data
real_cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")
fit <- readRDS("outputs/02a_td_fit.RDS")

# remove individuals first captured on last occasion
T <- sum(str_starts(names(real_cmr_data), "yr"))
real_cmr_data <- filter( real_cmr_data, fc!=T )

# define dataset structure; replicate data will have the same structure
data_str <- list(
  T = sum(str_starts(names(real_cmr_data), "yr")),
  N = nrow(real_cmr_data),
  fc = real_cmr_data$fc,
  fc_code = real_cmr_data$fc_code
)

# simulate replicate data
n_reps <- 100
draws <- sample(1:4000, n_reps, replace = FALSE)
rep_data_list <- list()
for (i in 1:n_reps) {
  rep_data_list[[i]] <- sim_rep_cmr_data(fit, draws[i], data_str)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Test 3G.SR for transients ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set occasions and states for which we can calculate test components
occasions <- 2:(T-1)
states <- c(3, 7, 11) # adult states only

# for each replicate dataset, calculate table of discrepancy measures 
# for Test 3G.SR at occasions and states defined above (~ 10 mins)
stats_tables_3G.SR <- list()
for (r in 1:n_reps) {
  # format replicate data for get_3G.SR_table()
  rep_data <- rep_data_list[[r]]
  y <- rep_data$y
  y[y==13] <- 0
  colnames(y) <- str_c("yr", 2012 + 1:ncol(rep_data$y))
  rep_cmr_data <- as_tibble(y) %>%
    add_column(fc = data_str$fc)
  
  # initialise table of statistics for different occasions and states
  temp <- tibble(
    rep = r,
    occasion = rep(occasions, length(states)),
    state = rep(states, each = length(occasions)),
    ft_stat = NA
  )
  for (occ in occasions) {
    for (st in states) {
      observed <- get_3G.SR_table(rep_cmr_data, occasion = occ, state = st, nStates = 12)
      expected <- get_expected_frequencies(observed)
      temp$ft_stat[(temp$occasion == occ)&(temp$state == st)] <-
        get_ft(observed, expected)
    }
  }
  stats_tables_3G.SR[[r]] <- temp
}

# bind all tables into one
test_3G.SR_replicates <- bind_rows(stats_tables_3G.SR)

# for the real data, calculate table of discrepancy measures for Test 3G.SR 
# at occasions and states defined above
test_3G.SR_real <- tibble(
  rep = NA,
  occasion = rep(occasions, length(states)),
  state = rep(states, each = length(occasions)),
  ft_stat = NA
)
for (occ in occasions) {
  for (st in states) {
    observed <- get_3G.SR_table(real_cmr_data, occasion = occ, state = st, nStates = 12)
    expected <- get_expected_frequencies(observed)
    test_3G.SR_real$ft_stat[(test_3G.SR_real$occasion == occ) & (test_3G.SR_real$state == st)] <-
      get_ft(observed, expected)
  }
}

# Plot density of posterior replicate statistics and the real data statistics,
# faceted by occasion and state (code to add vertical lines to each facet from
# https://stackoverflow.com/questions/48593011/adding-grouped-geom-vline-to-multiple-facets)
test_3G.SR_plot <- test_3G.SR_replicates %>%
  ggplot(aes(x = ft_stat)) +
  theme_classic() +
  geom_density(fill="navyblue", alpha=0.3) +
  geom_vline(
    data = test_3G.SR_real, aes(xintercept = ft_stat), 
    colour = "red", linetype = "dashed"
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  facet_grid(occasion ~ state, scales = "free") +
  labs(
    title = "Components of Test 3G.SR for 'transience'",
    x = "Freeman-Tukey discrepancy measure"
  )
test_3G.SR_plot
# ggsave("figs/02c_ppc_test_3G.SR.png", test_3G.SR_plot, scale = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Test WBWA for memory ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set occasions and states for which we can calculate test components
occasions <- 2:(T-1)
states <- c(3, 7, 11) # adult states only

# for each replicate dataset, calculate table of discrepancy measures 
# for Test WBWA at occasions and states defined above (~ 18 mins)
stats_tables_WBWA <- list()
for (r in 1:n_reps) {
  # format replicate data for get_WBWA_table()
  rep_data <- rep_data_list[[r]]
  y <- rep_data$y
  y[y==13] <- 0
  colnames(y) <- str_c("yr", 2012 + 1:ncol(rep_data$y))
  rep_cmr_data <- as_tibble(y) %>%
    add_column(fc = data_str$fc)
  
  # initialise table of statistics for different occasions and states
  temp <- tibble(
    rep = r,
    occasion = rep(occasions, length(states)),
    state = rep(states, each = length(occasions)),
    ft_stat = NA
  )
  for (occ in occasions) {
    for (st in states) {
      observed <- get_WBWA_table(rep_cmr_data, occasion = occ, state = st, nStates = 12)
      expected <- get_expected_frequencies(observed)
      temp$ft_stat[(temp$occasion == occ)&(temp$state == st)] <-
        get_ft(observed, expected)
    }
  }
  stats_tables_WBWA[[r]] <- temp
}

# bind all tables into one
test_WBWA_replicates <- bind_rows(stats_tables_WBWA)

# for the real data, calculate table of discrepancy measures for Test WBWA 
# at occasions and states defined above
test_WBWA_real <- tibble(
  rep = NA,
  occasion = rep(occasions, length(states)),
  state = rep(states, each = length(occasions)),
  ft_stat = NA
)
for (occ in occasions) {
  for (st in states) {
    observed <- get_WBWA_table(real_cmr_data, occasion = occ, state = st, nStates = 12)
    expected <- get_expected_frequencies(observed)
    test_WBWA_real$ft_stat[(test_WBWA_real$occasion == occ) & (test_WBWA_real$state == st)] <-
      get_ft(observed, expected)
  }
}

# Plot density of posterior replicate statistics and the real data statistics,
# faceted by occasion and state 
test_WBWA_plot <- test_WBWA_replicates %>%
  ggplot(aes(x = ft_stat)) +
  theme_classic() +
  geom_density(fill="navyblue", alpha=0.3) +
  geom_vline(
    data = test_WBWA_real, aes(xintercept = ft_stat), 
    colour = "red", linetype = "dashed"
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  facet_grid(occasion ~ state, scales = "free") +
  labs(
    title = "Components of Test WBWA for 'memory'",
    x = "Freeman-Tukey discrepancy measure"
  )
test_WBWA_plot
# ggsave("figs/02c_ppc_test_WBWA.png", test_WBWA_plot, scale = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Test M.ITEC for immediate trap-dependence ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set occasions for which we can calculate a test statistic
occasions <- 2:(T-2)

# for each replicate dataset, calculate table of discrepancy measures 
# for Test M.ITEC at occasions defined above 
stats_tables_M.ITEC <- list()
for (r in 1:n_reps) {
  # format replicate data for get_M.ITEC_table()
  rep_data <- rep_data_list[[r]]
  y <- rep_data$y
  y[y==13] <- 0
  colnames(y) <- str_c("yr", 2012 + 1:ncol(rep_data$y))
  rep_cmr_data <- as_tibble(y) %>%
    add_column(fc = data_str$fc)

  # initialise table of statistics for different occasions
  temp <- tibble(
    rep = r,
    occasion = occasions,
    ft_stat = NA
  )
  for (occ in occasions) {
    # get observed tables and split into 'component' and 'mixture' parts
    observed <- get_M.ITEC_table(rep_cmr_data, occ, nStates = 12)
    obs_comp <- observed[str_starts(rownames(observed), "curr"), ]
    obs_mix <- observed[str_starts(rownames(observed), "prev"), ]
    
    # remove zero rows from obs_comp and obs_mix
    obs_comp <- obs_comp[rowSums(obs_comp) != 0, ]
    obs_mix <- obs_mix[rowSums(obs_mix) != 0, ]
    
    # estimate cell probabilities and expected cell frequences for mixture model 
    cell_probs <- get_cell_probs_and_mixture_weights(
      M = ncol(obs_comp), N_mix = nrow(obs_mix), N_comp <- nrow(obs_comp),
      Y_comp = obs_comp, Y_mix = obs_mix
    )
    expected <- get_expected_frequencies_mixtures(
      obs_comp = obs_comp, obs_mix = obs_mix, 
      cell_probs_comp = cell_probs$cell_probs_comp, 
      cell_probs_mix = cell_probs$cell_probs_mix 
    )
    
    # arrange observed and expected data in matrices such that the upper rows
    # containg the 'mixture' part and lower rows the 'component' part
    expected <- with(expected, rbind(expected_mix, expected_comp))
    observed <- rbind(obs_mix, obs_comp)
    if (!all(rownames(expected) == rownames(observed))) {
      stop("Expected and observed matrices do not correspond.")
    }
    
    # calculate freeman-tukey discrepancy
    temp$ft_stat[temp$occasion == occ] <- get_ft(observed, expected)
  }
  stats_tables_M.ITEC[[r]] <- temp
}

# bind all tables into one
test_M.ITEC_replicates <- bind_rows(stats_tables_M.ITEC)

# for the real data, calculate table of discrepancy measures for Test M.ITEC 
# at occasions and states defined above
test_M.ITEC_real <- tibble(
  rep = NA,
  occasion = occasions,
  ft_stat = NA
)
for (occ in occasions) {
  # get observed tables and split into 'component' and 'mixture' parts
  observed <- get_M.ITEC_table(real_cmr_data, occ, nStates = 12)
  obs_comp <- observed[str_starts(rownames(observed), "curr"), ]
  obs_mix <- observed[str_starts(rownames(observed), "prev"), ]
  
  # remove zero rows from obs_comp and obs_mix
  obs_comp <- obs_comp[rowSums(obs_comp) != 0, ]
  obs_mix <- obs_mix[rowSums(obs_mix) != 0, ]
  
  # estimate cell probabilities and expected cell frequences for mixture model 
  cell_probs <- get_cell_probs_and_mixture_weights(
    M = ncol(obs_comp), N_mix = nrow(obs_mix), N_comp <- nrow(obs_comp),
    Y_comp = obs_comp, Y_mix = obs_mix
  )
  expected <- get_expected_frequencies_mixtures(
    obs_comp = obs_comp, obs_mix = obs_mix, 
    cell_probs_comp = cell_probs$cell_probs_comp, 
    cell_probs_mix = cell_probs$cell_probs_mix 
  )
  
  # arrange observed and expected data in matrices such that the upper rows
  # containg the 'mixture' part and lower rows the 'component' part
  expected <- with(expected, rbind(expected_mix, expected_comp))
  observed <- rbind(obs_mix, obs_comp)
  if (!all(rownames(expected) == rownames(observed))) {
    stop("Expected and observed matrices do not correspond.")
  }
  
  # calculate freeman-tukey discrepancy
  test_M.ITEC_real$ft_stat[test_M.ITEC_real$occasion == occ] <- get_ft(observed, expected)
}

# Plot density of posterior replicate statistics and the real data statistics,
# faceted by occasion and state 
test_M.ITEC_plot <- test_M.ITEC_replicates %>%
  ggplot(aes(x = ft_stat)) +
  theme_classic() +
  geom_density(fill="navyblue", alpha=0.3) +
  geom_vline(
    data = test_M.ITEC_real, aes(xintercept = ft_stat), 
    colour = "red", linetype = "dashed"
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  facet_wrap(vars(occasion), scales = "free") +
  labs(
    title = "Components of Test M.ITEC for 'trap-dependence'",
    x = "Freeman-Tukey discrepancy measure"
  )
test_M.ITEC_plot
# ggsave("figs/02c_ppc_test_M.ITEC.png", test_M.ITEC_plot, scale = 1.5)
