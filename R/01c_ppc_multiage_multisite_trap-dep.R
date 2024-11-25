#' Posterior predictive checks of model fit for '01_multiage_multisite_trap-dep.stan' 
#' to the data 'data/00b_cmr_data_multisite_multiage_trapdep.RDS'
library(tidyverse)
source("R/01b_fn_sim_rep_data_multiage_multisite_trap-dep.R")
source("R/00_functions_get_gof_tables.R")
source("R/00_function_get_expected_frequencies.R")
source("R/00_function_get_ft.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Load real data and simulate replicate datasets ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data and model fit to data
real_cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")
fit <- readRDS("outputs/01a_multiage_multisite_trap-dep_fit.RDS")

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
states <- c(3,7,11) # adult states only

# for each replicate dataset, calculate table of discrepancy measures 
# for Test 3G.SR at occasions and states defined above
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
    observed <- get_3G.SR_table(rep_cmr_data, occasion = occ, state = st, nStates = 12)
    expected <- get_expected_frequencies(observed)
    test_3G.SR_real$ft_stat[(test_3G.SR_real$occasion == occ) & (test_3G.SR_real$state == st)] <-
      get_ft(observed, expected)
  }
}

# Plot density of posterior replicate statistics and the real data statistics,
# faceted by occasion and state
test_3G.SR_replicates %>%
  ggplot(aes(x = ft_stat)) +
  theme_classic() +
  geom_density(fill="navyblue", alpha=0.5) +
  geom_vline(data = test_3G.SR_real, aes(xintercept = ft_stat), colour = "red") +
  facet_grid(occasion ~ state, scales = "free")
# source for code to add vertical lines to each facet, see 
# https://stackoverflow.com/questions/48593011/adding-grouped-geom-vline-to-multiple-facets

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Test WBWA for memory ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Test M.ITEC for immediate trap-dependence ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~