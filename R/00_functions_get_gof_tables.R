# Functions to get contingency tables for tests 'WBWA', '3G.SR' and 'M.ITEC'
# see Choquet et al U-CARE 2.2 User's Manual for details

# Warning! this code is written to work with 
# "data/00b_cmr_data_multisite_multiage_trapdep.RDS"
# and may not work as expected using if called on other data.
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- 3G.SR test of transience ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_3G.SR_table <- function(cmr_data, occasion, state, nStates) {
  # checks
  T <- sum(str_starts(names(cmr_data),"yr"))
  if (!(occasion %in% 2:(T-1))) { 
    stop(str_c("Error: with ",T," years, occasion must be between 2 and ",T-1,"."))
  }
  if (!(state %in% 1:nStates)) {
    stop(str_c("Error: only ",nStates," states available, but ",state," requested."))
  }
  codes <- cmr_data %>%
    select(starts_with("yr")) %>%
    unlist() %>%
    unique()
  if ((nStates+1) %in% codes) {
    stop(str_c("Error: code ", nStates+1, " may have been used for non-encounters after first capture;",
               "please recode these as zeros."))
  }
  
  # filter for individuals in encountered in given state on given occasion
  yr_occ <- str_c("yr", 2012+occasion)
  cmr_data <- cmr_data %>%
    filter(
      .data[[yr_occ]] == state
    )
  
  if (nrow(cmr_data) > 0) {
    # add last capture if necessary
    if (!("lc" %in% colnames(cmr_data))) {
      cmr_data <- cmr_data %>%
        rowwise() %>%
        mutate(lc = max(which(c_across(starts_with("yr")) != 0))) %>%
        ungroup()
    }
    
    # count combinations of 'new'/'old', 'later re-encountered'/'never re-encountered'
    obs_summary <- cmr_data %>%
      mutate(
        new = ifelse(fc == occasion, 1, 0),
        reenc = ifelse(lc > occasion, 1, 0)
      ) %>%
      group_by(new, reenc) %>%
      summarise(n=n(), .groups = "drop") %>%
      ungroup() 
    
    # fill in missing cases if any
    zero_tib <- tribble(
      ~new, ~reenc, ~n,
      0, 0, 0,
      0, 1, 0,
      1, 0, 0,
      1, 1, 0
    )
    obs_summary <- left_join(zero_tib, obs_summary, by=c("new","reenc")) %>%
      mutate(n = ifelse(is.na(n.y), n.x, n.y)) %>%
      select(-contains("."))
  }
  
  if (nrow(cmr_data)==0) {
    obs_summary <- tribble(
      ~new, ~reenc, ~n,
      0, 0, 0,
      0, 1, 0,
      1, 0, 0,
      1, 1, 0
    )
  }
  
  # format as table/matrix
  observed <- obs_summary %>%
    arrange(desc(new), desc(reenc)) %>%
    pull(n) %>%
    matrix(nrow = 2, ncol = 2, byrow = TRUE) 
  rownames(observed) <- c("new", "old")
  colnames(observed) <- c("reenc", "never")
  
  # return list for testing/debugging
  # list(observed = observed, summary = obs_summary)
  
  # return
  observed
}

for (state in 3:3) {
  for (occasion in 10:12) {
    get_3G.SR_table(cmr_data, occasion, state, nStates = 12) %>%
      print()
  }
}







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- M.ITEC test of trap-dependence ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- WBWA test of memory ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


