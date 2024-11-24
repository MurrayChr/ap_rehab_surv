# Functions to get contingency tables for tests 'WBWA', '3G.SR' and 'M.ITEC'
# see Choquet et al U-CARE 2.2 User's Manual for details

# Warning! 
# This code is written to work with "data/00b_cmr_data_multisite_multiage_trapdep.RDS"
# and may not work as expected with other data.
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- 3G.SR test of transience ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Considers all individuals captured on a given occasion in a given state, and
#' partitions them according to whether or not they were previously captured
#' ('new'/'old') and whether or not they were ever re-encountered at some
#' later occasion ('reenc'/'never'). 

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
    stop(str_c("Error: code ", nStates+1, " may have been used for non-encounters after first capture; ",
               "please investigate and recode these as zeros if necessary."))
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
  # return(list(observed = observed, summary = obs_summary))
  
  # return
  observed
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- M.ITEC test of trap-dependence ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Considers all individuals known to be alive at occasion i and i+1, then 
#' partitions them into whether or not they were encountered on occasion i, 
#' and whether or not they were re-encountered at occasion i+1. Since they are 
#' known to be alive on occasion i and i+1, those not encountered at i
#' must have been encountered previously ('current'/'previous') and those not
#' re-encountered at i+1 must have been re-encountered at a later stage ('immediate'/'later'). 
#' The four cases resulting from this partition are then further decomposed 
#' according to the state at most recent release and next re-encounter. 

get_M.ITEC_table <- function(cmr_data, occasion, nStates) {
  # checks
  T <- sum(str_starts(names(cmr_data),"yr"))
  if (!(occasion %in% 2:(T-2))) {  # least one 'later' occasion after 'occasion+1' required
    stop(str_c("Error: with ",T," years, occasion must be between 2 and ",T-2,"."))
  }
  codes <- cmr_data %>%
    select(starts_with("yr")) %>%
    unlist() %>%
    unique()
  if ((nStates+1) %in% codes) {
    stop(str_c("Error: code ", nStates+1, " may have been used for non-encounters after first capture; ",
               "please investigate and recode these as zeros if necessary."))
  }
  
  # add last capture if necessary
  if (!("lc" %in% colnames(cmr_data))) {
    cmr_data <- cmr_data %>%
      rowwise() %>%
      mutate(lc = max(which(c_across(starts_with("yr")) != 0))) %>%
      ungroup()
  }
  
  # filter for individuals known to be alive on occasions 'occasion' and
  # 'occasion + 1' 
  cmr_data <- cmr_data %>%
    filter((fc <= occasion) & (lc >= (occasion + 1)))

  if (nrow(cmr_data)>0) {
    # create variables for occasion of most recent release (on or before 'occasion')
    # and next re-encounter (on or after 'occasion+1')
    yrs_release <- str_c("yr", 2013:(2012 + occasion))
    yrs_next_reenc <- str_c("yr", (2012 + occasion + 1):(2012 + T))
    cmr_data <- cmr_data %>%
      rowwise() %>%
      mutate(
        occ_release = max(which(c_across(all_of(yrs_release))!=0)),
        occ_next_reenc = occasion + 
          min(which(c_across(all_of(yrs_next_reenc))!=0))
      ) %>%
      ungroup()
    
    # create variables for states at most recent release and next re-encounter
    cmr_data <- cmr_data %>%
      rowwise() %>%
      mutate(
        state_release = c_across(starts_with("yr"))[occ_release],
        state_next_reenc = c_across(starts_with("yr"))[occ_next_reenc]
      ) 
    
    # create variables for 'current'/'previous' release and 'immediate'/'later'
    # next re-encounter and summarise
    obs_summary <- cmr_data %>%
      mutate(
        current = ifelse(occ_release == occasion, 1, 0),
        immediate = ifelse(occ_next_reenc == (occasion + 1), 1, 0)
      ) %>%
      group_by(
        current, immediate, state_release, state_next_reenc
      ) %>%
      summarise(n = n(), .groups="drop")
    
    # add the (very many) missing cases
    zero_tib <- tibble(
      current = rep(c(0,0,1,1), each=nStates^2),
      immediate = rep(c(0,1,0,1), each=nStates^2),
      state_release = rep(rep(1:nStates, nStates), 4),
      state_next_reenc = rep(rep(1:nStates, each=nStates), 4),
      n = 0
    )
    obs_summary <- left_join(
      zero_tib, obs_summary, 
      by=c("current","immediate", "state_release", "state_next_reenc")) %>%
      mutate(n = ifelse(is.na(n.y), n.x, n.y)) %>%
      select(-contains("."))
    
    # order obs_summary so that it can be arranged in a table by row
    obs_summary <- obs_summary %>%
      group_by(current, state_release) %>%
      arrange(desc(immediate), state_next_reenc) %>%
      ungroup() %>%
      arrange(current, state_release) %>%
      relocate( state_release, .after = current)
    
    # format as table/matrix
    observed <- obs_summary %>%
      pull(n) %>%
      matrix(nrow = 2*nStates, ncol = 2*nStates, byrow = TRUE) 
    rownames(observed) <- c(str_c("prev.",1:nStates), str_c("curr.",1:nStates))
    colnames(observed) <- c(str_c("immed.",1:nStates), str_c("later.",1:nStates))
  }
  
  if (nrow(cmr_data)==0) {
    observed <- matrix(0, nrow = 2*nStates, ncol = 2*nStates) 
    rownames(observed) <- c(str_c("prev.",1:nStates), str_c("curr.",1:nStates))
    colnames(observed) <- c(str_c("immed.",1:nStates), str_c("later.",1:nStates))
  }
  
  # for testing/debugging
  # return(list(observed = observed, summary = obs_summary))
  
  # return
  observed
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- WBWA test of memory ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Considers all individuals seen on occasion i in state l that have been seen
#' previously and will be seen again, and decomposes by last encounter state
#' and next re-encounter state.

get_WBWA_table <- function(cmr_data, occasion, state, nStates) {
  # checks
  T <- sum(str_starts(colnames(cmr_data),"yr"))
  if (!(occasion %in% 2:(T-1))) {  
    stop(str_c("Error: with ",T," years, occasion must be between 2 and ",T-1,"."))
  }
  codes <- cmr_data %>%
    select(starts_with("yr")) %>%
    unlist() %>%
    unique()
  if ((nStates+1) %in% codes) {
    stop(str_c("Error: code ", nStates+1, " may have been used for non-encounters after first capture; ",
               "please investigate and recode these as zeros if necessary."))
  }
  
  # filter all individuals seen at given occasion in given state
  yr_occ <- str_c("yr", 2012 + occasion)
  cmr_data <- cmr_data %>%
    filter(.data[[yr_occ]] == state)
  
  if (nrow(cmr_data) > 0) {
    # add last capture if necessary
    if (!("lc" %in% colnames(cmr_data))) {
      cmr_data <- cmr_data %>%
        rowwise() %>%
        mutate(lc = max(which(c_across(starts_with("yr")) != 0))) %>%
        ungroup()
    }
    
    # filter for individuals seen before and after occasion
    cmr_data <- cmr_data %>%
      filter((fc < occasion) & (lc > occasion))
    
    if (nrow(cmr_data) > 0) {
      # create variables for occasion of most recent release (strictly before 'occasion')
      # and next re-encounter (strictly after 'occasion+1')
      yrs_release <- str_c("yr", 2013:(2012 + occasion - 1))
      yrs_next_reenc <- str_c("yr", (2012 + occasion + 1):(2012 + T))
      cmr_data <- cmr_data %>%
        rowwise() %>%
        mutate(
          occ_release = max(which(c_across(all_of(yrs_release))!=0)),
          occ_next_reenc = occasion + 
            min(which(c_across(all_of(yrs_next_reenc))!=0)) 
        ) %>%
        ungroup()
      
      # create variables for states at most recent release and next re-encounter
      cmr_data <- cmr_data %>%
        rowwise() %>%
        mutate(
          state_release = c_across(starts_with("yr"))[occ_release],
          state_next_reenc = c_across(starts_with("yr"))[occ_next_reenc]
        )
      
      # summarise
      obs_summary <- cmr_data %>%
        group_by(state_release, state_next_reenc) %>%
        summarise(n = n(), .groups="drop") %>%
        ungroup()
      
      # add many missing (state-release, state_next_reenc) cases
      zero_tib <- tibble(
        state_release = rep(1:12, each = 12),
        state_next_reenc = rep(1:12, 12),
        n = 0
      )
      
      obs_summary <- left_join(
        zero_tib, obs_summary, 
        by=c("state_release", "state_next_reenc")) %>%
        mutate(n = ifelse(is.na(n.y), n.x, n.y)) %>%
        select(-contains("."))
      
      # order obs_summary so that it can be arranged in a table by row
      obs_summary <- obs_summary %>%
        arrange(state_release, state_next_reenc) 
      
      # format as table/matrix
      observed <- obs_summary %>%
        pull(n) %>%
        matrix(nrow = nStates, ncol = nStates, byrow = TRUE) 
      rownames(observed) <- str_c("prev.",1:nStates)
      colnames(observed) <- str_c("next.",1:nStates)
    }
  }
  
  if (nrow(cmr_data) == 0) {
    observed <- matrix(0, nrow = nStates, ncol = nStates) 
    rownames(observed) <- str_c("prev.",1:nStates)
    colnames(observed) <- str_c("next.",1:nStates)
  }
  
  # return
  observed
}
