# Create mark-recapture data sets from encounter data
library(tidyverse)

# load encounter data
enc <- readRDS("data/encounter_wca.RDS")

# select the colonies
colonies <- c("Robben", "Boulders", "Stony")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Identify and subset birds 'belonging' to colonies ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify all birds 'belonging' to selected colonies because they were either 
# (i) released at one of the colonies after being rehabilitated as chicks or
# (ii) marked at one of the colonies 

# filter marking encounters of birds belonging to selected colonies 
mark <- enc %>% 
  filter(
    type=="M",
    loc %in% colonies | ( str_detect(loc,"Rehab") & rhb_rls %in% colonies )
  )

# remove one bird marked as 'J' (juvenile)
mark <- filter( mark, age!="J" )

# extract the primary id's of these birds
pids <- unique( mark$primary_id ) # a few birds marked > 1, hence unique() call 

# filter encounters to these birds, and only at the selected colonies 
# (using 'site' variable)
enc_colonies <- enc %>%
  filter( 
    primary_id %in% pids,
    site %in% colonies
  )

# lump the "RW" code for 'resighted in the wild' under "HR" for handheld-reader
enc_colonies <- enc_colonies %>%
  mutate(
    type = ifelse( type=="RW", "HR", type)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              ---- Create mark-recapture datasets ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' TODO:
#' (i) make the function able to account for different years i.e. don't hard-code
#' years <- 2013:2023
#' (ii) ditto for the first and last months in which to bin encounters
#' (ii) this is quite a complicated function, could we have some testing of it?

# Function to create mark-recapture dataset from dataframe of encounters
get_cmr_data <- function( enc_df ) {
  # unique primary ids sorted by marking year
  pid <- enc_df %>%
    filter( type == "M" ) %>%
    arrange( year ) %>%
    pull( primary_id ) %>%
    unique()
  nbirds <- length( unique( enc_df$primary_id ) )
  if (length(pid) != nbirds) {
    stop("Some birds lack marking encounters.")
  }
  
  # study years
  years <- 2013:2023
  nyears <- length( years )
  
  # initialise tibble for capture histories
  cmr_data <- matrix(NA, nrow = nbirds, ncol = nyears) %>%
    as_tibble(.name_repair = ~ str_c("yr", years)) %>%
    add_column( primary_id = pid ) 
  
  # add age at (first) marking
  age_at_marking <- enc_df %>%
    filter( type == "M" ) %>%
    group_by( primary_id ) %>%
    arrange( year ) %>%
    filter( row_number()==1 ) %>%
    ungroup() %>%
    mutate(
      marking_age = case_when(
        age %in% c("C","B") ~ 1,
        age == "J" ~ 2,  # there shouldn't be any birds marked as 'J' here
        age == "A" ~ 3
      )
    ) %>%
    select( primary_id, marking_age ) 
  cmr_data <- left_join( cmr_data, age_at_marking, by = "primary_id" )
  
  # takes a few minutes
  for (i in 1:nrow(cmr_data)) {
    pid_i <- cmr_data$primary_id[i]
    det_by_year <- enc_df %>%
      filter( primary_id == pid_i ) %>%
      mutate( 
        encounter_in_season = between(month, 3, 10),
        eligible_encounter = encounter_in_season | type == "M"
      ) %>%  
      group_by( year )  %>%
      summarise( det = any( eligible_encounter ) ) %>%
      ungroup()
    det_all_years <- left_join( tibble(year = years), det_by_year, by = "year") %>%
      mutate( det = ifelse(is.na(det), FALSE, det) )  # change NAs to FALSE
    cmr_data[i,1:nyears] <-  as.list( ifelse(det_all_years$det, 1, 0) ) 
  }
  
  # add first capture
  cmr_data <- cmr_data %>%
    add_column(
      fc = apply( cmr_data, 1, function(x){ min(which(x[1:nyears]>0)) } ) 
    )
  
  # return
  cmr_data
}
