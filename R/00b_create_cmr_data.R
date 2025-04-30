#' Create mark-recapture data set from encounter data.
#' The mark-recapture dataset created here is encoded for a multi-site, 
#' multi-age model with trap-dependence in the adult state.
library(tidyverse)

# load encounter data across the region west of Cape Agulhas
enc_wca <- readRDS("data/encounter_wca.RDS")

# select the colonies
colonies <- c("Robben", "Boulders", "Stony")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Identify and subset birds 'belonging' to colonies ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify all birds 'belonging' to selected colonies because they were either 
# (i) released at one of the colonies after being rehabilitated as chicks or
# (ii) marked at one of the colonies 

# get primary_id of all birds belonging to selected colonies 
pids <- enc_wca %>%
  filter(
    type=="M",
    loc %in% colonies | ( str_detect(loc,"Rehab") & rhb_rls %in% colonies ),
  ) %>%
  filter(
    age!="J"     # remove one bird marked as 'J' (juvenile) 
  ) %>%
  pull(primary_id) %>%
  unique()

# filter encounters to these birds, and only at the selected colonies 
# (using 'site' variable not 'loc')
enc <- enc_wca %>%
  filter( 
    primary_id %in% pids,
    site %in% colonies
  )

# filter encounters between March and October (unless a marking encounter,
# which we allow from any month)
enc <- enc %>%
  filter( type == "M" | between(month, 3, 10) )

# lump the "RW" code for 'resighted in the wild' under "HR" for handheld-reader
enc <- enc %>%
  mutate(
    type = ifelse( type=="RW", "HR", type)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    ---- Create tibble of bird-years with site, age and marking info ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Step 1: add assigned site for years in which a bird was encountered at more 
## than one site (see R/00a_multisite_birdyears.R)

# load multi-site bird-years and keep distinct bird-years only
multisite_birdyears <- readRDS("outputs/00a_multisite_birdyears.RDS")
multisite_birdyears <- multisite_birdyears %>%
  select("primary_id", "year", "assigned_site") %>%
  distinct()

# join to encounters
enc <- left_join(
  enc, multisite_birdyears, by = c("primary_id", "year")
)

# modify assigned site to equal site whenever an encounter is not in a 
# multi-site bird-year
enc <- enc %>%
  mutate( assigned_site = ifelse( is.na(assigned_site), site, assigned_site ) )

## Step 2: create age-at-marking tibble with the age at earliest marking of 
## every bird

# note that a few birds have more than one mark:
enc %>%
  filter( type=="M" ) %>%
  group_by( primary_id ) %>%
  summarise( n_marks = n() ) %>%
  ungroup() %>%
  group_by( n_marks ) %>%
  summarise( n_birds = n() ) %>%
  ungroup()

# create tibble with age-at-marking information
age_at_marking <- enc %>%
  filter( type == "M" ) %>%
  group_by( primary_id ) %>%
  arrange( year ) %>%             # this and next line filters for earliest...
  filter( row_number()==1 ) %>%   # ...encounter of each bird
  ungroup() %>%
  mutate(
    marking_age = case_when(
      age %in% c("C","B") ~ 1,
      age == "J" ~ 2,  # there shouldn't be any birds marked as 'J' here
      age == "A" ~ 3
    )
  ) %>%
  select( primary_id, marking_age ) 

# all birds should have marking_age available; uncomment and run to confirm
# all(!is.na(age_at_marking$marking_age)) # should be TRUE

# join to encounters
enc <- left_join( enc, age_at_marking, by = "primary_id" )

## Step 3: create a tibble with one row for each bird-year, including assigned
## site, marking year and age-at-marking information
bird_year_site <- enc %>%
  select( primary_id, year, assigned_site, mark_year, marking_age ) %>%
  distinct()

# there should be as many rows in bird_year_site as their are unique bird-years 
# in enc; uncomment and run following to confirm - should be TRUE
# nrow(bird_year_site) == nrow( distinct( select( enc, primary_id, year ) ) )

# create age_class '1', '2', '3' in every year from marking age and mark year
bird_year_site <- bird_year_site %>%
  mutate(
    age_class = pmin( year-mark_year+marking_age, 3) 
  )

# add marking site to dataframe
marking_sites <- bird_year_site %>%
  filter( year==mark_year ) %>%
  select( primary_id, assigned_site ) %>%
  rename(
    marking_site = assigned_site
  )
bird_year_site <- left_join( bird_year_site, marking_sites, by="primary_id" ) %>%
  relocate( marking_site, .after = assigned_site )

# confirm that assigned_site and marking_site are the same in mark_year;
# uncomment and run following; should be TRUE
# bird_year_site %>%
#   filter( year==mark_year ) %>%
#   mutate( check = assigned_site==marking_site ) %>%
#   pull( check ) %>%
#   all()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Create multi-age, multi-site mark-recapture data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Here we create mark-recapture data for a multi-age, multi-site model (without
#' trap-dependence)
 
# encode age_class and site (assigned site) information 
bird_year_site <- bird_year_site %>%
  mutate(
    code = case_when(
      age_class==1 & assigned_site=="Robben" ~ 1,
      age_class==2 & assigned_site=="Robben" ~ 2,
      age_class==3 & assigned_site=="Robben" ~ 3,
      age_class==1 & assigned_site=="Boulders" ~ 4,
      age_class==2 & assigned_site=="Boulders" ~ 5,
      age_class==3 & assigned_site=="Boulders" ~ 6,
      age_class==1 & assigned_site=="Stony" ~ 7,
      age_class==2 & assigned_site=="Stony" ~ 8,
      age_class==3 & assigned_site=="Stony" ~ 9
    )
  )

# Note: this encoding is specific to the model we will use to analyse the data,
# for a different model, a different encoding will usually be required

# pivot wider to mark-recapture dataset 
cmr_data <- bird_year_site %>%
  select( primary_id, year, code, mark_year, marking_age, marking_site ) %>%
  group_by( primary_id ) %>%
  pivot_wider(
    names_from = year, names_prefix = "yr", names_sort = TRUE,
    values_from = code, values_fill = 0
  ) %>% 
  ungroup() %>% 
  arrange( mark_year, marking_age, marking_site, primary_id )

# add first capture 
cmr_data <- cmr_data %>%
  rowwise() %>%
  mutate( 
    fc = min( which( c_across( starts_with("yr") )!=0 ) ),
    .after = marking_site
  ) %>%
  ungroup()

# uncomment to check that mark_year and first capture agree; assuming first year 
# is 2013, this should return TRUE:
# all(cmr_data$fc + 2013 - 1 == cmr_data$mark_year)

# add code at first capture 
# (used when simulating replicate data for posterior predictive checks)
cmr_data <- cmr_data %>%
  rowwise() %>%
  mutate(
    fc_code = c_across(starts_with("yr"))[fc],
    .after=fc
  )

##  Add hand-rearing information, checks, and save 
##

peng <- readRDS("data/penguins.RDS")
cmr_data <- left_join( cmr_data, peng, by="primary_id" ) %>%
  relocate( hr, .after=marking_age )

# check everyone has hand-reared information; following should return TRUE
# all(!is.na(cmr_data$hr))

# check hand-reared birds are all birds marked as chicks; should return TRUE
# all( filter(cmr_data, hr==1)$marking_age == 1)

# save
# saveRDS(cmr_data, "data/00b_cmr_data_multisite_multiage.RDS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Create multi-age, multi-site, trap-dep mark-recapture data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Here we create mark-recapture data for a multi-age, multi-site model that
#' that incorporates trap-dependence in the adult age class at each colony.

# encode age_class and site (assigned site) information 
bird_year_site <- bird_year_site %>%
  mutate(
    code = case_when(
      age_class==1 & assigned_site=="Robben" ~ 1,
      age_class==2 & assigned_site=="Robben" ~ 2,
      age_class==3 & assigned_site=="Robben" ~ 3,
      age_class==1 & assigned_site=="Boulders" ~ 5,
      age_class==2 & assigned_site=="Boulders" ~ 6,
      age_class==3 & assigned_site=="Boulders" ~ 7,
      age_class==1 & assigned_site=="Stony" ~ 9,
      age_class==2 & assigned_site=="Stony" ~ 10,
      age_class==3 & assigned_site=="Stony" ~ 11
    )
  )

# Note 1: codes 4, 8, 12 are for adults in a trap 'unaware' state which, by 
# definition, they're never encountered in, so these codes don't appear above.

# Note 2: this encoding is specific to the model we will use to analyse the data,
# for a different model, a different encoding will usually be required

# pivot wider to mark-recapture dataset 
cmr_data <- bird_year_site %>%
  select( primary_id, year, code, mark_year, marking_age, marking_site ) %>%
  group_by( primary_id ) %>%
  pivot_wider(
    names_from = year, names_prefix = "yr", names_sort = TRUE,
    values_from = code, values_fill = 0
  ) %>% 
  ungroup() %>% 
  arrange( mark_year, marking_age, marking_site, primary_id )

# add first capture 
cmr_data <- cmr_data %>%
  rowwise() %>%
  mutate( 
    fc = min( which( c_across( starts_with("yr") )!=0 ) ),
    .after = marking_site
  ) %>%
  ungroup()

# uncomment to check that mark_year and first capture agree; assuming first year 
# is 2013, this should return TRUE:
# all(cmr_data$fc + 2013 - 1 == cmr_data$mark_year)

# add code at first capture 
# (used when simulating replicate data for posterior predictive checks)
cmr_data <- cmr_data %>%
  rowwise() %>%
  mutate(
    fc_code = c_across(starts_with("yr"))[fc],
    .after=fc
  )

##  Add hand-rearing information, checks, and save 
##

peng <- readRDS("data/penguins.RDS")
cmr_data <- left_join( cmr_data, peng, by="primary_id" ) %>%
  relocate( hr, .after=marking_age )

# check everyone has hand-reared information; following should return TRUE
# all(!is.na(cmr_data$hr))
 
# check hand-reared birds are all birds marked as chicks; should return TRUE
# all( filter(cmr_data, hr==1)$marking_age == 1)

# save
# saveRDS(cmr_data, "data/00b_cmr_data_multisite_multiage_trapdep.RDS")
