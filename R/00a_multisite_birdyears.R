#' Assigning sites to birds in years where they were encountered at more than one
#' site (which we call 'multi-site bird-years'). 
#' The steps taken are:
#' (i) select focal colonies and filter out encounters elsewhere
#' (ii) keep only marking encounters or encounters between March and October
#' These preliminary steps ensure that every assigned site is in the focal 
#' colonies, and within the breeding season (or a marking encounter).
library(tidyverse)

# load encounters across the region west of Cape Agulhas
enc <- readRDS("data/encounter_wca.RDS")

# look at the sites present
table( enc$site )

# remove encounters outside of selected colonies
colonies <- c("Robben", "Boulders", "Stony")
enc <- filter( enc, site %in% colonies )

# keep only marking encounters and those between March and October
enc <- filter( enc, type=="M" | between(month, 3, 10) )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             ---- Multi-site bird-years ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Note! it is only a small percentage of bird-years in which a bird was
#' recorded at more than one site:
enc %>%
  group_by( primary_id, year ) %>%
  summarise( n_sites = length(unique(site)) ) %>% 
  ungroup() %>%
  group_by( n_sites ) %>% 
  summarise( n_birdyears = n() )
# 219 / (219 + 14677) ~ 0.014, less than 2%

# gather all the multi-site bird-years
multisite_birdyears <- enc %>%
  group_by( primary_id, year ) %>%
  reframe( unique_sites = unique(site) ) %>%  # returns ungrouped df (see ?reframe)
  group_by( primary_id, year ) %>%
  mutate( n_sites = length(unique_sites) ) %>% # constant within a bird-year
  ungroup() %>%
  filter( n_sites>1 )

# gather primary ids of birds with multi-site bird-years
pids_msby <- unique(multisite_birdyears$primary_id)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        ---- Assigning a single site to a multi-site bird-year ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' We assign a single site to a multi-site bird-year according to the 
#' following hierarchy:
#' (i) if the bird has ever been confirmed breeding at one of the bird-year 
#' sites, we choose that site, otherwise,
#' (ii) if the bird has been recorded at one of the bird-year sites more than
#' any other site in its life, we choose that site, otherwise,
#' (iii) we assign randomly from the bird-year sites.

# find breeding site(s), if any, for birds with multi-site bird-years
breeding_sites <- enc %>%
  filter(
    primary_id %in% pids_msby,
    activity %in% c("at nest", "breeding")
  ) %>%
  group_by( primary_id ) %>%
  reframe( br_site = unique(site) ) %>%  # returns ungrouped df (see ? reframe ) 
  group_by( primary_id ) %>%
  mutate( n_br_sites = length(br_site) ) %>%
  ungroup()

# if any bird has more than one breeding site we'd need to choose one, but 
# in this data none does - uncomment and run following to confirm:
# all( breeding_sites$n_br_sites == 1) # should be TRUE

# add breeding site information to 'multisite_birdyears'
multisite_birdyears <- left_join(
  multisite_birdyears, select(breeding_sites, -n_br_sites), by="primary_id"
)

# find 'preferred' site(s) for birds with multi-site bird-years, defined as
# sites (if any) at which they were encountered in more years than any other
years_by_site <- enc %>%
  filter( primary_id %in% pids_msby ) %>%
  select( primary_id, year, site ) %>% 
  distinct() %>%
  group_by( primary_id, site ) %>%
  summarise( n_years = n() ) %>%
  ungroup() 

preferred_sites <- years_by_site %>%
  group_by( primary_id ) %>%
  mutate( 
    max_years = n_years == max(n_years),   
    n_max = sum(max_years),              
    pref_site = ifelse( (n_max==1)&(max_years==TRUE), site, NA )
  ) %>%
  ungroup() %>%
  filter( !is.na(pref_site) ) %>%
  select(primary_id, pref_site) 
  
# add preferred site information to multisite_birdyears
multisite_birdyears <- left_join(
  multisite_birdyears, preferred_sites, by="primary_id"
)

# now use the info on breeding sites and preferred sites to assign a single site
# to each multi-site bird-year according to the hierarchy described above
# (see code to test this by inspection below, after saving)
multisite_birdyears <- multisite_birdyears %>%
  group_by( primary_id, year ) %>%
  mutate(
    assigned_site =
      case_when(
        any(!is.na(br_site)) & any(unique_sites==br_site) ~ br_site,   
        any(!is.na(pref_site)) & any(unique_sites==pref_site) ~ pref_site,
        .default = sample(unique_sites,1),
      )
  ) %>%
  ungroup() 

# to check that every bird-year has one and only one assigned site; uncomment 
# and run following code, should be TRUE
# multisite_birdyears %>%
#   select( primary_id, year, assigned_site ) %>%
#   distinct() %>%
#   group_by( primary_id, year ) %>%
#   summarise( exactly_one_assigned_site = n()==1 ) %>%
#   ungroup() %>%
#   pull( exactly_one_assigned_site ) %>%
#   all()

# save multisite_birdyears
# saveRDS(multisite_birdyears, "outputs/00a_multisite_birdyears.RDS")

# inspect the output of case_when() code used to create 'assigned_site' in
# a small example:
# test_msby <- tribble(
#   ~primary_id, ~year, ~unique_sites, ~br_site, ~pref_site,
#   "001", 2017, "Stony", NA, NA,
#   "001", 2017, "Boulders", NA, NA,
#   "002", 2013, "Stony", "Robben", NA,
#   "002", 2013, "Boulders", "Robben", NA,
#   "003", 2015, "Stony", NA, "Robben", 
#   "003", 2015, "Boulders", NA, "Robben", 
#   "004", 2016, "Stony", "Robben", "Robben", 
#   "004", 2016, "Boulders", "Robben", "Robben", 
#   "005", 2016, "Stony", "Robben", "Stony", 
#   "005", 2016, "Boulders", "Robben", "Stony", 
#   "006", 2020, "Stony", "Stony", NA,
#   "006", 2020, "Boulders", "Stony", NA,
#   "007", 2019, "Stony", NA, "Boulders",
#   "007", 2019, "Boulders", NA, "Boulders",
#   "008", 2022, "Stony", "Stony", "Boulders",
#   "008", 2022, "Boulders", "Stony", "Boulders"
# )
# test_msby %>%
#   group_by( primary_id, year ) %>%
#   mutate(
#     assigned_site =
#       case_when(
#         any(!is.na(br_site)) & any(unique_sites==br_site) ~ br_site,
#         any(!is.na(pref_site)) & any(unique_sites==pref_site) ~ pref_site,
#         .default = str_c("Randomly chose ", sample(unique_sites,1))
#       )
#   ) %>%
#   ungroup()
