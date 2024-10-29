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
multisite_birdyears <- multisite_birdyears %>%
  group_by( primary_id, year ) %>%
  mutate(
    assigned_site =
      case_when(
        any(!is.na(br_site)) & any(unique_sites==br_site) ~ br_site,
        any(!is.na(pref_site)) & any(unique_sites==pref_site) ~ pref_site,
        all(is.na(br_site)) & all(is.na(pref_site)) ~ sample(unique_sites,1)
      )
  ) %>%
  ungroup() 

# inspect the output of case_when() code used to create 'assigned_site' in
# a small example:
# test_msby <- tibble(
#   primary_id = c("001", "001", "002", "002", "003", "003"),
#   year = c(2021, 2021, 2019, 2019, 2023, 2023),
#   unique_sites = rep(c("Stony", "Boulders"),3),
#   n_sites = rep(2,6),
#   br_site = c(NA, NA, "Boulders", "Boulders", NA, NA),
#   pref_site = c(rep("Stony",4), rep(NA,2))
# )
# test_msby %>%
#   group_by( primary_id, year ) %>%
#   mutate(
#     assigned_site =
#       case_when(
#         any(!is.na(br_site)) & any(unique_sites==br_site) ~ br_site,
#         any(!is.na(pref_site)) & any(unique_sites==pref_site) ~ pref_site,
#         all(is.na(br_site)) & all(is.na(pref_site)) ~ sample(unique_sites,1)
#       )
#   ) %>%
#   ungroup()

