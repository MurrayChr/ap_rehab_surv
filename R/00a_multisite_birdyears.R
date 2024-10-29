# Resolving the 'site' variable for birds that were encountered at more
# than one site in a year
library(tidyverse)

# load encounters across the region west of Cape Agulhas
enc <- readRDS("data/encounter_wca.RDS")

# look at sites
table(enc$site)

# remove encounters outside of the region
enc <- filter( enc, site!="outside_wca" )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    ---- How many times were birds recorded at >1 site in a year? ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

