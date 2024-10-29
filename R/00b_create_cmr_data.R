# Create mark-recapture data sets from encounter data
library(tidyverse)

# load encounter data
enc_wca <- readRDS("data/encounter_wca.RDS")

# select the colonies
colonies <- c("Robben", "Boulders", "Stony")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Identify and subset birds 'belonging' to colonies ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify all birds 'belonging' to selected colonies because they were either 
# (i) released at one of the colonies after being rehabilitated as chicks or
# (ii) marked at one of the colonies 

# filter marking encounters of birds belonging to selected colonies 
mark <- enc_wca %>% 
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
enc_colonies <- enc_wca %>%
  filter( 
    primary_id %in% pids,
    site %in% colonies
  )

# lump the "RW" code for 'resighted in the wild' under "HR" for handheld-reader
enc_colonies <- enc_colonies %>%
  mutate(
    type = ifelse( type=="RW", "HR", type)
  )


