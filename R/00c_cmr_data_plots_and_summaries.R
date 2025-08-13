# Data plots and summaries
library(tidyverse)

# Load the data
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Number of newly marked birds ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Number of newly marked individuals by year, site, age and hand-rearing status
newly_marked <- cmr_data %>%
  group_by(mark_year, marking_site, marking_age, hr) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  complete(mark_year, marking_site, marking_age, hr, fill = list(n = 0)) 

hr_names <- c("0" = "Wild-raised", "1" = "Hand-reared")
site_names <- c("Boulders" = "Boulders", "Robben" = "Robben", "Stony" = "Stony")
newly_marked %>%
  filter(!(hr == 1 & marking_age == 3)) %>%  # no hand-reared birds marked as adults
  ggplot(aes(x = mark_year, colour = as.factor(marking_age), y = n)) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1, alpha = 0.1) +
  scale_x_continuous(
    breaks = seq(2013, 2023, 2),
    labels = str_c("'",seq(13,23,2))
  ) +
  scale_colour_manual(
    values = c("1" = "blue", "3"= "black"),
    labels = c("1" = "juvenile", "3" = "adult"),
    name = "marking age"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(colour = "grey95"),
    legend.position.inside = c(0.3, 0.8),
    legend.position = "inside"
  ) +
  facet_grid(
    hr ~ marking_site, 
    labeller = labeller(
      hr = as_labeller(hr_names),
      marking_site = as_labeller(site_names)
    )
  ) +
  labs(x = "Year", y = "Number of newly marked birds")
# ggsave("figs/00c_number_newly_marked_birds.png", height = 4, width = 7)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- How many wild-raised birds marked as chicks were observed as adults? ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmr_data %>%
  filter(
    marking_age == 1,
    hr == 0
  ) %>%
  mutate(
    obs_as_ad = fc + 2 <= lc  # marked as chicks, so adult two years later
  ) %>%
  group_by(
    obs_as_ad
  ) %>%
  summarise(
    n = n()
  )
# so 413 / (413 + 1291) ~ 0.24, about a quarter of them

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      ---- Numbers known to be alive and in a given state each year ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Here we use the cmr data that is encoded for the multistate model without
# trap-dependence. 
rm(list = ls())
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage.RDS")

# decant
y <- cmr_data %>%
  select(starts_with("yr")) %>%
  as.matrix()
fc <- cmr_data$fc
lc <- cmr_data$lc
T <- ncol(y)
N <- nrow(y)

#' Strategy: make a states x years matrix, loop through individuals incrementing
#' the appropriate entries. Add also two partially observed states for known age
#' but unknown site

# function to get the age class ("juv", "imm", "ad") of an individual each year
get_ages <- function(
    h, # capture history - a row of 'y'
    f # integer - first capture in h
  ) {
  # checks
  stopifnot(f == min(which(h!=0)))
  stopifnot(!(h[f] %in% c(2,5,8)))
  
  # create age vector
  T <- length(h)
  age <- rep(NA, T)
  age[f] <- case_when(
    h[f] %in% c(1,4,7) ~ "juv",
    h[f] %in% c(3,6,9) ~ "ad"
  )
  if (f < T) {
    for (t in f:(T-1)) {
      age[t+1] <- case_when(
        age[t] == "juv" ~ "imm",
        age[t] == "imm" ~ "ad",
        age[t] == "ad" ~ "ad"
      )
    }
  }
  age
}

# inspect get_ages function
# i <- sample(1:N, 1)
# y[i,]
# get_ages(y[i,], fc[i])

# construct matrix of totals
total <- matrix(0, nrow = 9+2, ncol = T)
rownames(total) <- c(
  "Robben, juv",
  "Robben, imm",
  "Robben, ad",
  "Boulders, juv",
  "Boulders, imm",
  "Boulders, ad",
  "Stony, juv",
  "Stony, imm",
  "Stony, ad",
  "Unkown, imm",
  "Unknown, ad"
)
colnames(total) <- 2013:2024

for (i in 1:N) {
  if (any(y[i,fc[i]:lc[i]] == 0)) {
    age <- get_ages(y[i,],fc[i])
  }
  for (t in fc[i]:lc[i]) {
    if (y[i,t] != 0) {
      total[y[i,t], t] <- total[y[i,t], t] + 1
    }
    if (y[i,t] == 0) {
      row_index <- case_when(
        age[t] == "imm" ~ 10,
        age[t] == "ad" ~ 11
      )
      total[row_index, t] <- total[row_index, t] + 1
    }
  }
}
total

# check that the total number observed in known states is the same as a
# direct count from y; uncommenting and running the following should give TRUE
# all(rowSums(total)[1:9] == table(c(y))[str_c(1:9)])

## Repeat separately for hand-reared and wild-raised birds
hr <- cmr_data$hr
y_hr <- y[hr==1,]
y_wr <- y[hr==0,]
fc_hr <- fc[hr==1]
fc_wr <- fc[hr==0]
lc_hr <- lc[hr==1]
lc_wr <- lc[hr==0]

get_totals <- function(y, fc, lc) {
  N <- nrow(y)
  T <- ncol(y)
  total <- matrix(0, nrow = 9+2, ncol = T)
  rownames(total) <- c(
    "Robben, juv",
    "Robben, imm",
    "Robben, ad",
    "Boulders, juv",
    "Boulders, imm",
    "Boulders, ad",
    "Stony, juv",
    "Stony, imm",
    "Stony, ad",
    "Unkown, imm",
    "Unknown, ad"
  )
  colnames(total) <- 2013:2024
  for (i in 1:N) {
    if (any(y[i,fc[i]:lc[i]] == 0)) {
      age <- get_ages(y[i,],fc[i])
    }
    for (t in fc[i]:lc[i]) {
      if (y[i,t] != 0) {
        total[y[i,t], t] <- total[y[i,t], t] + 1
      }
      if (y[i,t] == 0) {
        row_index <- case_when(
          age[t] == "imm" ~ 10,
          age[t] == "ad" ~ 11
        )
        total[row_index, t] <- total[row_index, t] + 1
      }
    }
  }
  total
}  

(total_hr <- get_totals(y_hr, fc_hr, lc_hr))
(total_wr <- get_totals(y_wr, fc_wr, lc_wr))

total_hr_tib <- as_tibble(total_hr) %>%
  add_column(
    site = c(rep(c("Robben", "Boulders", "Stony"), each = 3),rep("Unknown", 2)),
    age = c(rep(c("Juvenile", "Immature", "Adult"), 3), c("Immature", "Adult")),
    hr = 1
  )
total_wr_tib <- as_tibble(total_wr) %>%
  add_column(
    site = c(rep(c("Robben", "Boulders", "Stony"), each = 3),rep("Unknown", 2)),
    age = c(rep(c("Juvenile", "Immature", "Adult"), 3), c("Immature", "Adult")),
    hr = 0
  )

# Plot the balance
balance_plt_colours <- viridisLite::mako(2, begin = 0.3, end = 0.7)
total_tib <- rbind(total_hr_tib, total_wr_tib)
total_tib %>%
  pivot_longer(-c("site","age","hr"), names_to = "year", values_to = "count") %>%
  mutate(
    signed_count = ifelse(hr==1, 1, -1)*count,
    year = as.integer(year)
  ) %>%
  ggplot(aes(x = year, fill = as.factor(hr))) +
  geom_col(aes(y = signed_count), position = "identity") +
  theme_classic() +
  coord_cartesian(ylim = c(-500, 500)) +
  scale_y_continuous(
    breaks = seq(-500,500,length.out = 5), 
    labels = abs(seq(-500,500,length.out = 5))
  ) +
  scale_x_continuous(
    # breaks = 2013:2024,
    # labels = str_c("'", 13:24),
    breaks = seq(2013, 2024, by = 2),
    labels = str_c("'", seq(13, 24, by = 2))
  ) +
  scale_fill_manual(
  # scale_fill_discrete(
    breaks = c("1", "0"),
    labels = c("1" = "hand-reared", "0" = "wild-raised"),
    # values = c("navyblue", "skyblue")
    values = balance_plt_colours
  ) +
  # scale_fill_viridis_d(name = "magma") +
  theme(
    panel.grid.major = element_line(linewidth = 0.25),
    legend.position.inside = TRUE,
    legend.position = c(0.875, 0.2),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    aspect.ratio = 0.7
  ) +
  labs(
    y = "Number of birds known to be alive",
    x = "Year"
  ) +
  facet_grid(age ~ site)
# ggsave("figs/00c_numbers_known_to_be_alive.png", scale = 1.8)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        ---- Observed lifespans ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Here we use the cmr data that is encoded for the multistate model without
# trap-dependence. 
rm(list = ls())
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage.RDS")

# what is the distribution of observed lifespans?
# we assume that birds marked as adults were >= 2 years when marked
cmr_data <- cmr_data %>%
  mutate(
    obs_ls = case_when(
      marking_age == 1 ~ lc - fc,
      marking_age == 3 ~ lc - fc + 2,
    ),
    .after = "lc"
  ) 
obs_lifespan <- cmr_data$obs_ls

# all individuals
barplot(table(obs_lifespan))

# individuals observed survive to 'adult' (age >= 2)
barplot(table(obs_lifespan[obs_lifespan >= 2]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  ---- Observations of immatures ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' What proportion of individuals marked as juveniles and known to survive to
#' adulthood were resighted as immatures?
rm(list = ls())

# again, use the data encoded for the multistate model without trap-dependence
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage.RDS")

# add observed lifespan
# we assume that birds marked as adults were >= 2 years when marked
cmr_data <- cmr_data %>%
  mutate(
    obs_ls = case_when(
      marking_age == 1 ~ lc - fc,
      marking_age == 3 ~ lc - fc + 2,
    ),
    .after = "lc"
  ) 

# filter for birds marked as juveniles and known to survive to adulthood,
# then query if they were seen the year after marking
target_cmr_data <- cmr_data %>%
  filter(
    marking_age == 1,
    obs_ls >= 2
  ) %>%
  rowwise() %>%
  mutate(
    fc_plus_one_code = c_across(starts_with("yr"))[fc + 1], # zero if not recaptured as immature, else non-zero
    .after=fc_code
  ) %>%
  ungroup() 

# number of birds marked as juveniles and known to survive to adulthood
nrow(target_cmr_data)

# proportion seen as immatures
mean(target_cmr_data$fc_plus_one_code != 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          ---- Adult age distributions in the Stony data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Here we compare the age distributions for adults between hand-reared and
#' wild-raised birds in the Stony dataset, which consists only of birds marked
#' as fledglings. 

# load, subset and prepare Stony data as in 07a_fit_stony_td_hr.R
# load data with states coded with trap-dependence
cmr_data <- readRDS("data/00b_cmr_data_multisite_multiage_trapdep.RDS")

# filter for Stony birds marked as juveniles
cmr_data <- cmr_data %>%
  filter(
    marking_site == "Stony",
    marking_age == 1
  ) 

# set number of years
T <- sum(str_starts(colnames(cmr_data),"yr"))

# extract capture history matrix
y <- cmr_data %>%
  select( starts_with("yr") ) %>%
  as.matrix()
fc <- cmr_data$fc
hr <- cmr_data$hr

# recode entries for get_marray()
# the codes in y are:
table(c(y))
# first, all codes relating to Robben and Boulders changed to zero
y[y %in% 1:8] <- 0
# second, codes for Stony mapped to 1:4
y[y==9] <- 1    # juvenile
y[y==10] <- 2   # immature
y[y==11] <- 3   # adult, trap-aware
y[y==12] <- 4   # adult, trap-unaware

# change zeros after first capture to nStates + 1 = 5
for (i in 1:nrow(y)) {
  y[i,fc[i]:T][y[i, fc[i]:T]==0] <- 5
}

# last capture occasion
lc <- apply(y, 1, function(x) (max(which(x!=5))))

# create matrix with age of each individual in the years they're known to be alive
ages <- matrix(NA, nrow = nrow(y), ncol = T)
for (i in 1:nrow(y)) {
  ages[i, fc[i]:lc[i]] <- (0:T)[1:(lc[i]-fc[i]+1)]
}

# summarise over individuals
ages_hr <- ages[hr == 1,]
ages_wr <- ages[hr == 0,]
age_distr_hr_all_years <- table(ages_hr, useNA = "no") %>%
  as_tibble() %>%
  rename(age = ages_hr) %>%
  mutate(age = as.integer(age)) %>%
  add_column(hr = "hr")
age_distr_wr_all_years <- table(ages_wr, useNA = "no") %>%
  as_tibble() %>%
  rename(age = ages_wr) %>%
  mutate(age = as.integer(age)) %>%
  add_column(hr = "wr")
age_distr_all_years <- rbind(age_distr_hr_all_years, age_distr_wr_all_years)

# create age distributions for adults pooled over all years
ad_age_distr_all_years <- age_distr_all_years %>%
  filter(age >= 1) %>%  # immature survival (1yo) lumped with adult survival (2+) in Stony model
  complete(age, hr, fill = list(n = 0)) # complete age x hr combinations

# add within-group proportions for each age
ad_age_distr_all_years <- ad_age_distr_all_years %>%
  group_by(hr) %>%
  mutate(
    p = n/sum(n), .after = n
  ) %>%
  ungroup() 

# plot (mirror)
# balance_plt_colours <- viridisLite::mako(2, begin = 0.3, end = 0.7)
# ad_age_distr_all_years %>%
#   mutate(signed_proportion = ifelse(hr == "hr", 1, -1)*p) %>%
#   ggplot(aes(x = age, y = signed_proportion, fill = as.factor(hr))) +
#   geom_col(position = "identity") +
#   scale_x_continuous(breaks = 1:10) + 
#   scale_y_continuous(
#     breaks = seq(-0.3, 0.3, length.out = 7),
#     labels = abs(seq(-0.3, 0.3, length.out = 7))
#   ) +
#   scale_fill_manual(
#     breaks = c("hr", "wr"),
#     labels = c("hr" = "hand-reared", "wr" = "wild-raised"),
#     values = balance_plt_colours
#   ) +
#   labs(
#     title = "'Adult' age distributions in Stony dataset, pooled over years",
#     y = "Proportion of 'adults' per group", x = "Age"
#   ) +
#   theme_classic() +
#   theme(
#     panel.grid.major = element_line(colour = "grey95"),
#     legend.position = "inside",
#     legend.position.inside = c(0.775,0.175),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 12)
#   )
  
# plot (dodged) - this one makes it easier to compare distributions
balance_plt_colours <- viridisLite::mako(2, begin = 0.3, end = 0.7)
ad_age_distr_all_years %>%
  mutate(signed_proportion = ifelse(hr == "hr", 1, -1)*p) %>%
  ggplot(aes(x = age, y = p, fill = as.factor(hr))) +
  geom_col(position = "dodge", alpha = 0.9) +
  scale_x_continuous(breaks = 1:10) + 
  scale_y_continuous(
    breaks = seq(-0.3, 0.3, length.out = 7)
  ) +
  scale_fill_manual(
    breaks = c("hr", "wr"),
    labels = c("hr" = "hand-reared", "wr" = "wild-raised"),
    values = balance_plt_colours
  ) +
  labs(
    title = "'Adult' age distributions in Stony dataset, pooled over years",
    y = "Proportion of 'adults' per group", x = "Age"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(colour = "grey95"),
    legend.position = "inside",
    legend.position.inside = c(0.775,0.475),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )
# ggsave("figs/00c_number_adult-imms_known_to_be_alive_stony.png", scale = 1.5)


