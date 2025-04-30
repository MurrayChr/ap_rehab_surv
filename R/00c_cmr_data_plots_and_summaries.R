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
    legend.position.inside = TRUE,
    legend.position = c(0.33, 0.8)
  ) +
  facet_grid(
    hr ~ marking_site, 
    labeller = labeller(
      hr = as_labeller(hr_names),
      marking_site = as_labeller(site_names)
    )
  ) +
  labs(x = "Year", y = "Number of newly marked birds")
# ggsave("figs/00c_number_newly_marked_birds.png", scale = 1.2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- How many wild-raised birds marked as chicks were observed as adults? ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract relevant capture histories
y <- cmr_data %>%
  filter(
    marking_age == 1,
    hr == 0
  ) %>%
  select(starts_with("yr")) %>%
  as.matrix()

# first and last capture
fc <- apply(y, 1, function(x) min(which(x!=0)))
lc <- apply(y, 1, function(x) max(which(x!=0)))

# first-captured as juveniles, so fc + 2 is youngest adult occasion
sum(fc + 2 <= lc) # number detected as adults
100*mean(fc + 2 <= lc)/length(fc) # as a proportion of total




