#' Compare estimates between models multi-age, multi-site, trap-dependence models
#' with and without transients
library(tidyverse)

# Load model fits
fit_01 <- readRDS("outputs/01a_multiage_multisite_trap-dep_fit.RDS")
fit_02 <- readRDS("outputs/02a_multiage_multisite_trap-dep_trans_fit.RDS")

# function to create summary from fit for given pars
get_summary <- function(fit, pars) {
  fit$summary(pars) %>%
    mutate(
      par = str_extract(variable, "[a-zA-Z_]+(?=\\[)"),
      site_index = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
      t = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") ),
      year = 2012 + t,
      site = case_when(
        site_index == 1 ~ "Robben",
        site_index == 2 ~ "Boulders",
        site_index == 3 ~ "Stony"
      )
    ) 
}

# create combined summaries for common parameters
pars <- c("phi_jv", "phi_ad", "p_im", "p_ad_A", "p_ad_U", "p_ad_N")
precis <- rbind(
  get_summary(fit_01, pars) %>% add_column(model = "td"),
  get_summary(fit_02, pars) %>% add_column(model = "tdtr")
)

# plot select parameters
plot_pars <- function(precis, pars_to_plot) {
  precis %>%
    filter(par %in% pars_to_plot) %>%
    mutate(
      xshift = 0.1*case_when(
        model == "td" ~ -1,
        model == "tdtr" ~ 1
      )
    ) %>%
    ggplot(aes(x = year + xshift, colour = model)) +
    geom_pointrange(aes(y = median, ymin = q5, ymax = q95)) +
    theme_classic() +
    theme(
      panel.grid.major = element_line()
    ) +
    scale_x_continuous(breaks = seq(2013, 2023, 2)) +
    coord_cartesian(ylim = c(0,1)) +
    facet_grid(par ~ site) +
    labs(x="year", y="")
}

# adult and juvenile survival
plot_pars(precis, c("phi_ad", "phi_jv"))
# ggsave("figs/02d_survival_estimates_01_vs_02.png", height = 5, width = 10)

# detection probabilities
plot_pars(precis, c("p_im", "p_ad_A", "p_ad_U", "p_ad_N"))
# ggsave("figs/02d_detection_estimates_01_vs_02.png", height = 8, width = 10)


 