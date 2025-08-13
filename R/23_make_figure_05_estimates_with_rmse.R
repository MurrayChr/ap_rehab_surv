# Plot survival estimates from model 05 with an indication of their reliability,
# as measured by root-mean-square-error in the simulations for the model
library(tidyverse)
library(cowplot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            ---- Load model fit and simulation results ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit <- readRDS("outputs/05a_td_tr_hr_phi-im_fit.RDS")
sim_results <- readRDS("outputs/sim_05/05c_survival_bias_coverage_rmse.RDS")

# combine real-data estimates and simulation results on bias/coverage/rmse 
pars <- unique(sim_results$var_name)
est <- fit$summary(pars) %>%
  select(variable, median, q5, q95) 
est <- left_join(est, sim_results, by = "variable") %>%
  rename(site_index = site)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Functions to plot estimates (modified from 05a) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set colour scheme for colonies
# tried to get these of equal luminance (darkness/lightness) so that changes in 
# transparency affect each equally (cf https://ggplot2-book.org/scales-colour.html)
colony_colours <- RColorBrewer::brewer.pal(4, "Set1")[2:4]

# plot survival that vary by site and year
plot_estimates_by_site_year <- function(
    est, 
    par, 
    rmse_min,   # rmse value at which transparency starts to increase as rmse increases
    rmse_max,   # rmse value at which transparency stops increasing as rmse increases
    alpha_min = 0.1,  # maximum transparency (achieved for rmse >= rmse_max) 
    alpha_max = 1  # minimum transparency (achieved for rmse < rmse_min) 
    ) {
  tib <- est %>%
    filter(var_name == par) %>%
    mutate(
      xshift = 0.1*case_when(
        site_index == 1 ~ -1,
        site_index == 2 ~ 0,
        site_index == 3 ~ 1
      ),
      site = case_when(
        site_index == 1 ~ "Robben",
        site_index == 2 ~ "Boulders",
        site_index == 3 ~ "Stony"
      )
    ) %>%
    mutate(
      rmse_capped = case_when(
        rmse > rmse_max ~ rmse_max,
        rmse < rmse_min ~ rmse_min,
        (rmse_min <= rmse) & (rmse <= rmse_max) ~ rmse
      ),
      .after = rmse
    )
  x_breaks <- 2012 + 1:max(tib$t)
  tib %>%
    ggplot( aes(x=year+xshift, colour=site) ) +
    geom_point(
      aes(y = median, alpha = rmse_capped), size = 3
    ) +
    geom_linerange(
      aes(ymin=q5, ymax=q95, alpha = rmse_capped), linewidth = 1
    ) +
    scale_x_continuous( breaks=x_breaks ) +
    scale_y_continuous( breaks=seq(0,1,length.out=6) ) +
    scale_colour_manual(values = colony_colours) +
    scale_alpha_continuous(
      limits = c(rmse_min, rmse_max),
      range = c(alpha_max, alpha_min)
    ) +
    guides(colour = guide_legend(title = "Colony")) +
    coord_cartesian( ylim=c(0,1) ) +
    theme_classic() +
    guides(alpha = "none") +
    theme(
      panel.grid.major = element_line(colour = "grey95"),
      legend.text = element_text(size = 10)
    ) +
    labs( x= "Year", y="Estimate")
}

# function to plot parameters that only vary by year (just phi-im)
plot_estimates_by_year <- function(
    est, 
    par, 
    rmse_min,   # rmse value at which transparency starts to increase as rmse increases
    rmse_max,   # rmse value at which transparency stops increasing as rmse increases
    alpha_min = 0.1,  # maximum transparency (achieved for rmse >= rmse_max) 
    alpha_max = 1  # minimum transparency (achieved for rmse < rmse_min) 
    ) {
  if (par != "phi_im") {
    stop("plot_estimates_by_year only accepts 'phi_im' for the 'par' argument")
  }
  tib <- est %>%
    filter(var_name == par) %>%
    mutate(
      rmse_capped = case_when(
        rmse > rmse_max ~ rmse_max,
        rmse < rmse_min ~ rmse_min,
        (rmse_min <= rmse) & (rmse <= rmse_max) ~ rmse
      ),
      .after = rmse
    )
  x_breaks <- 2012 + 1:max(tib$t)
  tib %>%
    ggplot( aes(x=year) ) +
    geom_pointrange( 
      aes(y=median, ymin=q5, ymax=q95, alpha = rmse_capped), 
      size=0.6, linewidth=0.8, colour = "navyblue"
    ) +
    scale_x_continuous( breaks=x_breaks ) +
    scale_y_continuous( breaks=seq(0,1,length.out=6) ) +
    scale_alpha_continuous(
      limits = c(rmse_min, rmse_max),
      range = c(alpha_max, alpha_min)
    ) +
    coord_cartesian( ylim=c(0,1) ) +
    theme_classic() +
    guides(alpha = "none") +
    theme(
      panel.grid.major = element_line(colour = "grey95"),
      legend.text = element_text(size = 10)
    ) + 
    labs( x= "Year", y="Estimate")
}

# make survival plots
(plt_phi_ad_wr <- plot_estimates_by_site_year(est, "phi_ad_wr", rmse_min = 0.05, rmse_max = 0.1) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.25)
    ) +
    labs(title = "Adult survival, wild-raised birds"))
(plt_phi_jv_wr <- plot_estimates_by_site_year(est, "phi_jv_wr", rmse_min = 0.05, rmse_max = 0.1) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.35, 0.8)
    ) +
    labs(title = "Juvenile survival, wild-raised birds"))
(plt_phi_im <- plot_estimates_by_year(est, "phi_im", rmse_min = 0.05, rmse_max = 0.1) +
    labs(title = "Immature survival, all birds and colonies"))

# make movement plot
(plt_mv <- fit$draws(format = "df") %>%
    select( starts_with("m_") ) %>%
    pivot_longer( everything(), names_to = "variable", values_to = "draw" ) %>%
    mutate(
      site_from_index = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
      site_to_index = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") ),
      site_from = case_when(
        site_from_index == 1 ~ "Robben",
        site_from_index == 2 ~ "Boulders",
        site_from_index == 3 ~ "Stony"
      ),
      site_to = case_when(
        site_to_index == 1 ~ "Robben",
        site_to_index == 2 ~ "Boulders",
        site_to_index == 3 ~ "Stony"
      ),
      age_class = str_extract(variable, "m_[a-z]+")
    ) %>% 
    ggplot( aes(x = draw, fill = age_class ) ) +
    geom_density( bounds = c(0,1), alpha = 0.7, linewidth = 0.2 ) +
    scale_x_continuous(breaks = seq(0,1,length.out=6)) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(colour = "grey95"),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.9, 0.5)
    ) +
    scale_fill_manual(
      values=c("m_ad"="navyblue", "m_jv"="skyblue"),
      labels=c("m_ad"="Adult", "m_jv"="Juv. & imm.")
    ) +
    coord_cartesian( xlim = c(0,1), ylim = c(0,30) ) +
    facet_grid(rows = vars(site_from), cols = vars(site_to), switch = "y" ) +
    labs(
      title = "Movement probabilities, all birds",
      x = "Probability",
      fill = "Age class"
    ))

# combine plots of survival and movement
plot_grid(
  plt_phi_ad_wr, plt_phi_jv_wr, plt_phi_im, plt_mv,
  nrow = 2, ncol = 2, labels = "auto"
)
# ggsave("figs/23_model_05_survival_movement_estimates.png", height = 8, width = 12)
