# Compare estimates of hand-rearing effects across models
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    ---- Load fitted models ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m_hr <- readRDS("outputs/04a_td_tr_hr_fit.RDS")
m_hr_imm <- readRDS("outputs/05a_td_tr_hr_phi-im_fit.RDS")
m_hr_imm_bss <- readRDS("outputs/06a_td_tr-bss_hr_phi-im_fit.RDS")
m_stony <- readRDS("outputs/07a_stony_td_hr_fit.RDS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Logit-scale plot: Extract posterior draws and plot ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prep_draws <- function(model, model_name) {
  model$draws(c("hr_ad", "hr_jv"), format = "df") %>%
    select(-starts_with(".")) %>%
    pivot_longer(everything(), names_to = "par") %>%
    add_column(model = model_name)
}

model_list <- list(
  "m_astony" = m_stony,
  "m_hr" = m_hr,
  "m_hr_imm" = m_hr_imm,
  "m_hr_imm_bss" = m_hr_imm_bss
)

post <- NULL
for (i in 1:length(model_list)) {
  post <- rbind(post, prep_draws(model_list[[i]], names(model_list)[i]))
}

model_colours <- c(
  "grey80", viridisLite::viridis(3, alpha = 0.7, begin = 0.1, end = 1, direction = -1)
)
model_labels <- c(
  "m_astony" = expression(M[Stony]^{hr}),
  "m_hr" = expression(M[td+tr]^{hr}),
  "m_hr_imm" = expression(M[td+tr+imm]^{hr}), 
  "m_hr_imm_bss" = expression(M[td+tr+imm+bss]^{hr})
)
facet_names <- c("hr_ad" = "Adults", "hr_jv" = "Juveniles")
logit_plot <- post %>%
  ggplot(aes(x = value, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  stat_density(linewidth = 1.2, geom = "line", position = "identity") +
  theme_classic() +
  scale_x_continuous(breaks = seq(-1,1 , length.out = 9)) +
  scale_colour_manual(labels = model_labels, values = model_colours) +
  facet_wrap(vars(par), labeller = as_labeller(facet_names)) +
  guides(colour = guide_legend(title = "Model")) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11),
    aspect.ratio = 0.7,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8) 
  ) +
  labs(
    y = "Posterior density",
    x = "Hand-rearing differences (logit scale)",
  )
logit_plot
# ggsave("figs/22_logit_hr_differences_compared.png", plot = logit_plot, scale = 1.4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Probability-scale plot: Extract posterior draws and plot ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Our model stipulates that hand-rearing differences are constant (in time and space) 
#' on the logit scale, but they are not constant on the probability scale. Therefore 
#' we average over space and time on the probability scale as follows:
#' For each posterior draw of phi_wr and phi_hr (these being site x year matrices),
#' we calculate the element-wise difference (i.e. the differences at corresponding 
#' sites and years), then average over the matrix of differences to obtain a single 
#' draw of the 'probability scale hand-rearing difference'. 

prep_draws_prob_scale <- function(model, model_name) {
  # extract posterior draws in draw x variable matrix
  model_draws <- model$draws(
    c("phi_ad_wr", "phi_ad_hr", "phi_jv_wr", "phi_jv_hr"),
    format = "df"
  ) %>%
    # create visible 'draw' variable, remove hidden variables, pivot longer
    mutate(draw = .draw) %>%
    select(-starts_with(".")) %>%
    pivot_longer(-"draw")
  #' check if this model is the Stony model or one of the multi-site models by 
  #' checking for the comma that separates the site and time indices in the variable 
  #' names of multi-site models 
  single_site <- !str_detect(model_draws$name[1], ",")
  if (single_site) {
    # pivot longer and extract age, time and hr status info from variable name
    model_draws <- model_draws %>%
      mutate(
        age = str_extract(name, "(?<=_)[a-z]+(?=_)"),
        t = str_extract(name, "(?<=\\[)[0-9]+(?=\\])"),
        hr_status = str_extract(name,"(?<=_)[a-z]+(?=\\[)")
      ) %>%
      select(-name) %>%
      # for each age, draw and t there is one row for hr and one for wr,
      # we want to pivot wider to calculate the difference between these
      group_by(age, draw, t) %>%
      pivot_wider(values_from = value, names_from = hr_status) %>%
      ungroup()
  }
  if (!single_site) {
    # extract age, site, time and hr status info from variable name
    model_draws <- model_draws %>%
      mutate(
        age = str_extract(name, "(?<=_)[a-z]+(?=_)"),
        site = str_extract(name, "(?<=\\[)[0-9]+(?=,)"),
        t = str_extract(name, "(?<=,)[0-9]+(?=\\])"),
        hr_status = str_extract(name,"(?<=_)[a-z]+(?=\\[)")
      ) %>%
      select(-name) %>%
      # for each age, draw, site and t there is one row for hr and one for wr,
      # we want to pivot wider to calculate the difference between these
      group_by(age, draw, site, t) %>%
      pivot_wider(values_from = value, names_from = hr_status) %>%
      ungroup()
  }
  # calculate hr differences and average over time (and site, if applicable), 
  # for each draw and age
  model_draws %>%
    mutate(diff = hr - wr) %>%
    group_by(draw, age) %>%
    summarise(
      mean_diff = mean(diff)
    ) %>%
    ungroup() %>%
    # drop draw info 
    select(-draw) %>%
    # add model name
    add_column(model = model_name)
}

model_list <- list(
  "m_astony" = m_stony,
  "m_hr" = m_hr,
  "m_hr_imm" = m_hr_imm,
  "m_hr_imm_bss" = m_hr_imm_bss
)

post <- NULL
for (i in 1:length(model_list)) {
  post <- rbind(post, prep_draws_prob_scale(model_list[[i]], names(model_list)[i]))
}

model_colours <- c(
  "grey80", viridisLite::viridis(3, alpha = 0.7, begin = 0.1, end = 1, direction = -1)
)
model_labels <- c(
  "m_astony" = expression(M[Stony]^{hr}),
  "m_hr" = expression(M[td+tr]^{hr}),
  "m_hr_imm" = expression(M[td+tr+imm]^{hr}), 
  "m_hr_imm_bss" = expression(M[td+tr+imm+bss]^{hr})
)
facet_names <- c("ad" = "Adults", "jv" = "Juveniles")
prob_plot <- post %>%
  ggplot(aes(x = mean_diff, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  stat_density(linewidth = 1.2, geom = "line", position = "identity") +
  theme_classic() +
  scale_x_continuous(breaks = seq(-0.3,0.3 , length.out = 13)) +
  scale_colour_manual(labels = model_labels, values = model_colours) +
  facet_wrap(vars(age), labeller = as_labeller(facet_names)) +
  guides(colour = guide_legend(title = "Model")) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11),
    aspect.ratio = 0.7,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8) 
  ) +
  labs(
    y = "Posterior density",
    x = "Hand-rearing differences (probability scale)",
  )
prob_plot
# ggsave("figs/22_prob_hr_differences_compared.png", plot = prob_plot, scale = 2)

# Calculate summaries for each posterior
post %>%
  group_by(age, model) %>%
  summarise(
    median = median(mean_diff),
    q5 = quantile(mean_diff, 0.05),
    q95 = quantile(mean_diff, 0.95)
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 ---- Combine plots and save ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(cowplot)

plot_grid(logit_plot, prob_plot, nrow = 2)
ggsave("figs/22_hr_differences_compared.png", scale = 2)
