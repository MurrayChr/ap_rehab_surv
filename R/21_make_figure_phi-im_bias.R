#' Make the composite figure illustrating the potential bias arising if 
#' immature survival were in fact lower than adult survival, but we assumed
#' they were equal, which could create a differential survival bias between 
#' hand-reared and wild-raised groups due to their composition (with more immatures
#' in the hand-reared group), which the model would presumably attribute to hand-rearing
#' differences
library(cmdstanr)
library(tidyverse)
library(svglite)

# set seed to reproduce figure
set.seed(2)

# Set some true parameter values to be used throughout
phi_jv <- 0.4
phi_im <- 0.6
phi_ad <- 0.8

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Confounding in juvenile and immature survival ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fit a mock model
N <- 200
stan_data <- list(N = N, y = rbinom(1, N, phi_jv*phi_im))
stan_file <- "
data {
  int N;
  array[1] int y;
}

parameters {
  real<lower=0, upper = 1> phi_im;
  real<lower=0, upper = 1> phi_jv;
}

model {
  phi_im ~ beta(2,2);
  phi_jv ~ beta(2,2);
  y ~ binomial(N, phi_jv*phi_im);
}
"
mod <- cmdstan_model(write_stan_file(stan_file))
fit <- mod$sample(stan_data, parallel_chains = 4)
post <- fit$draws(variables = c("phi_jv", "phi_im"), format = "df")

# now fit with biased information on phi_im
stan_file2 <- "
data {
  int N;
  array[1] int y;
}

parameters {
  real<lower=0, upper = 1> phi_im;
  real<lower=0, upper = 1> phi_jv;
}

model {
  phi_im ~ beta(70,30);
  phi_jv ~ beta(2,2);
  y ~ binomial(N, phi_jv*phi_im);
}
"
mod2 <- cmdstan_model(write_stan_file(stan_file2))
fit2 <- mod2$sample(stan_data, parallel_chains = 4)
post2 <- fit2$draws(variables = c("phi_jv", "phi_im"), format = "df")

#
# plot with just the confounded estimates for phi_jv and phi_im
#
(plot_a <- post %>%
  ggplot(aes(x = phi_im, y = phi_jv)) +
  geom_point(colour = "grey", alpha = 0.1) +
  #
  # code for vertical line and its annotation
  #
  geom_vline(xintercept = phi_im, linetype = "dashed", alpha = 0.3) +
  annotate(
    "text", x = phi_im - 0.1, y = 0.85, 
    label = expression(varphi[imm]^"truth"),
    size = 5, colour = "grey20"
  ) +
  #
  # code for horizontal line and it's annotation 
  #
  geom_hline(yintercept = phi_jv, linetype = "dashed", alpha = 0.3) +
  annotate(
    "text", x = 0.2, y = phi_jv + 0.075, 
    label = expression(varphi[juv]^"truth"),
    size = 5, colour = "grey20"
  ) +
  # 
  # code for global formatting
  #
  theme_classic() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = expression("Immature survival estimate"~hat(varphi)[imm]),
    y = expression("Juvenile survival estimate"~hat(varphi)[juv])
  ))
# ggsave("figs/21_phi-im_bias_plot_a.svg", plot_a, device = svg)

# plot with the estimates biased, assuming phi_imm = phi_ad if in fact
# phi_im < phi_ad
(plot_b <- post %>%
  ggplot(aes(x = phi_im, y = phi_jv)) +
  #
  # code for vertical lines and their annotations
  #
  geom_vline(xintercept = phi_im, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = phi_ad, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0.5*phi_im + 0.5*phi_ad, linetype = "dashed", 
             alpha = 0.4, colour = "darkred") +
  annotate(
    "text", x = phi_im - 0.1, y = 0.85, 
    label = expression(varphi[imm]^"truth"),
    size = 5, colour = "grey20"
  ) +
  annotate(
    "text", x = phi_ad + 0.1, y = 0.85, 
    label = expression(varphi[ad]^"truth"),
    size = 5, colour = "grey20"
  ) +
  annotate(
    "text", x = 0.885, y = 0.05,
    label = expression(hat(varphi)[imm]~"="~hat(varphi)[ad]),
    size = 5, colour = "darkred", alpha = 0.9
  ) +
  annotate(
    "segment", x = phi_im, y = 0.1, xend = phi_im + 0.09, yend = 0.1, 
    size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  #
  # code for horizontal lines and their annotataions 
  #
  geom_hline(yintercept = phi_jv, linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = phi_jv - 0.08, linetype = "dashed", alpha = 0.4, colour = "darkred") +
  annotate(
    "text", x = 0.2, y = phi_jv + 0.075, 
    label = expression(varphi[juv]^"truth"),
    size = 5, colour = "grey20"
  ) +
  annotate(
    "text", x = 0.2, y = phi_jv - 0.175,
    label = expression(hat(varphi)[juv]),
    size = 5, colour = "darkred", alpha = 0.9
  ) +
  annotate(
    "segment", x = 0.2, y = phi_jv, xend = 0.2, yend = phi_jv - 0.07, 
    size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  # 
  # code for the scatterplots of posterior draws and global formatting
  #
  geom_point(colour = "grey", alpha = 0.1) +
  geom_point(data = post2, colour = "darkred", alpha = 0.02) +
  theme_classic() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = expression("Immature survival estimate"~hat(varphi)[imm]),
    y = expression("Juvenile survival estimate"~hat(varphi)[juv])
  ))
# ggsave("figs/21_phi-im_bias_plot_b.svg", plot_b, device = svg)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              ---- Differential bias in adults ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(plot_d <- ggplot() +
   geom_point(
     data = tibble(x = runif(100, 0, 1.2), y = runif(100, 0, 1)), 
     mapping = aes(x = x, y = y),
     colour = "white"
   ) + 
   #
   # global formatting
   #
   theme_classic() +
   theme(
     # aspect.ratio = 0.7/1.2,
     # aspect.ratio = 0.5,
     aspect.ratio = 0.29,
     axis.ticks = element_blank(),
     axis.text = element_blank()
   ) +
   coord_cartesian(ylim = c(0.5, 0.9)) +
   labs(x = "Hand-rearing status", y = "Survival probability") +
   # 
   # code for 'hand-reared' column
   #
   annotate(
     "text", x = 0.35, y = phi_ad + 0.1, label = "Hand-reared"
   ) +
   # true adult survival
   geom_line(
     data = tibble(x = c(0.2, 0.5), y = rep(phi_ad, 2)),
     mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
   ) +  
   annotate(
     "text", x = 0.125, y = phi_ad, size = 4,
     label = expression(varphi[ad]^"truth")
   ) +
   # true immature survival
   geom_line(
     data = tibble(x = c(0.2, 0.5), y = rep(phi_im, 2)),
     mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
   ) +
   annotate(
     "text", x = 0.125, y = phi_im, size = 4,
     label = expression(varphi[imm]^"truth")
   ) +
   # estimate phi_im = phi_ad
   geom_line(
     data = tibble(x = c(0.2, 0.7), y = rep(0.6*phi_im + 0.4*phi_ad, 2)),
     mapping = aes(x = x, y = y), linetype = "dashed", colour = "darkred"
   ) +
   annotate(
     "text", x = 0.075, y = 0.6*phi_im + 0.4*phi_ad + 0.02, size = 4,
     label = expression(hat(varphi)[imm]^"hr"~"="~hat(varphi)[ad]^"hr"),
     colour = "darkred"
   ) +
   # arrow for bias in phi_im
   # annotate(
   #   "segment", x = 0.35, y = phi_im, xend = 0.35, yend = 0.6*phi_im + 0.4*phi_ad - 0.01, 
   #   size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
   #   arrow = arrow(type = "closed", length = unit(0.015, "npc"))
   # ) +
   # 
   # code for 'wild-raised' column
   #
   annotate(
     "text", x = 0.85, y = phi_ad + 0.1, label = "Wild-raised"
   ) +
   # true adult survival
   geom_line(
     data = tibble(x = c(0.7, 1), y = rep(phi_ad, 2)),
     mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
   ) +  
   # true immature survival
   geom_line(
     data = tibble(x = c(0.7, 1), y = rep(phi_im, 2)),
     mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
   ) +
   # estimate phi_im = phi_ad
   geom_line(
     data = tibble(x = c(0.5, 1), y = rep(0.2*phi_im + 0.8*phi_ad, 2)),
     mapping = aes(x = x, y = y), linetype = "dashed", colour = "darkred"
   ) +
   annotate(
     "text", x = 1.1, y = 0.2*phi_im + 0.8*phi_ad, size = 4,
     label = expression(hat(varphi)[imm]^"wr"~"="~hat(varphi)[ad]^"wr"),
     colour = "darkred"
   ) +
   # arrow for bias in phi_im
   # annotate(
   #   "segment", x = 0.85, y = phi_im, xend = 0.85, yend = 0.2*phi_im + 0.8*phi_ad - 0.01, 
   #   size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
   #   arrow = arrow(type = "closed", length = unit(0.015, "npc"))
   # ) +
   #
   # arrows for differential bias
   #
   annotate(
     "segment", x = 0.6, xend = 0.6,
     y = 0.2*phi_im + 0.8*phi_ad, yend = 0.6*phi_im + 0.4*phi_ad + 0.01,
     size = 1.5, linejoin = "mitre", colour = "darkorange", alpha = 0.8,
     arrow = arrow(type = "closed", length = unit(0.015, "npc"))
   ))
# ggsave("figs/21_phi-im_bias_plot_d.svg", plot_d, device = svg)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Differential bias in adults and juveniles ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(plot_c <- ggplot() +
  geom_point(
    data = tibble(x = runif(100, 0, 1.2), y = runif(100, 0, 1)), 
    mapping = aes(x = x, y = y),
    colour = "white"
  ) + 
  #
  # global formatting
  #
  theme_classic() +
  theme(
    # aspect.ratio = 0.7/1.2,
    aspect.ratio = 0.5,
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  coord_cartesian(ylim = c(0.2, 0.9)) +
  labs(x = "Hand-rearing status", y = "Survival probability") +
  # 
  # code for 'hand-reared' column
  #
  annotate(
    "text", x = 0.35, y = phi_ad + 0.1, label = "Hand-reared"
  ) +
  # true adult survival
  geom_line(
    data = tibble(x = c(0.2, 0.5), y = rep(phi_ad, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
  ) +  
  annotate(
    "text", x = 0.125, y = phi_ad, size = 4,
    label = expression(varphi[ad]^"truth")
  ) +
  # true immature survival
  geom_line(
    data = tibble(x = c(0.2, 0.5), y = rep(phi_im, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
  ) +
  annotate(
    "text", x = 0.125, y = phi_im, size = 4,
    label = expression(varphi[imm]^"truth")
  ) +
  # true juvenile survival
  geom_line(
    data = tibble(x = c(0.2, 0.5), y = rep(phi_jv, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
  ) +
  annotate(
    "text", x = 0.125, y = phi_jv, size = 4,
    label = expression(varphi[juv]^"truth")
  ) +
  # estimate phi_im = phi_ad
  geom_line(
    data = tibble(x = c(0.2, 0.7), y = rep(0.6*phi_im + 0.4*phi_ad, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", colour = "darkred"
  ) +
  annotate(
    "text", x = 0.075, y = 0.6*phi_im + 0.4*phi_ad + 0.02, size = 4,
    label = expression(hat(varphi)[imm]^"hr"~"="~hat(varphi)[ad]^"hr"),
    colour = "darkred"
  ) +
  # estimate phi_jv
  geom_line(
    data = tibble(x = c(0.2, 0.7), y = rep(phi_jv - 0.4*(phi_ad - phi_im), 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", colour = "darkred"
  ) +
  annotate(
    "text", x = 0.125, y = phi_jv - 0.4*(phi_ad - phi_im), size = 4,
    label = expression(hat(varphi)[juv]^"hr"),
    colour = "darkred"
  ) +
  # arrow for bias in phi_im
  annotate(
    "segment", x = 0.35, y = phi_im, xend = 0.35, yend = 0.6*phi_im + 0.4*phi_ad - 0.01,
    size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  # arrow for bias in phi_jv
  annotate(
    "segment", x = 0.35, y = phi_jv, xend = 0.35, yend = phi_jv - 0.4*(phi_ad - phi_im) + 0.01,
    size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  # 
  # code for 'wild-raised' column
  #
  annotate(
    "text", x = 0.85, y = phi_ad + 0.1, label = "Wild-raised"
  ) +
  # true adult survival
  geom_line(
    data = tibble(x = c(0.7, 1), y = rep(phi_ad, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
  ) +  
  # true immature survival
  geom_line(
    data = tibble(x = c(0.7, 1), y = rep(phi_im, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
  ) +
  # true juvenile survival
  geom_line(
    data = tibble(x = c(0.7, 1), y = rep(phi_jv, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", alpha = 0.6
  ) +
  # estimate phi_im = phi_ad
  geom_line(
    data = tibble(x = c(0.5, 1), y = rep(0.2*phi_im + 0.8*phi_ad, 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", colour = "darkred"
  ) +
  annotate(
    "text", x = 1.125, y = 0.2*phi_im + 0.8*phi_ad, size = 4,
    label = expression(hat(varphi)[imm]^"wr"~"="~hat(varphi)[ad]^"wr"),
    colour = "darkred"
  ) +
  # estimate phi_jv
  geom_line(
    data = tibble(x = c(0.5, 1), y = rep(phi_jv - 0.8*(phi_ad - phi_im), 2)),
    mapping = aes(x = x, y = y), linetype = "dashed", colour = "darkred"
  ) +
  annotate(
    "text", x = 1.05, y = phi_jv - 0.8*(phi_ad - phi_im), size = 4,
    label = expression(hat(varphi)[juv]^"wr"),
    colour = "darkred"
  ) +
  # arrow for bias in phi_im
  annotate(
    "segment", x = 0.85, y = phi_im, xend = 0.85, yend = 0.2*phi_im + 0.8*phi_ad - 0.01, 
    size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  # arrow for bias in phi_jv
  annotate(
    "segment", x = 0.85, y = phi_jv, xend = 0.85, yend = phi_jv - 0.8*(phi_ad - phi_im) + 0.01, 
    size = 1, linejoin = "mitre", colour = "darkred", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  #
  # arrows for differential bias
  #
  annotate(
    "segment", x = 0.6, xend = 0.6,
    y = phi_jv - 0.8*(phi_ad - phi_im), yend = phi_jv - 0.4*(phi_ad - phi_im) - 0.01, 
    size = 1.5, linejoin = "mitre", colour = "navyblue", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ) +
  annotate(
    "segment", x = 0.6, xend = 0.6,
    y = 0.2*phi_im + 0.8*phi_ad, yend = 0.6*phi_im + 0.4*phi_ad + 0.01,
    size = 1.5, linejoin = "mitre", colour = "darkorange", alpha = 0.8,
    arrow = arrow(type = "closed", length = unit(0.015, "npc"))
  ))
# ggsave("figs/21_phi-im_bias_plot_c.svg", plot_c, device = svg)

