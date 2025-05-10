# Engineering a 'spike-and-slab' prior for probability parameters as a mixture
# of beta distributions. We dub the result a 'beta-spike-and-slab' prior, or
# 'bss' prior for short. 
library(tidyverse)

# Compare histograms and analytically derived density
alpha_for_spike <- 99
beta_for_spike <- 1
alpha_for_slab <- 2
beta_for_slab <- 1
slab_weight <- 5 # 'alpha' for mixture distribution
spike_weight <- 25 # 'beta' for mixture distribution

# Large sample
N <- 1e4
y <- rep(NA, N)
for (i in 1:N) {
  rho <- rbeta(1,spike_weight, slab_weight)
  zeta <- rbinom(1,1,rho)
  if (zeta == 1) {
    y[i] <- rbeta(1, alpha_for_spike, beta_for_spike)
  }
  if (zeta == 0) {
    y[i] <- rbeta(1, alpha_for_slab, beta_for_slab)
  }
}

hist(y, probability = TRUE, xlim = c(0,1), breaks = seq(0,1,length.out = 101))
curve(
  (spike_weight / (slab_weight + spike_weight)) * dbeta(x, alpha_for_spike, beta_for_spike) +
    (slab_weight / (slab_weight + spike_weight)) * dbeta(x, alpha_for_slab, beta_for_slab),
  0, 1, add = TRUE, lwd = 2
)

# Moderate sample
N <- 30
y <- rep(NA, N)
for (i in 1:N) {
  rho <- rbeta(1,spike_weight, slab_weight)
  zeta <- rbinom(1,1,rho)
  if (zeta == 1) {
    y[i] <- rbeta(1, alpha_for_spike, beta_for_spike)
  }
  if (zeta == 0) {
    y[i] <- rbeta(1, alpha_for_slab, beta_for_slab)
  }
}

hist(y, xlim = c(0,1), breaks = seq(0,1,length.out = 41), ylim = c(0,30))



