//
// Multinomial mixture model
// Given two sets of multinomial data, the model specifies each multinomial in
// one set as a mixture of those in the other set.
// see Choquet et al U-CARE User Guide, section 5.5 for details
//

data {
  int<lower=2> M;           // number of cells in all multinomial distributions
  int<lower=2> N_comp;      // number of component multinomials
  int<lower=1> N_mix;       // number of multinomials modelled as mixtures
  array[N_comp, M] int<lower=0> Y_comp;  // data for component multinomials
  array[N_mix, M] int<lower=0> Y_mix;    // data for mixture multinomials
}

parameters {
  array[N_comp] simplex[M] p_comp;       // cell probabilities of component multinomials
  array[N_mix] simplex[N_comp] mix_wts;   // mixture weights
}

transformed parameters {
  // arrange cell probabilities in the rows of a matrix
  matrix[N_comp, M] p_comp_mat;
  for (i in 1:N_comp) {
    p_comp_mat[i] = p_comp[i]';
  }
  
  // define cell probabilities for mixtures
  array[N_mix] vector[M] p_mix;
  for (i in 1:N_mix) {
    p_mix[i] = (mix_wts[i]' * p_comp_mat)';
  }
}

model {
  // likelihood
  for (i in 1:N_comp) {
    Y_comp[i] ~ multinomial(p_comp[i]);
  }
  for (i in 1:N_mix) {
    Y_mix[i] ~ multinomial(p_mix[i]);
  }
}
