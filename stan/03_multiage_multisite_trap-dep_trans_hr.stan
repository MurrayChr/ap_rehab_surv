//
// Multi-age, multi-site model incorporating trap-dependence and transience in the
// adult age class. The treatment of trap-dependence follows Pradel and Sanz-Aguilar (2012). 
//
// Here we incorporate hand-rearing effects on adult and juvenile survival.
// 
// The model's states are defined to be:
// 1 - Robben, 0yrs old, 
// 2 - Robben, 1yrs old, 
// 3 - Robben, 2+ yrs old, trap-aware, 
// 4 - Robben, 2+ yrs old, trap-unaware, 
// 5 - Boulders, 0yrs old, 
// 6 - Boulders, 1yrs old, 
// 7 - Boulders, 2+ yrs old, trap-aware, 
// 8 - Boulders, 2+ yrs old, trap-unaware, 
// 9 - Stony, 0yrs old, 
// 10 - Stony, 1yrs old, 
// 11 - Stony, 2+ yrs old, trap-aware, 
// 12 - Stony, 2+ yrs old, trap-unaware, 
//
// Important! the m-array for wild-raised birds must be constructed after first *excluding* all
// birds marked as adults that were never recaptured after first capture.
// (This is not an issue for hand-reared birds, because they are all marked as juveniles 
//  and our residence/transience probabilities only refer to newly marked adults).
//

functions { 
  matrix get_multinomial_probs(
    data int T,
    data real eps,
    array[] vector phi_jv,
    array[] vector phi_ad,
    array[] vector m_jv,
    array[] vector m_ad,
    array[] vector p_im,
    array[] vector p_ad_N,
    array[] vector p_ad_A,
    array[] vector p_ad_U
    ) {
    // matrix to store multinomial probabilities
    matrix[12*(T-1), 12*(T-1)+1] pr = rep_matrix(0.0, 12*(T-1), 12*(T-1)+1);  

    // auxillary matrices to define multinomial probabilities
    array[T-1] matrix[12,12] gamma;  // transition matrices
    array[T] matrix[12,2] omega;     // observation matrices
    // define gamma
    for (t in 1:(T-1)) {
      // define 4x4 site blocks
      for (I in 1:3) {     // I indexes source site
        int first_i = 4*(I-1) + 1;
        int last_i = 4*I;
        for (J in 1:3) {   // J indexes destination site
          int first_j = 4*(J-1) + 1;
          int last_j = 4*J;
          if (I==J) {  // remain at same site
            gamma[t][first_i:last_i, first_j:last_j] = [
              [ eps, phi_jv[I][t]*m_jv[I][J],                                    eps,                                        eps ],
              [ eps,                     eps, phi_ad[I][t]*m_jv[I][J]*p_ad_N[J][t+1], phi_ad[I][t]*m_jv[I][J]*(1-p_ad_N[J][t+1]) ],
              [ eps,                     eps, phi_ad[I][t]*m_ad[I][J]*p_ad_A[J][t+1], phi_ad[I][t]*m_ad[I][J]*(1-p_ad_A[J][t+1]) ],
              [ eps,                     eps, phi_ad[I][t]*m_ad[I][J]*p_ad_U[J][t+1], phi_ad[I][t]*m_ad[I][J]*(1-p_ad_U[J][t+1]) ]
            ];
          } 
          if (I!=J) {  // move to different site: only permitted if detected i.e. no transition into 'trap-unaware' state
            gamma[t][first_i:last_i, first_j:last_j] = [
              [ eps, phi_jv[I][t]*m_jv[I][J],                                    eps, eps ],
              [ eps,                     eps, phi_ad[I][t]*m_jv[I][J]*p_ad_N[J][t+1], eps ],
              [ eps,                     eps, phi_ad[I][t]*m_ad[I][J]*p_ad_A[J][t+1], eps ],
              [ eps,                     eps, phi_ad[I][t]*m_ad[I][J]*p_ad_A[J][t+1], eps ]    // does p_ad_A make sense here?
            ];
          }
        }
      }
    }
    // define omega
    for (t in 1:T) {
      for (I in 1:3) {
        int first_i = 4*(I-1) + 1;
        int last_i = 4*I;
        // define 4x2 site blocks
        omega[t][first_i:last_i, 1:2] = [
          [    1 - eps,          eps ],  // juv
          [ p_im[I][t], 1-p_im[I][t] ],  // immature
          [    1 - eps,          eps ],  // adult, trap-aware
          [        eps,      1 - eps ]  // adult, trap-unaware
        ];
      }
    }
    // define entries of pr using matrix multiplication
    for (I in 1:(T-1)) {  // block row index
      for (J in I:(T-1)) {  // block column index
        // row and column indices to define the 12x12 block
        int i1 = 1 + 12*(I-1);
        int i2 = 12*I;
        int j1 = 1 + 12*(J-1);
        int j2 = 12*J;
        if (I == J) {  // diagonal block
          pr[i1:i2, j1:j2] = diag_post_multiply(gamma[I], omega[I+1][:,1]);
        }
        else if (I < J) {  // above-diagonal block
          matrix[12,12] temp = diag_post_multiply(gamma[I], omega[I+1][:,2]);
          if ( (I+1) < J ) {
            for (t in (I+1):(J-1)) {
              temp = temp * diag_post_multiply(gamma[t], omega[t+1][:,2]);
            }
          }
          pr[i1:i2, j1:j2] = temp * diag_post_multiply(gamma[J], omega[J+1][:,1]);
        }
      }
    }
    // enforce row-sum-to-one constraint
    for (i in 1:(12*(T-1))){
      pr[i, 12*(T-1)+1] = 1 - sum(pr[i, 1:12*(T-1)]);
    }
      return pr;
    }
}

data {
  int<lower=1> T;                                  // number of years
  array[12*(T-1), 12*(T-1)+1] int<lower=0> marr_hr;  // multi-state m-array for hand-reared birds
  array[12*(T-1), 12*(T-1)+1] int<lower=0> marr_wr;  // multi-state m-array for wild-raised birds (see above)
  // number of birds marked as adults (thus wild-raised) per site and cohort that were...
  array[3, T-1] int<lower=0> N_1;                  // ... ever recaptured
  array[3, T-1] int<lower=0> N_0;                  // ... never recaptured
}

transformed data {
  real<lower=0> eps = pow(10, -15);
}

parameters {
  // time- and site-varying survival probabilities for wild-raised birds  
  array[3] vector<lower=0, upper=1>[T-1] phi_jv_wr;   // for juveniles (0 yrs old)
  array[3] vector<lower=0, upper=1>[T-1] phi_ad_wr;   // for immatures, adults (1, 2+ yrs old)
  // logit-scale additive effect for hand-rearing, shared across sites and years
  real hr_jv;
  real hr_ad;
  // proportion 'residents' among birds marked as adults per cohort and site
  array[3] vector<lower=0,upper=1>[T-1] pi_r;    
  // time- and site-varying detection for immatures
  array[3] vector<lower=0, upper=1>[T] p_im;
  // time-and site-varying, trap-dependent detection probabilities for adults
  array[3] vector<lower=0, upper=1>[T] p_ad_A;  // 'trap-aware'
  array[3] vector<lower=0, upper=1>[T] p_ad_U;  // 'trap-unaware'
  // time- and site-varying prob. of det. for newly matured adults (exactly 2 yrs old) 
  array[3] vector<lower=0, upper=1>[T] p_ad_N;
  // time-constant, age-dependent movement probabilities
  array[3] simplex[3] m_jv;  // juvs, imms (0, 1 yrs old)
  array[3] simplex[3] m_ad;  // adults (2+ yrs old)
}

transformed parameters {
  array[3] vector<lower=0, upper=1>[T-1] phi_jv_hr; 
  array[3] vector<lower=0, upper=1>[T-1] phi_ad_hr; 
  for (s in 1:3) {
    phi_jv_hr[s] = inv_logit( logit( phi_jv_wr[s] ) + hr_jv );
    phi_ad_hr[s] = inv_logit( logit( phi_ad_wr[s] ) + hr_ad );
  }
}

model {
  // default uniform priors on all probabilities
  // vague priors on hand-rearing effects (on logit scale)
  hr_jv ~ normal(0, 1);
  hr_ad ~ normal(0, 1);
  // calculate multinomial probabilities
  matrix[12*(T-1), 12*(T-1)+1] pr_wr;
  matrix[12*(T-1), 12*(T-1)+1] pr_hr;
  pr_wr = get_multinomial_probs( 
    T, eps, phi_jv_wr, phi_ad_wr, m_jv, m_ad, p_im, p_ad_N, p_ad_A, p_ad_U
  ); 
  pr_hr = get_multinomial_probs( 
    T, eps, phi_jv_hr, phi_ad_hr, m_jv, m_ad, p_im, p_ad_N, p_ad_A, p_ad_U
  ); 
  // likelihood
  // m-array for all birds except those marked as adults and never seen again
  for (i in 1:(12*(T-1))) {
    marr_wr[i] ~ multinomial(pr_wr[i]');
    marr_hr[i] ~ multinomial(pr_hr[i]');
  }
  for (i in 1:3) {
    for (t in 1:(T-1)) {
      target += N_1[i,t] * log(pi_r[i][t]);
    }
  }
  //  for birds marked as adults and never seen again
  for (i in 1:3) {
    for (t in 1:(T-1)) {
      int ind = 3 + 4*(i-1) + 12*(t-1); 
      target += N_0[i,t] * log( pr_wr[ind][12*(T-1)+1] * pi_r[i][t] + (1 - pi_r[i][t]) );
    }
  }
}

