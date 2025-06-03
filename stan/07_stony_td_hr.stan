//
// Multi-age, single-site model incorporating hand-rearing effect on adult and 
// juvenile survival. Trap-dependence is not modelle and immature survival is
// not separated from adult survival
// 
// We use informative priors on movement.
// 
// The model's states are defined to be:
// 1 - Stony, 0yrs old, 
// 2 - Stony, 1yrs old, 
// 3 - Stony, 2+ yrs old, trap-aware
// 4 - Stony, 2+ yrs old, trap-unaware
//

functions { 
  matrix get_multinomial_probs(
    data int T,
    data real eps,
    vector phi_jv,
    vector phi_ad,
    real m_jv,
    real m_ad,
    vector p_im,
    vector p_ad_N,
    vector p_ad_A,
    vector p_ad_U
    ) {
    // matrix to store multinomial probabilities
    matrix[4*(T-1), 4*(T-1)+1] pr = rep_matrix(0.0, 4*(T-1), 4*(T-1)+1);  

    // auxillary matrices to define multinomial probabilities
    array[T-1] matrix[4,4] gamma;  // transition matrices
    array[T] matrix[4,2] omega;     // observation matrices
    // define gamma
    for (t in 1:(T-1)) {
      gamma[t] = [
        [eps, phi_jv[t]*m_jv,                      eps,                            eps],
        [eps,            eps, phi_ad[t]*m_jv*p_ad_N[t], phi_ad[t]*m_jv*(1 - p_ad_N[t])],
        [eps,            eps, phi_ad[t]*m_ad*p_ad_A[t], phi_ad[t]*m_ad*(1 - p_ad_A[t])],
        [eps,            eps, phi_ad[t]*m_ad*p_ad_U[t], phi_ad[t]*m_ad*(1 - p_ad_U[t])]
      ];
    }
    // define omega
    for (t in 1:T) {
      omega[t] = [
        [1 - eps,       eps],  // juv
        [p_im[t], 1-p_im[t]],  // immature
        [1 - eps,       eps],  // trap-aware adult
        [    eps,     1-eps]   // trap-unaware adult
      ];
    }
    // define entries of pr using matrix multiplication
    for (I in 1:(T-1)) {  // block row index
      for (J in I:(T-1)) {  // block column index
        // row and column indices to define the 12x12 block
        int i1 = 1 + 4*(I-1);
        int i2 = 4*I;
        int j1 = 1 + 4*(J-1);
        int j2 = 4*J;
        if (I == J) {  // diagonal block
          pr[i1:i2, j1:j2] = diag_post_multiply(gamma[I], omega[I+1][:,1]);
        }
        else if (I < J) {  // above-diagonal block
          matrix[4, 4] temp = diag_post_multiply(gamma[I], omega[I+1][:,2]);
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
    for (i in 1:(4*(T-1))){
      pr[i, 4*(T-1)+1] = 1 - sum(pr[i, 1:4*(T-1)]);
    }
      return pr;
    }
}

data {
  int<lower=1> T;                                   // number of years
  array[4*(T-1), 4*(T-1)+1] int<lower=0> marr_wr;   // m-array wild-raised
  array[4*(T-1), 4*(T-1)+1] int<lower=0> marr_hr;   // m-array hand-raised
  // prior parameters for probability of not moving
  real<lower=0> m_jv_alpha;
  real<lower=0> m_jv_beta;
  real<lower=0> m_ad_alpha;
  real<lower=0> m_ad_beta;
}

transformed data {
  real<lower=0> eps = pow(10, -15);
}

parameters {
  // time-varying survival probabilities  
  vector<lower=0, upper=1>[T-1] phi_jv_wr;   // for juveniles (0 yrs old)
  vector<lower=0, upper=1>[T-1] phi_ad_wr;   // for immatures, adults (1, 2+ yrs old)
  // logit-scale additive effect for hand-rearing, shared across years
  real hr_jv;
  real hr_ad;
  // time-varying detection for immatures and adults
  vector<lower=0, upper=1>[T] p_im;
  vector<lower=0, upper=1>[T] p_ad_N;
  vector<lower=0, upper=1>[T] p_ad_A;
  vector<lower=0, upper=1>[T] p_ad_U;
  // probabilities of not moving to Robben or Boulders
  real<lower=0,upper=1> m_jv;  // juvs, imms (0, 1 yrs old)
  real<lower=0,upper=1> m_ad;  // adults (2+ yrs old)
}

transformed parameters {
  vector<lower=0, upper=1>[T-1] phi_jv_hr; 
  vector<lower=0, upper=1>[T-1] phi_ad_hr; 
  phi_jv_hr = inv_logit( logit( phi_jv_wr ) + hr_jv );
  phi_ad_hr = inv_logit( logit( phi_ad_wr ) + hr_ad );
}

model {
  // default uniform priors on all probabilities
  // probabilities of not moving from Stony given informative prior based on
  // model 
  m_jv ~ beta(m_jv_alpha, m_jv_beta);
  m_ad ~ beta(m_ad_alpha, m_ad_beta);
  
  // calculate multinomial probabilities
  matrix[4*(T-1), 4*(T-1)+1] pr_wr;
  matrix[4*(T-1), 4*(T-1)+1] pr_hr;
  pr_wr = get_multinomial_probs( 
    T, eps, phi_jv_wr, phi_ad_wr, m_jv, m_ad, p_im, p_ad_N, p_ad_A, p_ad_U
  ); 
  pr_hr = get_multinomial_probs( 
    T, eps, phi_jv_hr, phi_ad_hr, m_jv, m_ad, p_im, p_ad_N, p_ad_A, p_ad_U
  ); 
  // m-array capture-recapture likelihood
  for (i in 1:(4*(T-1))) {
    marr_wr[i] ~ multinomial(pr_wr[i]');
    marr_hr[i] ~ multinomial(pr_hr[i]');
  }
}
