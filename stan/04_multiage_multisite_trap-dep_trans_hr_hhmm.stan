//
// Multi-age, multi-site model incorporating trap-dependence and transience in the
// adult age class, and hand-rearing effects on adult and juvenile survival. 
// The treatment of trap-dependence follows Pradel and Sanz-Aguilar (2012). 
// 
// Here we use a 'multievent' or hierarchical hidden Markov model likelihood. 
//
// The model's states are defined to be:
// 1 - Robben, 0yrs old, 
// 2 - Robben, 1yrs old, 
// 3 - Robben, 2+ yrs old, transient
// 4 - Robben, 2+ yrs old, resident, trap-aware, 
// 5 - Robben, 2+ yrs old, resident, trap-unaware, 
// 6 - Boulders, 0yrs old, 
// 7 - Boulders, 1yrs old, 
// 8 - Boulders, 2+ yrs old, transient
// 9 - Boulders, 2+ yrs old, resident, trap-aware, 
// 10 - Boulders, 2+ yrs old, resident, trap-unaware, 
// 11 - Stony, 0yrs old, 
// 12 - Stony, 1yrs old, 
// 13 - Stony, 2+ yrs old, transient 
// 14 - Stony, 2+ yrs old, resident, trap-aware, 
// 15 - Stony, 2+ yrs old, resident, trap-unaware, 
// 16 - dead
//

// The capture histories are coded as:
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
// 13 - not captured
//
// Important note: in calculating unique capture histories and their
// multiplicities, individuals with the same history, but different hand-rearing 
// covariate values must be treated as distinct because their likelihoods differ
//

functions {
  // create state transition matrices
  array[] matrix get_Gamma_matrices(
    data int T,
    array[] vector phi_jv,
    array[] vector phi_ad,
    array[] vector pi_r,
    array[] vector m_jv,
    array[] vector m_ad,
    array[] vector p_ad_A,
    array[] vector p_ad_U,
    array[] vector p_ad_N
  ) {
    // convenience detection parameters
    array[3] vector[T] q_ad_A;
    array[3] vector[T] q_ad_U;
    array[3] vector[T] q_ad_N;
    for (s in 1:3){
      q_ad_A[s] = 1 - p_ad_A[s];
      q_ad_U[s] = 1 - p_ad_U[s];
      q_ad_N[s] = 1 - p_ad_N[s];
    }
    // construct transition matrices 
    array[T-1] matrix[16, 16] Gamma;
    for (t in 1:(T-1)) {
      // define 5x5 site blocks
      for (I in 1:3) {         // I indexes source site
        int first_i = 5*(I-1) + 1;
        int last_i = 5*I;
        for (J in 1:3) {       // J indexes destination site
          int first_j = 5*(J-1) + 1;
          int last_j = 5*J;
          Gamma[t][first_i:last_i, first_j:last_j] = [
            [0, phi_jv[I][t]*m_jv[I][J], 0,                                      0,                                      0],
            [0,                       0, 0, phi_ad[I][t]*m_jv[I][J]*p_ad_N[J][t+1], phi_ad[I][t]*m_jv[I][J]*q_ad_N[J][t+1]],
            [0,                       0, 0,                                      0,                                      0],
            [0,                       0, 0, phi_ad[I][t]*m_ad[I][J]*p_ad_A[J][t+1], phi_ad[I][t]*m_ad[I][J]*q_ad_A[J][t+1]],
            [0,                       0, 0, phi_ad[I][t]*m_ad[I][J]*p_ad_U[J][t+1], phi_ad[I][t]*m_ad[I][J]*q_ad_U[J][t+1]]
          ];
        }
      }
      // last row (except last entry)
      Gamma[t][16, 1:15] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
      // last-column (rows sum to one) 
      Gamma[t][1:16,16] = [1-phi_jv[1][t], 1-phi_ad[1][t], 1, 1-phi_ad[1][t], 1-phi_ad[1][t],
                           1-phi_jv[2][t], 1-phi_ad[2][t], 1, 1-phi_ad[2][t], 1-phi_ad[2][t],
                           1-phi_jv[3][t], 1-phi_ad[3][t], 1, 1-phi_ad[3][t], 1-phi_ad[3][t],
                           1]';
      
    }
    return Gamma;
  }
  
  // create observation matrices
  array[] matrix get_Omega_matrices(
    data int T,
    array[] vector p_im
  ) {
    array[T] matrix[16, 13] Omega; // there are 13 codes in capture histories
    for (t in 1:T) {
      // define in 5x4 site blocks
      for (I in 1:3) {
        int first_i = 5*(I-1) + 1;
        int last_i = 5*I;
        for (J in 1:3) {
          int first_j = 4*(J-1) + 1;
          int last_j = 4*J;
          if (I==J) {  // non-zero block
            Omega[t][first_i:last_i, first_j:last_j] = [
              [1,          0, 0, 0],
              [0, p_im[I][t], 0, 0],
              [0,          0, 0, 0],
              [0,          0, 1, 0],
              [0,          0, 0, 0]
            ];
          }
          if (I!=J) {  // zero block
            Omega[t][first_i:last_i, first_j:last_j] = rep_matrix(0, 5, 4);
          }
        }
      }
      // last row (except last entry)
      Omega[t][16, 1:12] = rep_row_vector(0, 12);
      // last column
      Omega[t][:,13] = [0, 1-p_im[1][t], 1, 0, 1,
                        0, 1-p_im[2][t], 1, 0, 1,
                        0, 1-p_im[3][t], 1, 0, 1,
                        1]'; 
    }
    return Omega;
  }
}

data {
  int<lower=1> T;                       // number of years
  int<lower=1> N;                       // no. unique capture histories (see important note above)
  array[N,T] int<lower=0, upper=13> y;  // unique captures histories (see important note above)
  array[N] int<lower=1,upper=T-1> fc;   // first capture occasion
  array[N] int fc_state;                // first capture state
  array[N] int mult;                    // capture history multiplicity
  array[N] int hr;                      // hand-rearing covariate
}

transformed data {
  // check input data
  for (n in 1:N) {
    if (fc_state[n]%4 == 0) {
      reject("Cannot be trap-unaware at first capture (codes 4, 8, 12).");
    }
  }
  // create quotient and remainder for fc_state (once upfront, not in model block)
  // the capture history codes 1,...,12 repeat in length 4 'blocks' per site
  // so writing f=fc_state[n] as f = 4q + r for integers 0 <= q <= 2, 0 <= r <= 3
  // gives the site (equal to q+1) and the state within site:
  // r = 1 for juveniles
  // r = 2 for immatures
  // r = 3 for trap-aware adults
  // r = 0 for trap-unaware adults
  array[N] int<lower=0, upper=2> fc_state_q;
  array[N] int<lower=0, upper=3> fc_state_r;
  for (n in 1:N) {
    fc_state_q[n] = fc_state[n] %/% 4;
    fc_state_r[n] = fc_state[n] % 4;
  }
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
  // get transition and emission matrices
  array[T-1] matrix[16,16] Gamma_wr = get_Gamma_matrices(
    T, phi_jv_wr, phi_ad_wr, pi_r, m_jv, m_ad, p_ad_A, p_ad_U, p_ad_N );
  array[T-1] matrix[16,16] Gamma_hr = get_Gamma_matrices(
    T, phi_jv_hr, phi_ad_hr, pi_r, m_jv, m_ad, p_ad_A, p_ad_U, p_ad_N );
  array[T] matrix[16,13] Omega = get_Omega_matrices( T, p_im );
  
  // default uniform priors on all probabilities
  // informative Dirichlet priors on movement simplices
  m_jv[1] ~ dirichlet([7, 1.5, 1.5]);
  m_jv[2] ~ dirichlet([1.5, 7, 1.5]);
  m_jv[3] ~ dirichlet([1.5, 1.5, 7]);
  m_ad[1] ~ dirichlet([8, 1, 1]);
  m_ad[2] ~ dirichlet([1, 8, 1]);
  m_ad[3] ~ dirichlet([1, 1, 8]);
  // vague priors on hand-rearing effects (on logit scale)
  hr_jv ~ normal(0, 1);
  hr_ad ~ normal(0, 1);
  
  // likelihood
  for (n in 1:N) {
    // initialise forward algorithm
    array[T] row_vector[16] fwd;
    fwd[fc[n]] = rep_row_vector(0, 16);
    if (fc_state_r[n] == 1) { // juvenile
      fwd[fc[n]][5*fc_state_q[n] + fc_state_r[n]] = 1;
    }
    if (fc_state_r[n] == 2) {  // immature  
      fwd[fc[n]][5*fc_state_q[n] + fc_state_r[n]] = 1;
    }
    if (fc_state_r[n] == 3) { // adult ('trap-aware')
      real res_prob = pi_r[fc_state_q[n]+1][fc[n]];
      fwd[fc[n]][5*fc_state_q[n]+3] = 1-res_prob;  // transient adult
      fwd[fc[n]][5*fc_state_q[n]+5] = res_prob;    // resident adult
    }
    // recursion
    if (hr[n] == 1) {  // hand-reared
      for (t in fc[n]:(T-1)) {
        fwd[t+1] = fwd[t]*diag_post_multiply(Gamma_hr[t], Omega[t+1][:,y[n,t+1]]);
      }
    }
    if (hr[n] == 2) { // wild-raised
      for (t in fc[n]:(T-1)) {
        fwd[t+1] = fwd[t]*diag_post_multiply(Gamma_wr[t], Omega[t+1][:,y[n,t+1]]);
      }
    }
    // increment log-likelihood 
    target += mult[n]*log(sum(fwd[T]));
  }
}
