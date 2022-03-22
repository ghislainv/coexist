data {
  int<lower=0> N;              // num individuals (observations)
  int<lower=1> K;              // num ind predictors (= nb of species)
  int<lower=1> J;              // num groups (= nb of species)
  int<lower=1> L;              // num group predictors (1 intercept)
  int<lower=1,upper=J> jj[N];  // group for individual (species)
  matrix[N, K] x;              // individual predictors (all abundances at time t)
  matrix[L, J] u;          // group predictors (1 everywhere)
  vector[N] y;                 // outcomes (abundance of focal species at time t)
}
parameters {
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;  // prior scale
  matrix[K, L] gamma;                        // group coeffs
  real<lower=0> sigma;                       // prediction error scale
}
transformed parameters {
  vector<lower=0>[K] tau = 2.5 * tan(tau_unif);
  matrix[K, J] beta = gamma * u + diag_pre_multiply(tau, L_Omega) * z;
}
model {
  vector[N] mu;
  for(n in 1:N) {
    mu[n] = x[n, ] * beta[, jj[n]];
  }
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 5);
  y ~ normal(mu, sigma);
}
