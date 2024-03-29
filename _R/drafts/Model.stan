data {
  int<lower=0> N;              // num individuals (observations)
  int<lower=1> K;              // num ind predictors (= nb of species)
  int<lower=1> J;              // num groups (= nb of species)
  int<lower=1> L;              // num group predictors (1?)
  int<lower=1,upper=J> jj[N];  // group for individual (species)
  matrix[N, K] x;              // individual predictors (all abundances at time t)
  row_vector[L] u[J];          // group predictors (1 everywhere)
  vector[N] y;                 // outcomes (abundance of focal species at time t)
}
parameters {
  corr_matrix[K] Omega;        // prior correlation
  vector<lower=0>[K] tau;      // prior scale
  matrix[L, K] gamma;          // group coeffs
  vector[K] beta[J];           // indiv coeffs by group
  real<lower=0> sigma;         // prediction error scale
}
model {
  tau ~ cauchy(0, 2.5);
  Omega ~ lkj_corr(1);
  to_vector(gamma) ~ normal(0, 5);
  {
    row_vector[K] u_gamma[J];
    for (j in 1:J)
      u_gamma[j] = u[j] * gamma;
    beta ~ multi_normal(u_gamma, quad_form_diag(Omega, tau));
  }
  {
    vector[N] x_beta_jj;
    for (n in 1:N)
      x_beta_jj[n] = x[n] * beta[jj[n]];
    y ~ normal(x_beta_jj, sigma);
  }
}
