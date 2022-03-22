data {
  int<lower=0> N;              // num individuals (observations)
  int<lower=1> K;              // num ind predictors (= nb of species)
  int<lower=1> J;              // num groups (= nb of species)
  int<lower=1,upper=J> jj[N];  // group for individual (species)
  matrix[N, K] x;              // individual predictors (all abundances at time t)
  vector[N] y;                 // outcomes (abundance of focal species at time t)
}
parameters {
  vector[K] beta[J];           // indiv coeffs by group
  real<lower=0> sigma;         // prediction error scale
}
model {
  {
    vector[N] x_beta_jj;
    for (n in 1:N)
      x_beta_jj[n] = x[n] * beta[jj[n]];
    y ~ normal(x_beta_jj, sigma);
  }
}
