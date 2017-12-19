data {
  int<lower = 1> num_cov;             //number covariate categories
  int<lower = num_cov> N;             //number of individuals
  vector[N] phenotype;                //phenotype
  int cov[N];                         //covariate
  corr_matrix[N] R;                   //kernel
}

transformed data {
  matrix[N, N] L;
  L = cholesky_decompose(R);
}

parameters {
  real grand_mean;
  real<lower = 0.01> grand_sig2;
  real<lower = 0.01> sig2_mef;
  real<lower = 0.01> sig2_vef;
  real<lower = 0.01> sig2_genef; 
  
  vector [num_cov] cov_mef;
  vector [num_cov] cov_vef;
  vector[N] u;
}

transformed parameters {
  vector[num_cov] h2;
  real<lower = 0.01> sum_sig2_1;
  real<lower = 0.01> sum_sig2_2;
  sum_sig2_1 = (grand_sig2 * exp(cov_vef[1])) + sig2_genef;
  sum_sig2_2 = (grand_sig2 * exp(cov_vef[2])) + sig2_genef;
  h2[1] = sig2_genef/sum_sig2_1;
  h2[2] = sig2_genef/sum_sig2_2;
}

model {
  // priors
  grand_mean   ~ normal(0, 100);
  grand_sig2   ~ inv_gamma(1, 1);
  sig2_mef     ~ inv_gamma(1, 1);
  sig2_vef     ~ inv_gamma(1, 1);
  sig2_genef   ~ inv_gamma(1, 1);

  // sex and genetic effects
  cov_mef ~ normal(0, sqrt(sig2_mef));
  cov_vef ~ normal(0, sqrt(sig2_vef));
  u       ~ normal(0, sqrt(sig2_genef));

  // likelihood
  phenotype ~ normal(grand_mean + cov_mef[cov] + L*u, grand_sig2 * exp(cov_vef[cov]));
}
//