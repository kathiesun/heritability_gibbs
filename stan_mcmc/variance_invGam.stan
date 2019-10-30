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
  real<lower = 0> sigma2;
  real<lower = 0.01> sigma2_covmef;
  real<lower = 0.01> sigma2_covvef;
  real<lower = 0.01> sigma2_genef; 

  vector[num_cov] cov_mef;
  vector[num_cov] cov_vef;
  vector[N] u;
}

transformed parameters {
  vector[num_cov] h2;
  real<lower = 0.01> sum_sig2_1;
  real<lower = 0.01> sum_sig2_2;
  sum_sig2_1 = sigma2 + cov_vef[1] + sigma2_genef;
  sum_sig2_2 = sigma2 + cov_vef[2] + sigma2_genef;
  h2[1] = sigma2_genef/sum_sig2_1;
  h2[2] = sigma2_genef/sum_sig2_2;
}

model {
  // priors
  grand_mean     ~ normal(0, 1);
  sigma2         ~ inv_gamma(1, 1);
  sigma2_covmef  ~ inv_gamma(1, 1);
  sigma2_covvef  ~ inv_gamma(1, 1);
  sigma2_genef   ~ inv_gamma(1, 1);

  // sex and genetic effects
  cov_mef ~ normal(0, sqrt(sigma2_covmef));
  cov_vef ~ normal(0, sqrt(sigma2_covvef));
  u       ~ normal(0, sqrt(sigma2_genef));

  // likelihood
  phenotype ~ normal(grand_mean + cov_mef[cov] + L*u, sqrt(sigma2) + cov_vef[cov]);
}
//
//




