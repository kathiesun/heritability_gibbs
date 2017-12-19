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
  real<lower = 0> grand_log_sig;
  real<lower = 0.01> sigma_covmef;
  real<lower = 0.01> sigma_covvef;
  real<lower = 0.01> sigma2_genef; 
  //real<lower = 0.01> sigma_genef; 


  vector[num_cov] cov_mef;
  vector[num_cov] cov_vef;
  vector[N] u;
}

transformed parameters {
  vector[num_cov] h2;
  real<lower = 0.01> sum_sig2_1;
  real<lower = 0.01> sum_sig2_2;
  sum_sig2_1 = pow(exp(grand_log_sig + cov_vef[1]), 2) + sigma2_genef;
  sum_sig2_2 = pow(exp(grand_log_sig + cov_vef[2]), 2) + sigma2_genef;
  //sum_sig2_1 = pow(exp(grand_log_sig + cov_vef[1]), 2) + pow(sigma_genef, 2);
  //sum_sig2_2 = pow(exp(grand_log_sig + cov_vef[2]), 2) + pow(sigma_genef, 2);
  h2[1] = sigma2_genef/sum_sig2_1;
  h2[2] = sigma2_genef/sum_sig2_2;
  //h2[1] = pow(sigma_genef,2)/sum_sig2_1;
  //h2[2] = pow(sigma_genef,2)/sum_sig2_2;
}

model {
  // priors
  grand_mean     ~ normal(0, 1);
  grand_log_sig  ~ normal(0, 1);
  sigma_covmef  ~ cauchy(0, 5);
  sigma_covvef  ~ cauchy(0, 5);
  sigma2_genef   ~ inv_gamma(1, 1);
  //sigma_genef   ~ cauchy(0, 5);

  // sex and genetic effects
  cov_mef ~ normal(0, sigma_covmef);
  cov_vef ~ normal(0, sigma_covvef);
  u       ~ normal(0, sqrt(sigma2_genef));
  //u       ~ normal(0, sigma_genef);

  // likelihood
  phenotype ~ normal(grand_mean + cov_mef[cov] + L*u, exp(grand_log_sig + cov_vef[cov]));
}
//


