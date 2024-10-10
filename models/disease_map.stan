
data {
  int n_obs;
  int n_coordinates;
  array[n_obs] int y;
  vector[n_obs] ye;
  array[n_obs] vector[n_coordinates] x;
  real rho_location_prior;
  real rho_scale_prior;
}

transformed data {
  real delta = 1e-8;
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
  vector[n_obs] eta;
}

transformed parameters {
  vector[n_obs] theta;
   {
     matrix[n_obs, n_obs] L_Sigma;
     matrix[n_obs, n_obs] Sigma;
     Sigma = gp_exp_quad_cov(x, alpha, rho);
     for (n in 1:n_obs) Sigma[n, n] = Sigma[n,n] + delta;
     L_Sigma = cholesky_decompose(Sigma);
     theta = L_Sigma * eta;
   }
}

model {
  rho ~ inv_gamma(rho_location_prior, rho_scale_prior);
  alpha ~ inv_gamma(10, 10);  // CHECK -- is this a good prior?

  eta ~ normal(0, 1);
  y ~ poisson_log(log(ye) + theta);
}
