
data {
  int<lower=0> N;
  vector[N] x1;
  vector[N] x2;
  array[N] int y;
  int prior_only;
}

transformed data {
  int N_trial = 10;
}

parameters {
  real alpha;
  real beta1;
  real beta2;
}

model {
  // alpha ~ normal(0, 0.1);
  // beta ~ normal(0, 0.1);
  // alpha ~ cauchy(0, 0.1);
  // beta ~ cauchy(0, 0.1);
  // alpha ~ student_t(4, 0, 0.1);
  // beta ~ student_t(4, 0, 0.1);

  alpha ~ double_exponential(0, 0.5);
  beta1 ~ double_exponential(0, 0.5);
  beta2 ~ double_exponential(0, 0.5);

  // y ~ binomial_logit(N_trial, alpha + beta1 * x1 + beta2 * x2);
  
  if (!prior_only) y ~ bernoulli_logit(alpha + beta1 * x1 + beta2 * x2);
}
