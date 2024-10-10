
data {
  int <lower = 0> J;
  array[J] real y;
  array[J] real<lower = 0> sigma;
}

parameters {
  vector[J] theta;
  real mu;
  real<lower = 0> tau;
}

model {
  mu ~ normal(5, 3);
  tau ~ normal(0, 5);
  theta ~ normal(mu, tau);
  y ~ normal(theta , sigma);
}
