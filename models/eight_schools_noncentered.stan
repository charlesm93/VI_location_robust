
data {
  int <lower = 0> J;
  array[J] real y;
  array[J] real<lower = 0> sigma;
}

parameters {
  vector[J] z;
  real mu;
  real<lower = 0> tau;
}

transformed parameters{
  vector[J] theta;
  theta = z * tau + mu;
}

model {
  mu ~ normal(5, 3);
  tau ~ normal(0, 5);
  z ~ std_normal();
  y ~ normal(theta , sigma);
}
