
parameters {
  vector[2] z;
}

model {
  z[1] ~ normal(0, 10);
  z[2] ~ normal(0.03 * (z[1]^2 - 100), 1);
}
