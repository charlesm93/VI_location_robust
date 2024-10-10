
parameters {
  vector[2] y;
}

model {
  target += log_sum_exp(log(0.5) + normal_lpdf(y | -1, 2),
                        log(0.5) + normal_lpdf(y | 3, 1));
}
