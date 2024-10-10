
data {
  int nu;
}

transformed data {
  int n_dim = 100;
  vector[n_dim] mu = rep_vector(0.0, n_dim);

  matrix[n_dim, n_dim] Sigma;
  for (i in 2:n_dim) {
    for (j in 1:(i -1)) {
      Sigma[i, j] = 1 / 2^abs(i - j);
      Sigma[j, i] = Sigma[i, j];
    }
  }
  for (i in 1:n_dim) Sigma[i, i] = 1;
}

parameters {
  vector[n_dim] z;
}

model {
  z ~ multi_student_t(nu, mu, Sigma);
}
