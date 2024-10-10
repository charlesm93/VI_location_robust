
data {
  int nu;
}

transformed data {
  vector[2] mu = to_vector({0.0, 0.0});

  matrix[2, 2] Sigma;
  Sigma[1, 1] = 1;
  Sigma[1, 2] = 0.5;
  Sigma[2, 1] = 0.5;
  Sigma[2, 2] = 1;
}

parameters {
  vector[2] z;
}

model {
  z ~ multi_student_t(nu, mu, Sigma);
}
