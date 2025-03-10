
.libPaths("~/Rlib/")
library(ggplot2)
library(distributional)


y1 <- c(-1, 0, 20)
y2 <- c(-1, 0, 1)
y3 <- rnorm(50, mean = 0, sd = 15)
y4 <- rnorm(1000, mean = 0, sd = 15)


likelihood <- function (y, mu, likelihood_names = "Gaussian") {
  # Evaluate log likelihood based on provided name
  ll <- 0
  
  if (likelihood_names == "Gaussian") {
    ll <- sum(dnorm(y, mean = mu, log = TRUE))
  } else if (likelihood_names == "Cauchy") {
    ll <- sum(dcauchy(y, location = mu, log = TRUE))
  } else if (likelihood_names == "Laplace") {
    ll <- - 0.5 * sum(abs(mu - y))
  } else if (likelihood_names == "Logistic") {
    ll <- sum(dlogis(y, location = mu, log = TRUE))
  }
  
  ll
}

# NOTE: logistic is numerically unstable
likelihood_names <- c("Gaussian", "Cauchy", "Laplace", "Logistic")

mu_range <- seq(from = -30, 30, by = 0.01)
l_v <- array(NA, dim = c(length(mu_range), length(likelihood_names), 4))

for (j in 1:length(likelihood_names)) {
  for (i in 1:length(mu_range)) {
    l_v[i, j, 1] <- likelihood(y1, mu_range[i], likelihood_names = likelihood_names[j])
    l_v[i, j, 2] <- likelihood(y2, mu_range[i], likelihood_names = likelihood_names[j])
    l_v[i, j, 3] <- likelihood(y3, mu_range[i], likelihood_names = likelihood_names[j])
    l_v[i, j, 4] <- likelihood(y4, mu_range[i], likelihood_names = likelihood_names[j])
  }
}

name_v <- rep(rep(likelihood_names, each = length(mu_range)), 4)
mu_v <- rep(rep(mu_range, length(likelihood_names)), 4)
y_regime_names <- c("x = (1, 2, 20)", "x = (1, 2, 3)", "N = 10", "N = 1000")
y_regime <- rep(y_regime_names, each = length(likelihood_names) *
                  length(mu_range))

plot.data <- data.frame(likelihood = name_v, mu = mu_v, l = c(l_v), y_regime)

# Order fators to order plots
plot.data$likelihood <-
  factor(plot.data$likelihood, levels = likelihood_names)

plot.data$y_regime <-
  factor(plot.data$y_regime, levels = y_regime_names)



p <- ggplot(data = plot.data, aes(x = mu, y = l, color = likelihood)) +
  geom_line(linewidth = 1) + theme_bw() +
  facet_wrap(y_regime ~ likelihood, scales = "free", ncol = 4) +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 18))
p

