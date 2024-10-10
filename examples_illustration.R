
## R script for illustrative example.
# Adjust library path, cmdstan path, and working directory to your setting.

.libPaths("~/Rlib/")
library(cmdstanr)
library(ggplot2)
library(scales)
library(VGAM)
library(posterior)
library(bayesplot)
library(sn)
library(matrixStats)
library(latex2exp)

set_cmdstan_path("/Users/charlesm/.cmdstan/cmdstan-2.34.1")
setwd("~/Code/vi_robust")


#####################################################################
## KL for mixture of Gaussian example (Section 3.2, Figure 3)

KL_eval <- function(nu, mu = 10, mc_iter = 1000) {
  integrand = rep(NA, mc_iter)
  z = rnorm(mc_iter, mean = nu)
  
  for (i in 1:mc_iter) {
    integrand[i] = dnorm(z[i], mean = nu, log = TRUE) -
      log(0.5 * dnorm(z[i], mean = mu) + 0.5 * dnorm(z[i], mean = -mu))
  }
  
  mean(integrand)
}

nu_range <- seq(from = -20, to = 20, by = 0.1)
KL_values <- rep(NA, length(nu_range))

for (i in 1:length(nu_range)) KL_values[i] <- KL_eval(nu_range[i])

plot.data <- data.frame(nu = nu_range, KL = KL_values)
p <- ggplot(data = plot.data, aes(x = nu, y = KL)) +
  geom_line() + theme_bw()
p

KL_values_close <- rep(NA, length(nu_range))
for (i in 1:length(nu_range)) {
  KL_values_close[i] <- KL_eval(nu_range[i], mu = 1)
}

# Plot KL divergences
plot.data <- data.frame(nu = rep(nu_range, 2),
                        KL = c(KL_values, KL_values_close),
                        mu = paste0("m=",
                                    rep(c(10, 1), each = length(nu_range))))

cbPalette <- c("#E69F00", "#56B4E9")

p <- ggplot(data = plot.data, aes(x = nu, y = KL, color = mu)) +
  geom_line(linewidth = 1.5) + theme_bw() +
  theme(text = element_text(size = 20)) + ylab(TeX("$KL(q_{\\nu}||p)$"))  +
  xlab(TeX("$\\nu$")) +
  theme(legend.position = "none") +
  scale_color_manual(values= cbPalette)
p


# Plot densities
dmixture <- function(z, mu) {
  0.5 * dnorm(z, mean = mu) + 0.5 * dnorm(z, mean = -mu)
}

z_grid <- seq(from = -15, to = 15, by = 0.1)
f_mixture_close <- rep(NA, length(z_grid))
f_mixture_far <- rep(NA, length(z_grid))

for (i in 1:length(z_grid)) {
  f_mixture_close[i] <- dmixture(z_grid[i], mu = 1)
  f_mixture_far[i] <- dmixture(z_grid[i], mu = 10)
}

plot.data <- data.frame(f = c(f_mixture_close, f_mixture_far),
                        z = rep(z_grid, 2),
                        m <- rep(c("m = 1", "m = 10"), each = length(z_grid)))


p <- ggplot(data = plot.data, aes(x = z, y = f, color = m)) +
  geom_line(linewidth=1.5) + theme_bw() + 
  theme(text = element_text(size = 20),
        axis.text.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) + ylab("p(z)") +
  scale_color_manual(values= cbPalette)
  
p

#####################################################################
## "Tail doesn't matter" plot (Section 3.2, Figure 2b)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
               "#999999")

x <- seq(from = -5, to = 5, by = 0.1)
p_gaussian <- dnorm(x)
p_cauchy <- dcauchy(x)
p_student <- dt(x, df = 10)
p_laplace <- dlaplace(x)

dist_name <- c("VI", "Laplace", "student-t", "Cauchy")
dist <- rep(dist_name, each = length(x))
dist <- factor(dist, levels = dist_name) 

plot.data <- data.frame(x = rep(x, length(dist)), 
                        y = c(p_gaussian, p_laplace, p_student, p_cauchy),
                        Distribution = dist,
                        alpha = c(rep(1, length(x)),
                                  rep(0.5, (length(dist) - length(x))))
                        )

p <- ggplot(data = plot.data, aes(x = x, y = y, color = Distribution,
                                  alpha = alpha)) +
  geom_line(linewidth = 1.25) + theme_bw() +
  scale_color_manual(values= cbPalette) +
  scale_alpha_identity() +
  ylab("p(z)") + xlab("z") + theme(legend.position = c(0.225, 0.7),
                                   legend.title=element_blank())  +
  theme(text = element_text(size = 20))
p

###############################################################################
## "Tail and factorization don't matter" (Section 3.2, Figure 2a)
stan_data <- "data/student_t.json"
mod <- cmdstan_model("models/multi_student_hc.stan")
fit_mcmc <- mod$sample(data = stan_data, 
                       chains = 4, parallel_chains = 4,
                       iter_warmup = 1e3, iter_sampling = 1e4, 
                       seed = 123, adapt_delta = 0.8)

mcmc_scatter(fit_mcmc$draws(), pars = c("z[1]", "z[2]"))

fit_vi <- mod$variational(data = stan_data, seed = 123,
                          algorithm = "meanfield",
                          elbo_samples = 1e3,
                          grad_samples = 1e3,
                          iter = 1e4, output_samples = 1e4)

draws_mcmc_raw <- fit_mcmc$draws()[, , 2:3]
draws_mcmc <- array(NA, dim = dim(draws_mcmc_raw)[c(1, 3)])
draws_mcmc[, 1] <- draws_mcmc_raw[, 1, 1]
draws_mcmc[, 2] <- draws_mcmc_raw[, 1, 2]

draws_vi <- fit_vi$draws()[, 3:4]

plot.data <- data.frame(
  z1 = c(draws_vi[, 1], draws_mcmc[, 1]),
  z2 = c(draws_vi[, 2], draws_mcmc[, 2]),
  algo = c(rep("q", nrow(draws_vi)), rep("p", nrow(draws_mcmc)))
)

plot.data2 <- data.frame(
  mean_z1 = c(mean(draws_vi[, 1]), mean(draws_mcmc[, 1])),
  mean_z2 = c(mean(draws_vi[, 2]), mean(draws_mcmc[, 2])),
  algo = c("q", "p")
)

p <- ggplot(data = plot.data, aes(x = z1, y = z2, color = algo)) +
  geom_density_2d(alpha = 0.4) + geom_point(alpha = 0.01) + theme_bw() +
  geom_point(data = plot.data2, aes(x = mean_z1, y = mean_z2,
                                    color = algo), size = 5,
             shape = 10) +
  xlim(-2, 2) + ylim(-2, 2) +
  xlab(TeX("$z_1$")) + ylab(TeX("$z_2$")) +
  theme(text = element_text(size = 20)) +
  scale_color_discrete(name = "")
  
p

###############################################################################
## Elliptical symmetry matters (Section 4.2, Figure 6)
mod <- cmdstan_model("models/crescent2.stan")
fit_mcmc <- mod$sample(chains = 1,
                       iter_warmup = 2e3, iter_sampling = 2e4, 
                       seed = 123, adapt_delta = 0.99)

mcmc_scatter(fit_mcmc$draws(), pars = c("z[1]", "z[2]"))

fit_vi <- mod$variational(seed = 123,
                          algorithm = "fullrank",
                          elbo_samples = 1e3,
                          grad_samples = 1e3,
                          iter = 1e4, output_samples = 1e4)


draws_mcmc_raw <- fit_mcmc$draws()[, , 2:3]
draws_mcmc <- array(NA, dim = dim(draws_mcmc_raw)[c(1, 3)])
draws_mcmc[, 1] <- draws_mcmc_raw[, 1, 1]
draws_mcmc[, 2] <- draws_mcmc_raw[, 1, 2]


draws_vi <- fit_vi$draws()[, 3:4]

stan_data <- "data/student_t.json"
mod2 <- cmdstan_model("models/multi_student_hc.stan")
fit2_mcmc <- mod2$sample(data = stan_data, chains = 1,
                       iter_warmup = 2e3, iter_sampling = 2e4, 
                       seed = 123, adapt_delta = 0.99)

mcmc_scatter(fit2_mcmc$draws(), pars = c("z[1]", "z[2]"))

fit2_vi <- mod2$variational(data = stan_data, seed = 123,
                            algorithm = "fullrank",
                            elbo_samples = 1e3,
                            grad_samples = 1e3,
                            iter = 1e4, output_samples = 1e4)

draws2_mcmc_raw <- fit2_mcmc$draws()[, , 2:3]
draws2_mcmc <- array(NA, dim = dim(draws2_mcmc_raw)[c(1, 3)])
draws2_mcmc[, 1] <- draws2_mcmc_raw[, 1, 1]
draws2_mcmc[, 2] <- draws2_mcmc_raw[, 1, 2]

draws2_vi <- fit2_vi$draws()[, 3:4]

plot.data <- data.frame(
  z1 = c(draws_vi[, 1], draws_mcmc[, 1], draws2_vi[, 1], draws2_mcmc[, 1]),
  z2 = c(draws_vi[, 2], draws_mcmc[, 2], draws2_vi[, 2], draws2_mcmc[, 2]),
  algo = rep(c(rep("q", nrow(draws_vi)), rep("p", nrow(draws_mcmc))), 2),
  target = rep(c("crescent", "student-t"), each = nrow(draws_vi) + nrow(draws_mcmc))
)

plot.data2 <- data.frame(
  mean_z1 = c(mean(draws_vi[, 1]), mean(draws_mcmc[, 1]),
              mean(draws2_vi[, 1]), mean(draws2_mcmc[, 1])
              ),
  mean_z2 = c(mean(draws_vi[, 2]), mean(draws_mcmc[, 2]),
              mean(draws2_vi[, 2]), mean(draws2_mcmc[, 2])
              ),
  algo = rep(c("q", "p"), 2),
  target = rep(c("crescent", "student-t"), each = 2)
)

plot.data.crescent <- plot.data[plot.data$target == "crescent", ]
plot.data2.crescent <- plot.data2[plot.data2$target == "crescent", ]


p <- ggplot(data = plot.data.crescent, aes(x = z1, y = z2, color = algo)) +
  geom_density_2d(alpha = 0.4) + geom_point(alpha = 0.005) + theme_bw() +
  geom_point(data = plot.data2.crescent, aes(x = mean_z1, y = mean_z2,
                                    color = algo), size = 5,
             shape = 10) +
  xlim(-20, 20) + ylim(-5, 5) +
  xlab(TeX("$z_1$")) + ylab(TeX("$z_2$")) +
  theme(text = element_text(size = 20)) +
  scale_color_discrete(name = "mean") +
  facet_wrap(~target, scale = "free") +
  theme(legend.position = "none") + ylab(" ")
p

plot.data.student <- plot.data[plot.data$target == "student-t", ]
plot.data2.student <- plot.data2[plot.data2$target == "student-t", ]

p <- ggplot(data = plot.data.student, aes(x = z1, y = z2, color = algo)) +
  geom_density_2d(alpha = 0.4) + geom_point(alpha = 0.005) + theme_bw() +
  geom_point(data = plot.data2.student, aes(x = mean_z1, y = mean_z2,
                                             color = algo), size = 5,
             shape = 10) +
  xlim(-2, 2) + ylim(-3, 3) +
  xlab(TeX("$z_1$")) + ylab(TeX("$z_2$")) +
  theme(text = element_text(size = 20)) +
  scale_color_discrete(name = "mean") +
  facet_wrap(~target, scale = "free") +
  theme(legend.position = c(0.1, 0.8))
p

###############################################################################
## Skewed normal target (Section 4.2, Figure 4)
# In this example, approximate a skewed normal distribution with a Laplace
# distribution.
# NOTE: code for skewed normal seems poorly written (log density is expensive!)

KL_eval <- function(nu, mu = 0, alpha = 0, mc_iter = 2e4) {
  integrand = rep(NA, mc_iter)
  z = rlaplace(mc_iter, location = nu, scale = 1)
  
  for (i in 1:mc_iter) {
    integrand[i] = dlaplace(z[i], location = nu, log = TRUE) -
      dsn(z[i], xi=mu, omega=1, alpha=alpha, log=TRUE)
  }
  
  mean(integrand)
}

nu_range <- seq(from = -0.5, to = 2.5, by = 0.025)
alpha_range <- seq(from = 0, to = 5, by = 1)


nu_star <- rep(NA, length(alpha_range))

for (j in 1:length(alpha_range)) {
  KL_values <- rep(NA, length(nu_range))
  alpha = alpha_range[i]
  for (i in 1:length(nu_range)) KL_values[i] <- KL_eval(nu_range[i],
                                                        alpha = alpha_range[j])
  nu_star[j] <- nu_range[which.min(KL_values)]
}

n_sim <- 1e5
z_p <- array(NA, c(n_sim, length(alpha_range)))
z_q <- array(NA, c(n_sim, length(alpha_range)))

for (j in 1:length(alpha_range)) {
  z_p[, j] <- rsn(n_sim, xi = 0, omega = 1, alpha = alpha_range[j])
  z_q[, j] <- rlaplace(n_sim, location = nu_star[j])
}

mean_p <- colMeans(z_p)
mean_q <- colMeans(z_q)

alpha_tex <- paste0(TeX("$\\alpha = $ "), alpha_range)

plot.data <- data.frame(z = c(z_p, z_q),
                        dist = rep(c("p", "q"), each = length(alpha_range) * n_sim),
                        alpha = as.character(rep(rep(factor(alpha_range), each = n_sim), 2)),
                        mean_p = rep(rep(mean_p, each = n_sim), 2),
                        mean_q = rep(rep(mean_q, each = n_sim), 2)
                        )


sub.plot.data = plot.data[plot.data$alpha %in% c(0, 3, 5), ]

appender <- function(string) {
  TeX(paste("$\\alpha = $", string))
}

p <- ggplot(data = sub.plot.data, 
            aes(x = z, color = dist)) +
  geom_density(alpha = 0.25, linewidth = 1) + theme_bw() +
  geom_vline(data = sub.plot.data, 
             aes(xintercept = mean_p), linetype = "solid",
             color = "#F8766D") +
  geom_vline(data = sub.plot.data,
             aes(xintercept = mean_q), linetype = "solid",
             color = "#00BFC4") +
  facet_wrap(~alpha, nrow = 1,
             labeller = as_labeller(appender,
                                    default = label_parsed),
             ) +
  theme(text = element_text(size = 22)) +
  xlim(-2.5, 6)
p


plot.data <- data.frame(z = c(z_p, z_q), dist = rep(c("p", "q"), each = n_sim))

p <- ggplot(data = plot.data, aes(x = z, fill = dist, color = dist)) +
  geom_histogram(color = "black", alpha = 0.5, position="identity", bins = 100) + 
  theme_bw()
p

p <- ggplot(data = plot.data, aes(x = z, color = dist)) +
  geom_density(alpha = 0.5, linewidth = 2) + theme_bw() +
  geom_vline(aes(xintercept = mean(z_p)), linetype = "dashed") +
  theme(text = element_text(size = 22))
p

mean(z_p)
mean(z_q)

n_mc_sim <- c(1e2, 1e3, 1e4)  #, 1e5)
KL_values <- array(NA, c(length(nu_range), length(n_mc_sim)))

for (j in 1:length(n_mc_sim)) {
  for (i in 1:length(nu_range)) KL_values[i, j] <- KL_eval(nu_range[i], alpha = 0,
                                                           mc_iter = n_mc_sim[j])  
}


## Additional plot to show how tricky noisy estimation of the KL can be.

plot.data <- data.frame(KL = c(KL_values),
                        n_MC = as.factor(rep(n_mc_sim, each = length(nu_range))),
                        nu = rep(nu_range, length(n_mc_sim))
)

p <- ggplot(data = plot.data, aes(x = nu, y = KL, color = n_MC)) +
  geom_line(linewidth=1.5, alpha = 0.85) + theme_bw() + 
  theme(text = element_text(size = 22)) +
  theme(legend.position = c(0.87, 0.8)) +
  ylab("Estimate of KL(q||p)")
p

#####################################################################
## Logistic regression (Section 2.3, Figure 1 and Appendix B2, Figure 9)

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

eval_log_joint <- function(z, stan_data) {
  log_joint = dlaplace(z[1], scale = 0.1, log = TRUE)
  log_joint = log_joint + dlaplace(z[2], scale = 0.1, log = TRUE)
  log_joint = log_joint + dlaplace(z[3], scale = 0.1, log = TRUE)
  
  log_joint = log_joint + sum(dbinom(stan_data$y, 1, 
                                     inv_logit(z[1] + stan_data$x1 * z[2] +
                                               stan_data$x2 * z[3]),
                                     log = TRUE))
  
  log_joint
}


check_sym <- function(eval_log_joint, symm_point, z) {
  # Returns the fractional difference between the log joint
  # evaluated at z and z'.
  eval_z <- eval_log_joint(z, stan_data)

  abs((eval_z - eval_log_joint(-z + 2 * symm_point, stan_data)) / eval_z)
}


n_sim_vec <- c(0, 1, 2, 4, 8, 16, 32, 128)
n_sim_vi <- 1
n_chains <- 4
n_sample <- 5e3

alpha_mc <- rep(NA, length(n_sim_vec))
beta1_mc <- rep(NA, length(n_sim_vec))
beta2_mc <- rep(NA, length(n_sim_vec))

alpha_vi <- array(NA, c(n_sim_vi, length(n_sim_vec)))
beta1_vi <- array(NA, c(n_sim_vi, length(n_sim_vec)))
beta2_vi <- array(NA, c(n_sim_vi, length(n_sim_vec)))
max_err <- array(NA, c(n_sim_vi, length(n_sim_vec)))
check_sym_eval <- array(NA, c(n_sim_vi, length(n_sim_vec)))
check_sym_mc <- rep(NA, length(n_sim_vec))

draws_mc_vec <- array(NA, c(length(n_sim_vec), n_chains * n_sample, 3))
draws_vi_vec <- array(NA, c(length(n_sim_vec), n_chains * n_sample, 3))


mod <- cmdstan_model("models/logistic2.stan")

set.seed(1958)
x1 <- rnorm(max(n_sim_vec), sd = 5)
x2 <- rnorm(max(n_sim_vec), sd = 5)
alpha <- -1
beta1 <- 0.5
beta2 <- 0.5
N <- 10

n_dim <- 3


y <- rbinom(max(n_sim_vec), 1, inv_logit(alpha + beta1 * x1 + beta2 * x2))
hist(inv_logit(alpha + beta1 * x1 + beta2 * x2), breaks = 30)

for (i in 1:length(n_sim_vec)) {
  print(paste0(i - 1, " / ", length(n_sim_vec)))
  n_sim <- n_sim_vec[i]
  n_sim_max <- max(c(n_sim, 1))
  
  stan_data <- list(N = n_sim_max, x1 = x1[1:n_sim_max], x2 = x2[1:n_sim_max],
                    y = y[1:n_sim_max], prior_only = 0)
  
  if (n_sim == 0) stan_data$prior_only = 1
  
  fit <- mod$sample(data = stan_data, seed = 123,
                    chains = 4, parallel_chains = 4,
                    iter_sampling = 5e3)
  mc_summary <- fit$summary()
  alpha_mc[i] <- mc_summary[2, 2]$mean
  beta1_mc[i] <- mc_summary[3, 2]$mean
  beta2_mc[i] <- mc_summary[4, 2]$mean
  
  alpha_sd <- mc_summary[2, 4]$sd
  beta1_sd <- mc_summary[3, 4]$sd
  beta2_sd <- mc_summary[4, 4]$sd
  
  if (FALSE) {
    mcmc_hist(fit$draws(c("alpha", "beta1", "beta2")))
  }
  
  draws_mc <- fit$draws()
  draws_mc <- array(draws_mc[, , 2:4], dim = c(n_sample * n_chains, 3))
  draws_mc_vec[i, , ] <- draws_mc
  
  check_sym_mc_local <- rep(NA, n_sample * n_chains)

  for (k in 1:(n_sample * n_chains)) {
    z <- c(draws_mc[k, ])
    check_sym_mc_local[k] <- check_sym(eval_log_joint,
                                       c(alpha_mc[i], beta1_mc[i], beta2_mc[i]),
                                       z)
  }
  check_sym_mc[i] <- sort(check_sym_mc_local)[0.5 * n_sample * n_chains]
  
  for (j in 1:n_sim_vi) {
    fit_vi <- mod$variational(data = stan_data, seed = 2004 + j, iter=1e4,
                              elbo_samples = 1e4, grad_samples = 1e4,
                              draws = n_chains * n_sample,
                              algorithm = "fullrank",
                              init = list(list(alpha = alpha_mc[i],
                                               beta1 = beta1_mc[i],
                                               beta2 = beta2_mc[i]))
                              )
    
    alpha_vi_local <- fit_vi$summary()[3, 2]$mean
    beta1_vi_local <- fit_vi$summary()[4, 2]$mean
    beta2_vi_local <- fit_vi$summary()[5, 2]$mean
    
    alpha_vi[j, i] <- abs(alpha_vi_local - alpha_mc[i]) / max(abs(alpha_mc[i]), alpha_sd)
    beta1_vi[j, i] <- abs(beta1_vi_local - beta1_mc[i]) / max(abs(beta1_mc[i]), beta1_sd)
    beta2_vi[j, i] <- abs(beta2_vi_local - beta2_mc[i]) / max(abs(beta2_mc[i]), beta2_sd)
    
    max_err[j, i] <- mean(c(alpha_vi[j, i], beta1_vi[j, i], beta2_vi[j, i]))
    
    draws <- fit_vi$draws()
  }
  
  
  draws_vi <- matrix(draws[, 3:5], nrow = nrow(draws))
  draws_vi_vec[i, , ] <- draws_vi
}

plot.data <- data.frame(n_sim = rep(as.factor(rep(n_sim_vec, each = n_sim_vi)), n_dim),
                        err = c(c(alpha_vi), c(beta1_vi), c(beta2_vi)),
                        parameter = rep(c("alpha", "beta1", "beta2"),
                                        each = length(n_sim_vec) * n_sim_vi),
                        check_sym = rep(c(check_sym_mc), n_dim))


p <- ggplot(data = plot.data, aes(x = n_sim, y = err, color = parameter)) + # , color = parameter)) +
  geom_boxplot() + theme_bw() +
  xlab("number of observations") +
  theme(legend.position = c(0.87, 0.8)) +
  theme(text = element_text(size = 22))
p


plot.data_sub = plot.data[plot.data$parameter == "alpha", ]
p <- ggplot(data = plot.data_sub, aes(x = n_sim, y = err)) +
  geom_point(shape = 19, size = 3) + theme_bw() +
  xlab("number of observations") +
  ylab("scaled error") +
  theme(legend.position = c(0.87, 0.8)) +
  theme(text = element_text(size = 22)) +
  theme(legend.position="none")
  
p



plot.data <- data.frame(
  draws = c(draws_mc_vec, draws_vi_vec),
  n_sim = rep(paste0("N = ", n_sim_vec), n_chains * n_sample * n_dim * 2),
  param = rep(rep(c("alpha", "beta1", "beta2"),
                  each = n_chains * n_sample * length(n_sim_vec)), 2),
  algo = rep(c("MCMC", "VI"), each = n_chains * n_sample * n_dim * length(n_sim_vec))
)

plot.data$n_sim <- factor(plot.data$n_sim, levels = paste0("N = ", n_sim_vec))

mean_draws <- c()
median_draws <- c()
mean_vi_draws <- c()
for (i in 1:length(n_sim_vec)) {
  mean_draws <- c(mean_draws, colMeans(draws_mc_vec[i, , ]))
  median_draws <- c(median_draws, colMedians(draws_mc_vec[i, , ]))
  mean_vi_draws <- c(mean_vi_draws, colMeans(draws_vi_vec[i, , ]))
}

plot.data2 <- data.frame(
  mean = mean_draws,
  mean_vi = mean_vi_draws,
  median = median_draws,
  param = rep(c("alpha", "beta1", "beta2"), length(n_sim_vec)),
  n_sim = rep(paste0("N = ", n_sim_vec), each = n_dim)
)

plot.data2$n_sim <- factor(plot.data2$n_sim, levels = paste0("N = ", n_sim_vec))

linewidth_mean = 0.5

p <- ggplot(data = plot.data, aes(x = draws, color = algo)) +
  geom_density(linewidth = 1) +
  geom_vline(data = plot.data2, linewidth = linewidth_mean, linetype = "solid",
             color = "#00BFC4",aes(xintercept = mean)) +
  geom_vline(data = plot.data2, linewidth = linewidth_mean, linetype = "solid",
             color = "#F8766D", aes(xintercept = mean_vi)) +
  facet_wrap(~n_sim+param, scales = "free", ncol = 3) + theme_bw() +
  theme(legend.title = element_blank())
p


# Use a subset of the data
lim_x <- list(
  First = scale_x_continuous(limits = c(-1.5, 1.5)),
  Second = scale_x_continuous(limits = c(-1.5, 1.5))
)

n_sub <- c(0, 4, 128)
parm_sub <- "alpha"

plot.data_sub <- plot.data[plot.data$n_sim %in% n_sub &
                             plot.data$param == parm_sub, ]
plot.data_sub$N <- factor(paste0("N = ", plot.data_sub$n_sim),
                          levels = paste0("N = ", n_sub))

plot.data2_sub <- plot.data2[plot.data2$n_sim %in% n_sub &
                               plot.data2$param == parm_sub, ]
plot.data2_sub$N <- factor(paste0("N = ", plot.data2_sub$n_sim),
                           levels = paste0("N = ", n_sub))

p <- ggplot(data = plot.data_sub, aes(x = draws, color = algo)) +
  geom_density(linewidth = 1) +
  geom_vline(data = plot.data2_sub, linewidth = linewidth_mean, linetype = "solid",
             color = "#F8766D",aes(xintercept = mean)) +
  geom_vline(data = plot.data2_sub, linewidth = linewidth_mean, linetype = "solid",
             color = "#00BFC4", aes(xintercept = mean_vi)) +
  facet_wrap(~N, scales = "free", nrow = 1) + theme_bw() +
  theme(text = element_text(size = 22)) + xlab(" ") + ylab(" ") +
  xlim(-1.5, 1.5) +
  theme(legend.position = c(0.9, 0.8),
        legend.title=element_blank()) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
p


###############################################################################
## Experiment on scale matrix (Section 4.2, Figure 5)

mod <- cmdstan_model("models/multi_student_highd.stan")

# nu_range <- seq(from = 2, to = 4, by = 2)
# nu_ramge <- seq(from = 10, to = 12, by = 2)
nu_range <- c(2, 4, 6, 8, 10, 20)
nu_range_l <- length(nu_range)

n_dim <- 100

cov_vi <- array(NA, c(nu_range_l, n_dim, n_dim))
cor_vi <- array(NA, c(nu_range_l, n_dim, n_dim))

for (i in 1:length(nu_range)) {
  stan_data <- list(nu = nu_range[i])

  fit_vi <-  mod$variational(
    data = stan_data,
    elbo_samples = 1e3, grad_samples = 1e3,
    draws = 1e5,
    algorithm = "fullrank")
  
  cov_vi[i, , ] <- cov(fit_vi$draws(variables = c("z")))
  cor_vi[i, , ] <- cor(fit_vi$draws(variables = c("z")))
}

S <- matrix(NA, nrow = n_dim, ncol = n_dim)
diag(S) <- 1
for (i in 1:n_dim) {
  for (j in 1:(i - 1)) {
    S[i, j] = 1 / 2^(abs(i - j))
    S[j, i] = S[i, j]
  }
}

R_vec <- rep(NA, nu_range_l)
lambda_max_vec <- rep(NA, nu_range_l)
lambda_min_vec <- rep(NA, nu_range_l)

for (i in 1:nu_range_l) {
  R_vec[i] <- 0.5 * (cov_vi[i, 1, 1] + cov_vi[i, 2, 2]) / cov_vi[i, 1, 2]
  lambda_max_vec[i] <- max(cov_vi[i, , ] / S)
  lambda_min_vec[i] <- min(cov_vi[i, , ] / S)
}

## New plots for higher dimensional target.

S_vec <- array(NA, c(nu_range_l, n_dim, n_dim))
for (i in 1:nu_range_l) {
  S_vec[i, , ] <- S
}


plot.data <- data.frame(
  cov_vi = c(cov_vi),
  cor_vi = c(cor_vi),
  k = rep(nu_range, nrow(S)^2),
  M = c(S_vec) 
)

plot.data_sub <- plot.data[plot.data$k %in% c(2, 6, 10, 20), ] 
plot.data_sub$k <- as.factor(plot.data_sub$k)

p <- ggplot(data = plot.data_sub, aes(x = M, y = cov_vi, color = k)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  geom_smooth(method = "lm", linetype = "dashed", linewidth = 0.5) + 
  theme_bw() +
  xlab(TeX("$M_{ij}$")) + ylab(TeX("$S_{ij}$")) +
  theme(text = element_text(size = 22)) +
  theme(legend.position = c(0.1, 0.7)) +
  geom_abline(slope=1, linetype="solid")
p

p <- ggplot(data = plot.data_sub, aes(x = M, y = cor_vi, color = k)) +
  geom_point(shape = 1, size = 5) + theme_bw() + # geom_smooth(method = "lm") + 
  theme_bw() +
  xlab(TeX("$\\rho_{ij}$ (True)")) + ylab(TeX("$\\rho_{ij}$ (VI)")) +
  theme(text = element_text(size = 22)) +
  theme(legend.position = c(0.1, 0.7)) +
  geom_abline(slope=1, linetype="solid")

p
