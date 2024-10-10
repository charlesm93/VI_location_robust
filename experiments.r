
## Numerical experiment (Section 5) for robust VI manuscript
# Adjust library path, cmdstan path and birdgestan path to your directory.

rm(list = ls())
gc()
.libPaths("~/Rlib/")
library(cmdstanr)
library(jsonlite)
library(bridgestan)
library(ggplot2)
library(latex2exp)
library(posterior)
library(bayesplot)

setwd("~/Code/vi_robust")
set_cmdstan_path("/Users/charlesm/.cmdstan/cmdstan-2.34.1")
bridge_stan_path <- "../bridgestan/"


#####################################################################
## Functions to run experiment

bridgestan_model <- function(stan_file, data_file, bridge_stan_path) {
  current_wd <- getwd()
  
  file_dir <- paste0(bridge_stan_path, "experiments/")
  model_dir <- paste0(file_dir, stan_file)
  data_dir <- paste0(file_dir, data_file)
  
  stan_file_name <- substr(stan_file, start = 8, stop = nchar(stan_file))
  file.copy(stan_file, model_dir, overwrite = TRUE)
  data_file_name <- substr(data_file, start = 5, stop = nchar(stan_file))
  file.copy(data_file, data_dir, overwrite = TRUE)

  bridge_model <- StanModel$new(model_dir, data_dir, 123)

  return(bridge_model)
}

check_sym <- function(bridge_model, symm_point, z) {
  # Returns the fractional difference between the log joint
  # evaluated at z and z'.
  eval_z <- bridge_model$log_density(theta_unc = z)
  abs((eval_z - bridge_model$log_density(theta_unc = -z + 2 * symm_point)) 
      / eval_z)
}

colMedian <- function(A) {
  A_col_median <- rep(NA, dim(A)[2])
  for (i in 1:dim(A)[2]) {
    A_col_median[i] <- median(A[, i])
  }
  A_col_median
}


#####################################################################

seed <- 123
batch_size <- 1e3
n_sample <- 5e3
n_chains <- 4

run_experiment <- function(data_file, stan_file, seed = 123, n_sample = 5e3,
                           n_chains = 4, sym_quant = 0.9,
                           vi_algorithm = "fullrank",
                           batch_size = 1e3,
                           draw_array = NULL,
                           run_fits = TRUE) {

  name <- substr(stan_file, start = 8, stop = nchar(stan_file) - 5)
  
  stan_data <- read_json(data_file, simplifyVector = TRUE)
  mod <- cmdstan_model(stan_file)
  
  # Construct bridgestan model to get transformation to unconstrained scale
  bridge_model <- bridgestan_model(stan_file, data_file, bridge_stan_path)
  n_dim <- bridge_model$param_num()
  
  if (is.null(draw_array)) {
    if (run_fits) {
    fit_mcmc <- mod$sample(data = stan_data, chains = n_chains,
                           parallel_chains = n_chains,
                           iter_warmup = 1e3,
                           iter_sampling = n_sample, 
                           seed = seed,
                           adapt_delta = 0.95)
    
    fit_mcmc$save_object(file = file.path("deliv", paste0(name, ".mcmc.fit.RDS")))
    }
    
    fit_mcmc <- readRDS(file = file.path("deliv", paste0(name, ".mcmc.fit.RDS")))
    
    # Get MCMC draws on the unconstrained scale
    draws <- fit_mcmc$draws()[, , 2:(1 + n_dim)]
    draw_array_const <- array(draws, dim = c(n_sample * n_chains, n_dim))
    draw_array <- array(NA, dim = c(n_sample * n_chains, n_dim))
    for (i in 1:dim(draw_array)[1]) {
      draw_array[i, ] <- bridge_model$param_unconstrain(draw_array_const[i, ])
    }
  }
  
  if (run_fits) {
    fit_advi <- mod$variational(data = stan_data, seed = seed,
                                algorithm = vi_algorithm,
                                elbo_samples = batch_size,
                                grad_samples = batch_size,
                                tol_rel_obj = NULL,
                                iter = 1e4, output_samples = 1e4)
    
    fit_advi$save_object(file = file.path("deliv", paste0(name, ".advi.fit.RDS")))
  }
  
  fit_advi <- readRDS(file = file.path("deliv", paste0(name, ".advi.fit.RDS")))
  
  # Now put ADVI draws on the unconstrained scale
  draw_advi_const <- fit_advi$draws()[, 3:(2 + n_dim)]
  draw_advi <- array(NA, dim = dim(draw_advi_const))
  for (i in 1:dim(draw_advi)[1]) {
    draw_advi[i, ] <- bridge_model$param_unconstrain(draw_advi_const[i, ])
  }
  
  cov_mcmc <- cov(draw_array)
  cor_mcmc <- cor(draw_array)
  cov_advi <- cov(draw_advi)
  cor_advi <- cor(draw_advi)
  
  mean_mcmc <- colMeans(draw_array)
  mean_advi <- colMeans(draw_advi)

  # test of symmetry using MCMC or reference draws
  check_sym_vec <- rep(NA, dim(draw_advi)[1])
  for (i in 1:dim(draw_array)[1]) {
      check_sym_vec[i] <- check_sym(bridge_model,
                                    symm_point = mean_advi,
                                    z = draw_array[i, ])
  }
 
  list(mean_err = abs((mean_mcmc - mean_advi) /
                        max(abs(mean_mcmc), sqrt(diag(cov_mcmc)))),
       cor_err = abs(cor_mcmc - cor_advi),
       cov_err = abs(cov_mcmc - cov_advi) / abs(cov_mcmc),
       check_sym = sort(check_sym_vec)[sym_quant * dim(draw_array)[1]],
       n_dim = n_dim)
}

#####################################################################
## Generate samples for illustrative models
draw_mixture <- matrix(NA, nrow = 2e4, ncol = 2)
for (i in 1:2e4) {
  for (j in 1:2) {
    u <- runif(1)
    if (u <= 0.5) {
      draw_mixture[i, j] <- rnorm(1, mean = -1, sd = 2)
    } else {
      draw_mixture[i, j] <- rnorm(1, mean = 3, sd = 1)
    }
  }
}

p <- ggplot(data = data.frame(x = draw_mixture[, 1],
                              y = draw_mixture[, 2]),
            aes(x = x, y = y)) +
  geom_point(alpha = 0.1)
p

draw_crescent <- matrix(NA, nrow = 2e4, ncol = 2)
for (i in 1:2e4) {
  draw_crescent[i, 1] = rnorm(1, sd = 10)
  draw_crescent[i, 2] = rnorm(1, mean = 0.03 * (draw_crescent[i, 1]^2 -100),
                              sd = 1)
}

p <- ggplot(data = data.frame(x = draw_crescent[, 1],
                              y = draw_crescent[, 2]),
            aes(x = x, y = y)) +
  geom_point(alpha = 0.1)
p

######################################################################

exp_crescent <-
  run_experiment(data_file = "data/student_t.json",
                 stan_file = "models/crescent2.stan",
                 draw_array = draw_crescent,
                 batch_size = 1e4,
                 run_fits = FALSE)

exp_student <-
  run_experiment(data_file = "data/student_t.json",
                 stan_file = "models/multi_student.stan",
                 batch_size = 1e2,
                 run_fits = FALSE)


exp_mixture <-
  run_experiment(data_file = "data/student_t.json",
                 stan_file = "models/mixture3.stan",
                 draw_array = draw_mixture,
                 run_fits = FALSE)


exp_eight_schools_nc <- 
  run_experiment(data_file = "data/eight_schools.json",
                 stan_file = "models/eight_schools_noncentered.stan",
                 run_fits = FALSE)

exp_eight_schools_c <-
  run_experiment(data_file = "data/eight_schools.json",
                 stan_file = "models/eight_schools_centered.stan",
                 run_fits = FALSE)

exp_glm <- 
  run_experiment(data_file = "data/GLM_Binomial_data.json",
                 stan_file = "models/GLM_Binomial_model.stan",
                 batch_size = 1e4,
                 run_fits = FALSE)

exp_skim <-
  run_experiment(data_file = "data/prostate_200.json",
                 stan_file = "models/skim_logit.stan",
                 seed = 123,
                 vi_algorithm = "fullrank",
                 batch_size = 50,
                 run_fits = FALSE)  


# this model gives the ADVI optimizer a hard time
exp_disease_map <-
  run_experiment(data_file = "data/disease_100.json",
                 stan_file = "models/disease_map.stan",
                 seed = 123,
                 vi_algorithm = "meanfield",
                 batch_size = 50,
                 run_fits = FALSE)


#####################################################################

model_names <- c("SKIM", "8schools_nc", "8schools", "GLM", "student",
                 "crescent", "mixture", "disease_map")

n_dim_vec <- c(exp_skim$n_dim,
               exp_eight_schools_nc$n_dim,
               exp_eight_schools_c$n_dim,
               exp_glm$n_dim,
               exp_student$n_dim,
               exp_crescent$n_dim,
               exp_mixture$n_dim,
               exp_disease_map$n_dim)

check_sym_vec <- c(exp_skim$check_sym,
                   exp_eight_schools_nc$check_sym,
                   exp_eight_schools_c$check_sym,
                   exp_glm$check_sym,
                   1e-4,  # for studen-t, set to a small value for log scale.
                   # exp_student$check_sym,
                   exp_crescent$check_sym,
                   exp_mixture$check_sym,
                   exp_disease_map$check_sym
                   )

plot.data <- data.frame(
  err = c(
          exp_skim$mean_err, c(exp_skim$cor_err), c(exp_skim$cov_err),
          exp_eight_schools_nc$mean_err, 
          c(exp_eight_schools_nc$cor_err),
          c(exp_eight_schools_nc$cov_err),
          exp_eight_schools_c$mean_err,
          c(exp_eight_schools_c$cor_err),
          c(exp_eight_schools_c$cov_err),
          exp_glm$mean_err, c(exp_glm$cor_err), c(exp_glm$cov_err),
          exp_student$mean_err,
          c(exp_student$cor_err),
          c(exp_student$cov_err),
          exp_crescent$mean_err,
          c(exp_crescent$cor_err),
          c(exp_crescent$cov_err),
          exp_mixture$mean_err,
          c(exp_mixture$cor_err),
          c(exp_mixture$cov_err),
          exp_disease_map$mean_err,
          c(exp_disease_map$cor_err),
          c(exp_disease_map$cov_err)
          ),
  model = as.factor(
    c(rep("SKIM", n_dim_vec[1] * (1 + 2 * n_dim_vec[1])),
      rep("8schools_nc", n_dim_vec[2] * (1 + 2 * n_dim_vec[2])),
      rep("8schools", n_dim_vec[3] * (1 + 2 * n_dim_vec[3])),
      rep("GLM", n_dim_vec[4] * (1 + 2 * n_dim_vec[4])),
      rep("student", n_dim_vec[5] * (1 + 2 * n_dim_vec[5])),
      rep("crescent", n_dim_vec[6] * (1 + 2 * n_dim_vec[6])),
      rep("mixture", n_dim_vec[7] * (1 + 2 * n_dim_vec[7])),
      rep("disease_map", n_dim_vec[8] * (1 + 2 * n_dim_vec[8]))
    )),
  quant = c(
    rep("mean", n_dim_vec[1]), rep("correlation", n_dim_vec[1]^2), rep("covariance", n_dim_vec[1]^2),
    rep("mean", n_dim_vec[2]), rep("correlation", n_dim_vec[2]^2), rep("covariance", n_dim_vec[2]^2),
    rep("mean", n_dim_vec[3]), rep("correlation", n_dim_vec[3]^2), rep("covariance", n_dim_vec[3]^2),
    rep("mean", n_dim_vec[4]), rep("correlation", n_dim_vec[4]^2), rep("covariance", n_dim_vec[4]^2),
    rep("mean", n_dim_vec[5]), rep("correlation", n_dim_vec[5]^2), rep("covariance", n_dim_vec[5]^2),
    rep("mean", n_dim_vec[6]), rep("correlation", n_dim_vec[6]^2), rep("covariance", n_dim_vec[6]^2),
    rep("mean", n_dim_vec[7]), rep("correlation", n_dim_vec[7]^2), rep("covariance", n_dim_vec[7]^2),
    rep("mean", n_dim_vec[8]), rep("correlation", n_dim_vec[8]^2), rep("covariance", n_dim_vec[8]^2)
    ),
  check_sym = c(
    rep(exp_skim$check_sym, each = n_dim_vec[1] * (1 + 2 * n_dim_vec[1])),
    rep(exp_eight_schools_nc$check_sym, each = n_dim_vec[2] * (1 + 2 * n_dim_vec[2])),
    rep(exp_eight_schools_c$check_sym, each = n_dim_vec[3] * (1 + 2 * n_dim_vec[3])),
    rep(exp_glm$check_sym, each = n_dim_vec[4] * (1 + 2 * n_dim_vec[4])),
    rep(1e-4, each = n_dim_vec[5] * (1 + 2 * n_dim_vec[5])),
    rep(exp_crescent$check_sym, each = n_dim_vec[6] * (1 + 2 * n_dim_vec[6])),
    rep(exp_mixture$check_sym, each = n_dim_vec[7] * (1 + 2 * n_dim_vec[7])),
    rep(exp_disease_map$check_sym, each = n_dim_vec[8] * (1 + 2 * n_dim_vec[8]))
  )
  )

# For student target, use analytical calculation of symmetry, rather than MC estimate.

plot.data$model <- factor(plot.data$model,
      levels = model_names[order(check_sym_vec, decreasing = FALSE)])


p <- ggplot(data = plot.data[plot.data$quant %in% c("correlation", "covariance"), ],
            aes(x = model, y = err, color = quant)) +
  geom_boxplot(outlier.alpha = 0.1) + theme_bw() +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 25),
                          text = element_text(size = 20)) +
  ylab("Scaled error") + xlab("") +
  theme(legend.title=element_blank())
p

plot.data$model <- factor(plot.data$model,
    levels = model_names[order(check_sym_vec, decreasing = TRUE)])

p <- ggplot(data = plot.data[plot.data$quant == "mean", ],
            aes(x = check_sym, y = err, color = model)) +
  geom_boxplot(outlier.alpha = 0.1, width = 0.25) + theme_bw() + 
  scale_x_log10() + scale_y_log10() + theme(legend.title=element_blank()) +
  theme(text = element_text(size = 20)) +
  xlab("Violation of symmetry") + ylab("Scaled error")
p


# Compute summaries for table
quant_oi <- "covariance"  # options: "mean", "correlation", "covariance"
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "student", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "disease_map", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "GLM", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "8schools_nc", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "mixture", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "SKIM", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "8schools", ]$err))
mean(abs(plot.data[plot.data$quant == quant_oi & plot.data$model == "crescent", ]$err))
