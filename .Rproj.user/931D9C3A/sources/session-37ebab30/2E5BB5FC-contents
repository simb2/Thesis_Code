library(tidyverse)
library(reshape2)
library(MASS)
library(scoringRules) # CRPS
library(here)

source(here("MCMC_Algorithms", "run_mcmc_UGLT.R"))
source(here("MCMC_Algorithms" ,"sim_data_3.R"))

# Prior specifications and Globals ----------------------------------------
set.seed(8)
n_subj <- c(2000)
n_vars <- c(30)
n_factors <- 10

# Helper Functions --------------------------------------------------------

# First a function to simulate the data

sim_data_UGLT <- function(N, V, q) {
  # priors for theta and sigma respecively
  # True value for sigma
  sigma2 <- rep(1, V)
  Sigma <- diag(sigma2)
  theta <- rep(1, q) # L2 shrinkage
  tau <- (1:q) / (2*q + 1) # L1 shrinkage / spike probability
  hyperparams <- list(aH = 2, bH = 2) # priors for tau
  pivots <- sample(1:(V-3), q)
  pivots <- pivots[order(pivots)]
  pivots <- 1:q # for simulating a plt structure. 
  delta <- matrix(0, nrow = V, ncol = q)
  Phi <- diag(q)
  for (j in 1:q) {
    for (i in pivots[j]:V) {
      if (i == pivots[j]) {
        delta[i, j] <- 1
      } else {
        if (runif(1) > tau[j]) {
          delta[i, j] <- 1
        }
      }
    }
  }
  # now sample Lambda according to the prior.
  
  Lambda <- matrix(0, nrow = V, ncol = q)
  # First computing
  for (i in 1:V) {
    # First we make sure the row is non zero
    if (sum(delta[i, ]) != 0) {
      filter <- which(delta[i, ] == 1)
      theta_a <- theta[filter]
      L_0 <- diag(length(theta_a)) * theta_a
      Lambda[i, filter] <- mvrnorm(1, mu = rep(0, sum(delta[i, ])), Sigma = L_0 * sigma2[i])
    }
  }
  Lambda <- Lambda%*%diag(diag(sign(Lambda[pivots, ])))
  data <- matrix(data = NA, nrow = dim(Lambda)[1], ncol = N)
  factors <- matrix(data = NA, nrow = dim(Lambda)[2], ncol = N)
  for (i in 1:N) {
    f <- MASS::mvrnorm(1, mu = rep(0, dim(Lambda)[2]), Sigma = Phi)
    factors[, i] <- f
    data[, i] <- MASS::mvrnorm(1, mu = Lambda %*% f, Sigma = Sigma)
  }
  return(list(
    data = data, factors = factors, Lambda = Lambda, Sigma = Sigma, N = N, V = V,
    delta = delta, pivots = pivots, theta = theta)
  )
}



avg_crps <- function(lambda_samps1, true_val) {
  # compute the avg crps
  
  V <- nrow(lambda_samps1[1, , ])
  q <- ncol(lambda_samps1[1, , ])
  
  crps <- 0
  for (i in 1:V) {
    for (j in 1:q) {
      crps <- crps + crps_sample(true_val[i, j], dat = lambda_samps1[, i, j])
    }
  }
  return(crps / (V * q))
}

# Now we make a helper function for getting the starting values to actually run the sampler.

get_starting_vals <- function(samp) {
  
  data <- samp$data
  V <- samp$V
  
  return(list(
    N = samp$N, q = nrow(samp$factors) + 3, n_runs = 6000, alpha = rep(1.5, V), beta = rep(1.5, V), 
    theta.shape = rep(1.5, V), theta.rate = rep(1.5, V), hyperparams = list(aH = 2, bH = 2), data = data, thin = 2, burn = 1000
  ))
}

compute_goodness_of_fit <- function(N, Lambda, Sigma, factors, Lambda_est, sigma_est, factors_est) {
  runs <- 20
  mse_vec <- numeric(runs)
  
  for (i in 1:runs) {
    samp <- sim_data_3(N, Lambda, Sigma, factors = factors)
    samp2 <- sim_data_3(N, Lambda_est, Sigma = diag(sigma_est), factors = factors_est)
    mse_vec[i] <- mean((samp - samp2)^2)
  }
  return(list(mean(mse_vec), mse_vec))
}


compute_posterior_estimates <- function(post_samp) {
  Lambda_estimates <- apply(post_samp$Lambda, c(2, 3), mean)
  sigma_est <- apply(post_samp$Variances, 2, mean)
  factors_est <- apply(post_samp$Factors, c(2, 3), mean)
  return(list(Lambda_estimates = Lambda_estimates, sigma_est = sigma_est, factors_est = factors_est))
}

make_heat_map <- function(Lambda) {
  Lambda_df <- reshape2::melt(Lambda)
  Lambda_df$Var1 <- factor(Lambda_df$Var1, levels = rev(unique(Lambda_df$Var1)))
  Lambda_df$Var2 <- factor(Lambda_df$Var2)
  plt <- ggplot(Lambda_df, aes(y = Var1, x = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(x = "j", y = "i") +
    theme(axis.text.y = element_text(size = 7), legend.position = "none")
  return(plt)
}

make_trace_acf_plot <- function(post_samp) {
  T_stat_df <- tibble(
    Value = post_samp$T_stat,
    Iteration = 1:length(post_samp$T_stat)
  )
  acf <- acf(post_samp$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  return(list(acf = acf, trace = trace))
}


##############################################################

# Here is where I actually run the simulation

settings <- tidyr::crossing(n_subj, n_vars, n_factors)
names(settings) <- c("N", "V", "q")
samples <- purrr::pmap(settings, sim_data_UGLT)
starting_vals <- purrr::map(samples, get_starting_vals)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_UGLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_ests <- purrr::map(post_samps, "estimates")
trace_acf <- purrr::map(post_draws, make_trace_acf_plot)


source(here("MCMC_Algorithms", "run_mcmc_sparse_PLT.R"))

post_samps_plt <- purrr::map(starting_vals, ~ do.call(run_mcmc_sparse_PLT, .x))
post_draws_plt <- purrr::map(post_samps_plt, "draws")
post_ests_plt <- purrr::map(post_samps_plt, "estimates")
trace_acf <- purrr::map(post_draws_plt, make_trace_acf_plot)

# Making trace and acf plots:

gof_inputs <- tibble::tibble(
  Lambda = purrr::map(samples, "Lambda"),
  Sigma = purrr::map(samples, "Sigma"),
  factors = purrr::map(samples, "factors"),
  Lambda_est = purrr::map(post_ests[[1]], "lambda_est"),
  sigma_est = purrr::map(post_ests[[1]], "sigma2_mean"),
  factors_est = purrr::map(post_ests[[1]], "factors_est")
)


# making an acf plot

gof_inputs$N <- 200
compute_goodness_of_fit2 <- function(lambda_est, sigma2_mean, factors_est) {
  runs <- 20
  mse_vec <- numeric(runs)
  N <- 200
  Lambda <- purrr::map(samples, "Lambda")[[1]]
  Sigma <- purrr::map(samples, "Sigma")[[1]]
  factors <- purrr::map(samples, "factors")[[1]]
  for (i in 1:runs) {
    samp <- sim_data_3(N = N, Lambda = Lambda, Sigma = Sigma, factors = factors)
    samp2 <- sim_data_3(N, lambda_est, Sigma = diag(sigma2_mean), factors = factors_est)
    mse_vec[i] <- mean((samp - samp2)^2)
  }
  sd_mse <- sqrt(var(mse_vec))
  return(list(mean(mse_vec), mse_vec))
}

mse <- numeric(length(post_ests[[1]]))
saveRDS(post_draws, "overfit_plt_setting_10_post_draws_UGLT.rds")
saveRDS(post_draws_plt, "overfit_plt_setting_10_post_draws_plt_sparse.rds")
saveRDS(samples, "overfit_plt_setting_10_init_val.rds")



