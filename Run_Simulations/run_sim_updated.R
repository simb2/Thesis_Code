library(tidyverse)
library(reshape2)
library(MASS)
library(scoringRules) # CRPS
library(here)

source(here("Run_Simulations" ,"sim_data_3.R"))
source(here("Run_Simulations" , 'sim_study_noisy_variables_helpers.R'))
source(here("MCMC_Algorithms" , "run_mcmc_UGLT_Static.R"))

# Prior specifications and Globals ----------------------------------------
set.seed(20)
n_subj <- c(1200)
n_vars <- c(20, 30, 40)
n_factors <- 3

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
  
  V <- nrow(lambda_samps1[[1]])
  q <- ncol(lambda_samps1[[1]])
  lambda_samps <- array(data = NA, dim = c(length(lambda_samps1), V, q))
  for (i in seq_along(lambda_samps1)) {
    lambda_samps[i, , ] <- lambda_samps1[[i]]
  }
  
  crps <- 0
  for (i in 1:V) {
    for (j in 1:q) {
      crps <- crps + crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
    }
  }
  return(crps / (V * q))
}
# Now we make a helper function for getting the starting values to actually run the sampler.

get_starting_vals <- function(samp) {
  
  data <- samp$data
  V <- samp$V
  
  return(list(
    N = samp$N, q = nrow(samp$factors), n_runs = 7000, alpha = rep(1.5, V), beta = rep(1.5, V), 
    theta.shape = rep(1.5, V), theta.rate = rep(1.5, V), hyperparams = list(aH = 2, bH = 2), data = data, thin = 2, ident = TRUE,
    burn = 1000
  ))
  
}


get_starting_vals_plt <- function(data) {
  pc <- princomp(t(data))
  factor_est <- t(pc$scores[, 1:n_factors])
  Lambda_start <- pc$loadings[, 1:n_factors]
  sigma2_start <- diag(cov(t(data - Lambda_start %*% factor_est)))
  s2 <- 1
  nu <- 2
  V <- dim(data)[1]
  return(list(
    Lambda_0 = Lambda_start, sigma2_0 = sigma2_start, n_runs = 7000, alpha = rep(1.5, V), beta = rep(1.5, V), theta.shape = rep(1.5, V), theta.rate = rep(1.5, V), hyperparams = list(aH = 2, bH = 2),
    data = data, q = n_factors, nu = nu, s2 = s2, c_0 = 1, thin = 2, burn = 1000
  ))
}

compute_goodness_of_fit <- function(N, V, q, Lambda, Sigma, factors, Lambda_est, sigma_est, factors_est) {
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

make_trace_acf_plot_UGLT <- function(post_draws) {
  T_stat_df <- tibble(
    Value = post_draws$T_stat,
    Iteration = 1:length(post_draws$T_stat)
  )
  acf <- acf(post_draws$T_stat)
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
samples <- pmap(settings, sim_data_UGLT)

starting_vals <- purrr::map(samples, get_starting_vals)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_UGLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_estimates <- purrr::map(post_samps, "estimates")
# Creating the overall data frame to save everything in
contr_noise_inputs <- tibble(
  L = map(post_estimates, function(x){
    x[[1]]$lambda_est
  }),
  y = map(samples, "data")
)

noise <- pmap(contr_noise_inputs, contributed_variance_noise)

make_trace_acf_plot <- function(p_draws) {
  
  T_stat_df <- tibble(
    Value = p_draws$T_stat,
    Iteration = 1:length(p_draws$T_stat)
  )
  acf <- acf(p_draws$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  return(list(acf = acf, trace = trace))
}

trace_acf <- purrr::map(post_draws, make_trace_acf_plot)
make_trace_acf_plot(post_draws[[1]])
# Computing CRPS 
estimates <- list()
for (i in 1:3) {
  estimates[[i]] <- post_estimates[[i]][[1]]
}




gof_inputs <- settings |> mutate(
  Lambda = map(samples, "Lambda"),
  Sigma = map(samples, "Sigma"),
  factors = map(samples, "factors"),
  Lambda_est = map(estimates, "lambda_est"),
  sigma_est = map(estimates, "sigma2_mean"),
  factors_est = map(estimates, "factors_est")
)

mse_vals <- pmap(gof_inputs, compute_goodness_of_fit)
mse_vals

gof_inputs$mse <- mse_vals

crps_inputs <- tibble(
  true_val = gof_inputs$Lambda,
  lambda_samps1 = map(post_draws, "Lambda_test")
)

avg_crps_list <- function(lambda_samps1, true_val) {
  # compute the avg crps
  
  V <- nrow(lambda_samps1[[1]])
  q <- ncol(lambda_samps1[[1]])
  lambda_samps <- array(data = NA, dim = c(length(lambda_samps1), V, q))
  for (i in seq_along(lambda_samps1)) {
    lambda_samps[i, , ] <- lambda_samps1[[i]]
  }
  
  crps <- 0
  for (i in 1:V) {
    for (j in 1:q) {
      crps <- crps + crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
    }
  }
  return(crps / (V * q))
}
# Creating the overall data frame to save everything in


crps <- pmap(crps_inputs, avg_crps_list)




gof_inputs$crps <- unlist(crps)
# Making plots

est_heat_map <- map(gof_inputs$Lambda_est, make_heat_map)
true_heat_map <- map(gof_inputs$Lambda, make_heat_map)
# trace_acf <- map(post_samps, make_trace_acf_plot)

# Creating the overall data frame to save everything in

tru_contr_noise_inputs <- tibble(
  L = map(samples, "Lambda"), 
  y = map(samples, "data")
)

tru_noise <- pmap(tru_contr_noise_inputs, contributed_variance_noise)

sim_results <- gof_inputs |> mutate(
  est_heat_map = est_heat_map,
  true_heat_map = true_heat_map,
  trace_acf = trace_acf,
  post_samps = post_samps,
  noise = noise, 
  tru_noise = tru_noise
)

gof_inputs$noise = noise

saveRDS(sim_results, "results_uglt_setting_UGLT.rds")
saveRDS(gof_inputs, "analysis_uglt_setting_UGLT.rds")

# Making trace and acf plots:

###########################################################
###########################################################
###########################################################

source(here("Run_Simulations", "sim_data_2.R"))
source(here("Run_Simulations", "sim_data_3.R"))
source(here("MCMC_Algorithms", 'run_mcmc_sparse_PLT_Static.R'))


settings <- tidyr::crossing(n_subj, n_vars, n_factors)
names(settings) <- c("N", "V", "q")
settings_sparse_plt <- settings
samples_sparse_plt <- samples
starting_vals <- purrr::map(samples_sparse_plt, get_starting_vals)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_sparse_PLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_estimates <- purrr::map(post_samps, "estimates")

make_trace_acf_plot <- function(post_draws) {
  T_stat_df <- tibble(
    Value = post_draws$T_stat,
    Iteration = 1:length(post_draws$T_stat)
  )
  acf <- acf(post_draws$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  return(list(acf = acf, trace = trace))
}

trace_acf <- map(post_draws, make_trace_acf_plot)
trace_acf
# Computing CRPS 
estimates <- list()
for (i in 1:3) {
  estimates[[i]] <- post_estimates[[i]][[1]]
}




gof_inputs <- settings_sparse_plt |> mutate(
  Lambda = map(samples_sparse_plt, "Lambda"),
  Sigma = map(samples_sparse_plt, "Sigma"),
  factors = map(samples_sparse_plt, "factors"),
  Lambda_est = map(estimates, "lambda_est"),
  sigma_est = map(estimates, "sigma2_mean"),
  factors_est = map(estimates, "factors_est")
)


mse_vals <- pmap(gof_inputs, compute_goodness_of_fit)
mse_vals

gof_inputs$mse <- mse_vals

crps_inputs <- tibble(
  true_val = gof_inputs$Lambda,
  lambda_samps1 = map(post_draws, "Lambda_test")
)

avg_crps_list <- function(lambda_samps1, true_val) {
  # compute the avg crps
  
  V <- nrow(lambda_samps1[[1]])
  q <- ncol(lambda_samps1[[1]])
  lambda_samps <- array(data = NA, dim = c(length(lambda_samps1), V, q))
  for (i in seq_along(lambda_samps1)) {
    lambda_samps[i, , ] <- lambda_samps1[[i]]
  }
  
  crps <- 0
  for (i in 1:V) {
    for (j in 1:q) {
      crps <- crps + crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
    }
  }
  return(crps / (V * q))
}
# Creating the overall data frame to save everything in


crps <- pmap(crps_inputs, avg_crps_list)


contr_noise_inputs <- tibble(
  L = map(estimates, "lambda_est"),
  y = map(samples_sparse_plt, "data")
)
noise <- pmap(contr_noise_inputs, contributed_variance_noise)

gof_inputs$crps <- unlist(crps)

gof_inputs <- gof_inputs %>% mutate(
  noise = noise
)
# Making plots

est_heat_map <- map(gof_inputs$Lambda_est, make_heat_map)
true_heat_map <- map(gof_inputs$Lambda, make_heat_map)
# trace_acf <- map(post_samps, make_trace_acf_plot)

# Creating the overall data frame to save everything in

sim_results <- gof_inputs |> mutate(
  post_samps = post_samps
)

saveRDS(sim_results, "results_uglt_setting_plt.rds")
saveRDS(gof_inputs, "analysis__uglt_setting_plt.rds")