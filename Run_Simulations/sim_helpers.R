library(MASS)
library(tidyverse)
library(reshape2)
library(scoringRules)

# Simulate data from the UGLT sparse factor model.
# Set plt_structure = TRUE to fix pivots at l_j = j (PLT arrangement).
sim_data_UGLT <- function(N, V, q, plt_structure = FALSE) {
  sigma2 <- rep(1, V)
  Sigma <- diag(sigma2)
  theta <- rep(1, q)
  tau <- seq_len(q) / (2 * q + 1)
  pivots <- sort(sample(seq_len(V - 3), q))
  if (plt_structure) pivots <- seq_len(q)

  delta <- matrix(0L, V, q)
  for (j in seq_len(q)) {
    delta[pivots[j], j] <- 1L
    below <- seq.int(pivots[j] + 1L, V)
    if (length(below)) {
      delta[below, j] <- rbinom(length(below), 1, 1 - tau[j])
    }
  }

  Lambda <- matrix(0, V, q)
  for (i in seq_len(V)) {
    active <- which(delta[i, ] == 1L)
    if (length(active)) {
      Lambda[i, active] <- MASS::mvrnorm(
        1, rep(0, length(active)),
        diag(theta[active], length(active)) * sigma2[i]
      )
    }
  }
  Lambda <- Lambda %*% diag(diag(sign(Lambda[pivots, ])), q)
  factors <- t(MASS::mvrnorm(N, rep(0, q), diag(q)))
  data <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), Sigma))

  list(
    data = data, factors = factors, Lambda = Lambda, Sigma = Sigma,
    N = N, V = V, delta = delta, pivots = pivots, theta = theta
  )
}

# Build the argument list for run_mcmc_UGLT / run_mcmc_sparse_PLT.
# theta.shape and theta.rate are scalars (hyperparameters for G^{-1} prior on theta_j).
# q_overfit adds extra factors beyond the true q (for overfitting experiments).
get_starting_vals_uglt <- function(samp, n_runs = 7000, thin = 2, burn = 1000,
                                   q_overfit = 0, ident = TRUE) {
  V <- nrow(samp$data)
  list(
    N = samp$N,
    q = nrow(samp$factors) + q_overfit,
    n_runs = n_runs,
    alpha = rep(1.5, V),
    beta = rep(1.5, V),
    theta.shape = 1.5,
    theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    data = samp$data,
    thin = thin,
    ident = ident,
    burn = burn
  )
}

# Average CRPS for Lambda over all (i, j) entries.
# lambda_samps: list of n draws, each a v x q matrix.
avg_crps <- function(lambda_samps, true_val) {
  V <- nrow(true_val)
  q <- ncol(true_val)
  n <- length(lambda_samps)
  samps <- array(unlist(lambda_samps), dim = c(V, q, n))
  total <- 0
  for (i in seq_len(V)) {
    for (j in seq_len(q)) {
      total <- total + scoringRules::crps_sample(true_val[i, j], dat = samps[i, j, ])
    }
  }
  total / (V * q)
}

# Average MSE between data simulated from the true and estimated models.
# ... absorbs extra pmap columns not used by this function (e.g. V, q).
compute_goodness_of_fit <- function(N, Lambda, Sigma, factors,
                                    Lambda_est, sigma_est, factors_est, ...) {
  mse_vec <- replicate(20, {
    samp <- sim_data_3(N, Lambda, Sigma, factors = factors)
    samp2 <- sim_data_3(N, Lambda_est, diag(sigma_est), factors = factors_est)
    mean((samp - samp2)^2)
  })
  list(mean(mse_vec), mse_vec)
}

# Fraction of total variance in y explained by L L'.
contributed_variance_noise <- function(L, y) {
  sum(diag(L %*% t(L))) / sum(diag(cov(t(y))))
}

make_heat_map <- function(Lambda) {
  Lambda_df <- reshape2::melt(Lambda)
  Lambda_df$Var1 <- factor(Lambda_df$Var1, levels = rev(unique(Lambda_df$Var1)))
  Lambda_df$Var2 <- factor(Lambda_df$Var2)
  ggplot(Lambda_df, aes(y = Var1, x = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(x = "j", y = "i") +
    theme(axis.text.y = element_text(size = 7), legend.position = "none")
}

make_trace_acf_plot <- function(post_draws) {
  T_stat_df <- tibble(
    Value     = post_draws$T_stat,
    Iteration = seq_along(post_draws$T_stat)
  )
  list(
    acf = acf(post_draws$T_stat),
    trace = ggplot(T_stat_df, aes(x = Iteration, y = Value)) +
      geom_line(linewidth = 0.4, alpha = 0.5) +
      theme_minimal() +
      labs(x = "Iteration", y = "Value")
  )
}
