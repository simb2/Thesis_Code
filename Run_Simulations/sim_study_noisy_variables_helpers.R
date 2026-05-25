sim_data_L <- function(N, V, q) {
  Phi<- diag(q)
  Sigma <- diag(V)
  Lambda <- matrix(0, nrow = V, ncol = 2)
  for (i in seq_len(V)) {
    for (j in seq_len(min(i, 2))) {
      Lambda[i, j] <- if (i == j) abs(rnorm(1, sd = sqrt(0.5))) else rnorm(1, sd = sqrt(0.5))
    }
  }
  factors <- t(MASS::mvrnorm(N, rep(0, 2), Phi))
  data <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), Sigma))
  list(data = data, factors = factors, Lambda = Lambda, Sigma = Sigma, N = N, V = V)
}

emp_crps <- function(samp1, samp2, true_val) {
  N <- length(samp1); M <- length(samp2)
  term1 <- term2 <- 0
  for (i in seq_len(N)) {
    for (j in seq_len(M))
      term1 <- term1 + abs(samp1[i] - samp2[j])
    term2 <- term2 + abs(samp1[i] - true_val)
  }
  2 / (N * M) * term1 - 1 / N * term2
}

# Array version of avg_crps for the PLT model, where lambda_samps is a [draws x v x q] array.
avg_crps_array <- function(lambda_samps, true_val) {
  V <- nrow(lambda_samps[1, , ])
  q <- ncol(lambda_samps[1, , ])
  total <- 0
  for (i in seq_len(V))
    for (j in seq_len(q))
      total <- total + scoringRules::crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
  total / (V * q)
}

# Starting values for run_mcmc (non-sparse PLT model). Takes the raw data matrix.
get_starting_vals <- function(data) {
  pc <- princomp(t(data))
  factor_est <- t(pc$scores[, 1:2])
  Lambda_start <- pc$loadings[, 1:2]
  sigma2_start <- diag(cov(t(data - Lambda_start %*% factor_est)))
  list(
    Lambda_0 = Lambda_start, sigma2_0 = sigma2_start, n_runs = 5000,
    data = data, q = n_factors, nu = 3, s2 = 0.5, c_0 = 1, thin = 10, burn = 2000
  )
}

compute_posterior_estimates <- function(post_samp) {
  list(
    Lambda_estimates = apply(post_samp$Lambda, c(2, 3), mean),
    sigma_est = apply(post_samp$Variances, 2, mean),
    factors_est = apply(post_samp$Factors, c(2, 3), mean)
  )
}
