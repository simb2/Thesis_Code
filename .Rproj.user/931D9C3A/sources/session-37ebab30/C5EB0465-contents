# Run MCMC UGLT
library(MASS)
library(cli)
library(sparvaride)
source(here("MCMC_Algorithms", "sample_column_shrinkage.R"))
source(here("MCMC_Algorithms", "sample_tau.R"))
source(here("MCMC_Algorithms", "sample_pivots.R"))
source(here("MCMC_Algorithms", "sample_factors_b.R"))
source(here("MCMC_Algorithms", "sample_sparsity.R"))
source(here("MCMC_Algorithms", "sample_loadings_delta.R"))
source(here("MCMC_Algorithms", "boost_uglt.R"))
source(here("MCMC_Algorithms", "filter_factors.R"))
source(here("MCMC_Algorithms", "compute_modes.R"))


run_mcmc_UGLT <- function(N, q, n_runs, alpha, beta, theta.shape, theta.rate, hyperparams, data, thin = 1, burn = 1,
                          ident = TRUE) {
  cli_progress_bar("Sampling from Posterior . . .", total = n_runs)
  y <- data
  V <- nrow(y)
  delta_start <- matrix(1, nrow = V, ncol = q)
  delta_start[upper.tri(delta_start, diag = FALSE)] <- 0
  pivots_start <- 1:q
  # sample size
  delta_test <- list()
  Lambda_test <- list()
  sigma_test <- list()
  tau_test <- list()
  theta_test <- list()
  W <- list()
  T_stat <- numeric(n_runs)
  pivot_test <- list()
  accepted <- 0
  delta_test[[1]] <- delta_start
  pivot_test[[1]] <- pivots_start
  tau_test[[1]] <- rep(0.5, q) # how to choose starting values for this?
  
  theta_test[[1]] <- rep(1, q) # how to choose starting values for this?
  pc <- princomp(t(y))
  scores <- pc$scores[, 1:q]
  W[[1]] <- t(scores)
  Lambda_est <- pc$loadings[, 1:q]
  sigma_test[[1]] <- diag(cov(t(y - Lambda_est %*% W[[1]])))
  
  
  for (i in 2:n_runs) {
    tau_test[[i]] <- sample_tau(hyperparams, delta_test[[i - 1]], pivot_test[[i - 1]])
    res <- sample_sparsity(y, W[[i - 1]], tau_test[[i]], theta_test[[i - 1]], delta_test[[i - 1]], alpha, beta)
    pivots_new <- apply(res$delta_new, 2, function(col) which(col != 0)[1])
    res2 <- update_pivots(
      delta = res$delta_new, theta = theta_test[[i - 1]], pivots = pivots_new, factors = W[[i - 1]], y = y,
      alpha = alpha, beta = beta, hyperparams = hyperparams, move_probs = list(
        pshift = 0.3, pswitch = 0.3,
        pa = 0.5
      )
    )
    delta_test[[i]] <- res2$delta
    pivot_test[[i]] <- res2$pivots
    pivot_test[[i]] <- unlist(pivot_test[[i]])
    res3 <- sample_loadings_variances(y, W[[i - 1]], delta_test[[i]], theta_test[[i - 1]], alpha, beta)
    Lambda_test[[i]] <- res3$Lambda_new
    sigma_test[[i]] <- res3$sigma2_new
    W[[i]] <- sample_factors(Lambda_test[[i]], sigma_test[[i]], y, q)
    theta_test[[i]] <- sample_column_shrinkage(
      theta.shape, theta.rate, Lambda_test[[i]],
      sigma_test[[i]], delta_test[[i]]
    )
    
    T_stat[[i]] <- sum(diag(0.5 * (Lambda_test[[i]] %*% t(Lambda_test[[i]]) + diag(sigma_test[[i]]))))
    boost <- boost_uglt(Lambda = Lambda_test[[i]], factors = W[[i]], sigma2 = sigma_test[[i]], theta = theta_test[[i]])
    Lambda_test[[i]] <- boost$Lambda_new
    W[[i]] <- boost$factors_new
    cli_progress_update()
  }
  
  thin_burn <- seq(burn, n_runs, by = thin)
  Lambda_test <- Lambda_test[thin_burn]
  sigma_test <- sigma_test[thin_burn]
  tau_test <- tau_test[thin_burn]
  theta_test <- theta_test[thin_burn]
  delta_test <- delta_test[thin_burn]
  pivot_test <- pivot_test[thin_burn]
  T_stat <- T_stat[thin_burn]
  W <- W[thin_burn]
  print(length(thin_burn))
  pivot_test <- as.matrix(pivot_test)
  
  if (ident) {
    delta_test_id <- lapply(delta_test, function(x) x[rowSums(x) > 0, ])
    
    identifiable <- unlist(lapply(delta_test_id, sparvaride::counting_rule_holds))
    Lambda_test <- Lambda_test[identifiable]
    sigma_test <- sigma_test[identifiable]
    tau_test <- tau_test[identifiable]
    theta_test <- theta_test[identifiable]
    delta_test <- delta_test[identifiable]
    pivot_test <- pivot_test[identifiable]
    T_stat <- T_stat[identifiable]
    W <- W[identifiable]
  }
  
  for (i in seq_along(Lambda_test)) {
    pivot_test[[i]] <- unlist(pivot_test[[i]])
    delta_test[[i]] <- delta_test[[i]][, order(pivot_test[[i]])]
    Lambda_test[[i]] <- Lambda_test[[i]][, order(pivot_test[[i]])]
    
    for (f in 1:dim(data)[2]) {
      W[[i]][, f] <- W[[i]][order(pivot_test[[i]]), f] %*% diag(diag(sign(Lambda_test[[i]][pivot_test[[i]][order(pivot_test[[i]])], ])))
    }
    Lambda_test[[i]] <- Lambda_test[[i]] %*% diag(diag(sign(Lambda_test[[i]][pivot_test[[i]][order(pivot_test[[i]])], ])))
    pivot_test[[i]] <- pivot_test[[i]][order(pivot_test[[i]])]
  }
  # Now we want to count the number of things with different sizes
  r <- numeric(length(Lambda_test))
  for (i in seq_along(Lambda_test)) {
    r[i] <- ncol(Lambda_test[[i]])
  }
  # first get everything into a data frame.
  
  draw_df <- tibble(
    pivot_test = pivot_test,
    W = W,
    Lambda_test = Lambda_test,
    delta_test = delta_test,
    theta_test = theta_test,
    tau_test = tau_test,
    sigma_test = sigma_test,
    index = seq_along(W),
    r = r,
    T_stat = T_stat
  )
  r_vals <- unique(r)
  
  estimate_results <- purrr::map(r_vals, function(val) {
    dfr <- draw_df |> filter(r == val)
    delta_mode <- post_mode_delta(dfr)
    append(
      list(
        delta_est = delta_mode,
        lambda_est = post_est_lambda(
          x = dfr,
          post_mode_delta = delta_mode
        )
      ),
      compute_post_means(x = dfr)
    )
  })
  estimate_results$r <- r_vals
  
  
  return(list(
    # Posterior summaries
    estimates = estimate_results,
    # Full MCMC draws
    draws = draw_df
  ))
}
