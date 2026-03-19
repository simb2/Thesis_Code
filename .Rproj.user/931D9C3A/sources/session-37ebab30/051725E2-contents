# finding posterior estimates:

post_mode_delta <- function(x) {
  delta_test <- array(NA, dim = c(length(x$delta_test), 
                                  nrow(x$delta_test[[1]]), # V
                                  ncol(x$delta_test[[1]])))
  
  for (i in seq_along(x$delta_test)) {
    delta_test[i, , ] <- x$delta_test[[i]]
  }
  
  counted_sparsities <- as_tibble(delta_test) |>
    group_by_all() |>
    summarise(n = n())
  post_mode_delta <- as.numeric(
    counted_sparsities[which.max(counted_sparsities$n), ]
  )
  matrix(post_mode_delta[1:(nrow(x$delta_test[[1]]) * ncol(x$delta_test[[1]]))], 
         nrow = nrow(x$delta_test[[1]]))
}

post_est_lambda <- function(x, post_mode_delta) {
  
  filtered_lambda <- list()
  true_vals <- c()
  post_mean_lambda <- matrix(0, nrow = nrow(x$Lambda_test[[1]]), ncol = ncol(x$Lambda_test[[1]]))
  for (n in seq_along(x$Lambda_test)) {
    if (all((x$Lambda_test[[n]] != 0) * 1 == post_mode_delta)) {
      print(TRUE)
      filtered_lambda[[n]] <- x$Lambda_test[[n]]
      true_vals[length(true_vals) + 1] <- n
    }
  }
  if (length(true_vals) > 1) {
    filtered_lambda <- filtered_lambda[true_vals]
    for (m in filtered_lambda) {
      post_mean_lambda = post_mean_lambda + m
    }
    post_mean_lambda = post_mean_lambda/(length(filtered_lambda))
  } else {
    post_mean_lambda <- filtered_lambda[[1]]
  }
  
  return(post_mean_lambda)
}

compute_post_means <- function(x) {
  sigma2_mean <- numeric(length = length(x$sigma_test[[1]]))
  tau_mean <- theta_mean <- numeric(length = length(x$tau_test[[1]]))
  factors_est <- matrix(0, nrow = nrow(x$W[[1]]), ncol = ncol(x$W[[1]]))
  for (i in seq_along(x$Lambda_test)) {
    sigma2_mean <- sigma2_mean + x$sigma_test[[i]]
    tau_mean <- tau_mean + x$tau_test[[i]]
    theta_mean <- theta_mean + x$theta_test[[i]]
    factors_est <- factors_est + x$W[[i]]
  }
  sigma2_mean <- (sigma2_mean)/length(x$Lambda_test)
  tau_mean <- (tau_mean)/length(x$Lambda_test)
  theta_mean <- (theta_mean)/length(x$Lambda_test)
  factors_est <- (factors_est)/length(x$Lambda_test)
  
  
  list(sigma2_mean = sigma2_mean, tau_mean = tau_mean, theta_mean = theta_mean, 
       factors_est = factors_est)
}
