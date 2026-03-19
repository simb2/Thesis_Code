sample_variances <- function(nu, s2, data, factors, Lambda) {
  sigma2_new <- numeric(length = dim(data)[1])
  N <- dim(data)[2]
  
  for (i in 1:dim(data)[1]) {
    shape_param <- (nu + N) / 2
    
    di <- 0
    for (t in 1:N) {
      residual <- data[i, t] - t(Lambda[i, ]) %*% factors[, t]
      di <- di + residual^2
    }
    rate_param <- (nu * s2 + di) / 2
    
    # Sample from Inverse Gamma using 1/gama
    sigma2_new[i] <- 1 / rgamma(1, shape = shape_param, rate = rate_param)
  }
  return(sigma2_new)
}