library(MASS)
sample_loadings_variances <- function(data, Factors, delta, theta, alpha, beta) {
  y <- data
  V <- dim(data)[1]  # number of variables
  q <- dim(Factors)[1]  # number of factors
  N <- dim(data)[2]  # number of observations
  factors <- t(Factors)
  
  Lambda_new <- matrix(0, V, q)
  sigma2_new <- numeric(V)
  # First computing 
  for (i in 1:V) {
    # First we make sure the row is non zero
    if (sum(delta[i, ]) != 0) {
      filter <- which(delta[i, ] == 1)
      theta_a <- theta[filter]
      L_0 <- diag(length(theta_a))*theta_a
      L_0_inv <- solve(L_0)
      X_i_delta <- factors[, filter]
      L_iN <- solve(t(X_i_delta)%*%X_i_delta + L_0_inv)
      m_iN <- t(X_i_delta) %*% y[i, ]
      mi <- L_iN %*% m_iN
      sigma2_new[i] <- 1/rgamma(1, shape = alpha[i] + N/2, rate = 0.5*(sum(y[i, ]^2) - t(m_iN) %*% L_iN %*% m_iN))
      Lambda_new[i, filter] <- mvrnorm(1, mu = mi, Sigma = L_iN * sigma2_new[i])
    } else {
      sigma2_new[i] <- 1/rgamma(1, shape = alpha[i] + N/2, rate =  0.5*(sum(y[i, ]^2)))
    }
  }
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  return(list(Lambda_new = Lambda_new, sigma2_new = sigma2_new))
}