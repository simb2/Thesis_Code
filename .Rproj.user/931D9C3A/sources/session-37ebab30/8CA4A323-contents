#' Compute marginal posterior Odds:
#'
#' @param delta The current sparsity matrix (V x q)
#' @param y the data matrix (V x N)
#' @param factors a matrix containing the factors (q x N)
#' @param theta a vector containing the column wise L2 shrinkage (q x 1)
#' @param alpha a vector containing the shape parameter for the priors on the sigma's (V x 1)
#' @param beta a vector for the scale parameters for the prior on the sigma's (V x 1)

compute_log_likelihood_ratio <- function(delta, l_new, j, factors, y, alpha, beta, theta) {
  i <- l_new
  N <- dim(y)[2]
  
  # First we gotta compute the likelihood under delta_ij = 1
  delta[i, j] <- 1
  filter <- which(delta[i, ] == 1)
  theta_a <- theta[filter]
  L_0 <- diag(length(theta_a))*theta_a # need to double check this (in the original paper it depends on i but I can't see it)
  L_0_inv <- solve(L_0)
  factors <- t(factors)
  # Filtering out the columns
  X_i_delta <- factors[, filter]
  L_iN_inv <- t(X_i_delta)%*%X_i_delta + L_0_inv
  L_iN <- solve(t(X_i_delta)%*%X_i_delta + L_0_inv)
  M_i <- L_iN %*% t(X_i_delta) %*% y[i, ]
  
  log_lik <- 0.5*(log(det(L_iN)) - log(det(L_0))) - 0.5*N*log(2*pi) +
    alpha[i]*log(beta[i]) - lgamma(alpha[i])+ lgamma(N/2 + alpha[i]) -
    (N/2 + alpha[i])*log(beta[i] * 0.5*(t(y[i, ])  %*% (y[i, ]) -
                                          t(M_i) %*% solve(L_iN) %*% M_i))
  
  delta[i, j] <- 0
  
  if (sum(delta[i, ]) == 0) {
    log_lik_null <- -N/2 * log(2*pi) + alpha[i]*log(beta[i]) - lgamma(alpha[i]) +
      lgamma(N/2 + alpha[i]) - (N/2 + alpha[i])*log(beta[i] + 0.5*t(y[i, ]) %*% y[i, ] )
  } else {
    delta[i, j] <- 0
    filter <- which(delta[i, ] == 1)
    theta_a <- theta[filter]
    L_0 <- diag(length(theta_a))*theta_a
    L_0_inv <- solve(L_0)
    
    # Filtering out the columns
    X_i_delta <- factors[, filter]
    L_iN <- solve(t(X_i_delta)%*%X_i_delta + L_0_inv)
    
    
    M_i <- L_iN %*% t(X_i_delta) %*% y[i, ]
    log_lik_null <- 0.5*(log(det(L_iN)) - log(det(L_0))) - 0.5*N*log(2*pi) +
      alpha[i]*log(beta[i]) - lgamma(alpha[i])+ lgamma(N/2 + alpha[i]) -
      (N/2 + alpha[i])*log(beta[i] * 0.5*(t(y[i, ])  %*% y[i, ] -
                                            t(M_i) %*% solve(L_iN) %*% M_i))
  }
  PO = log_lik - log_lik_null
  return(PO)
}
