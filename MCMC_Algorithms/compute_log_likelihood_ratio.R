#' Compute the log marginal posterior odds for a single sparsity entry
#'
#' Returns \eqn{\log p(y_i \mid \delta_{ij}=1, \delta_{i,-j}, F, \theta, \alpha, \beta) - \log p(y_i \mid \delta_{ij}=0, \delta_{i,-j}, F, \theta, \alpha, \beta)},
#' integrating out \eqn{\lambda_i} and \eqn{\sigma^2_i} under the Normal-Inverse-Gamma prior.
#' Used within the pivot move samplers.
#'
#' @param delta \eqn{v \times q} binary sparsity matrix \eqn{\delta} (modified locally; not altered on return).
#' @param l_new Row index i of the entry being evaluated.
#' @param j Column index j of the entry being evaluated.
#' @param factors \eqn{q \times N} factor matrix \eqn{F}.
#' @param y \eqn{v \times N} data matrix \eqn{Y}.
#' @param alpha Length-v shape parameters \eqn{\alpha} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param beta Length-v rate parameters \eqn{\beta} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param theta Length-q column shrinkage vector \eqn{\theta}.
#' @return Scalar log posterior odds.
compute_log_likelihood_ratio <- function(delta, l_new, j, factors, y, alpha, beta, theta) {
  i <- l_new
  N <- dim(y)[2]
  
  # First we gotta compute the likelihood under delta_ij = 1
  delta[i, j] <- 1
  filter <- (delta[i, ] == 1)
  theta_a <- theta[filter]
  L_0 <- diag(length(theta_a))*theta_a # need to double check this (in the original paper it depends on i but I can't see it)
  L_0_inv <- solve(L_0) # this needs to be reviesd. 
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
