library(MASS)

#' Sample loadings and idiosyncratic variances (UGLT model)
#'
#' Draws \eqn{\lambda_{k,\delta_k}} and \eqn{\sigma^2_k} jointly from their Normal-Inverse-Gamma conditional
#' posteriors for each row k. Entries with \eqn{\delta_{kj} = 0} remain zero.
#' Prior: \eqn{\lambda_{kj} \mid \delta_{kj}=1 \sim N(0, \theta_j)} and \eqn{\sigma^2_k \sim G^{-1}(\alpha_k, \beta_k)}.
#'
#' @param data \eqn{v \times N} data matrix \eqn{Y}.
#' @param Factors \eqn{q \times N} factor matrix \eqn{F}.
#' @param delta \eqn{v \times q} binary sparsity matrix \eqn{\delta}.
#' @param theta Length-q column shrinkage vector \eqn{\theta}.
#' @param alpha Length-v shape parameters \eqn{\alpha_1, \ldots, \alpha_v} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param beta Length-v rate parameters \eqn{\beta_1, \ldots, \beta_v} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param inner_prod_y Length-v vector of precomputed row inner products \eqn{\sum_t y_{kt}^2}.
#' @return List with \code{Lambda_new} (\eqn{v \times q} matrix) and \code{sigma2_new} (length-v vector).
sample_loadings_variances <- function(data, Factors, delta, theta, alpha, beta, inner_prod_y) {
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
      L_0_inv <- diag(1/theta_a, nrow = length(theta_a))
      X_i_delta <- factors[, filter, drop = FALSE]
      P_i <- crossprod(X_i_delta) + L_0_inv
      m_iN <- drop(crossprod(X_i_delta, y[i, ]))
      L_chol <- t(chol(P_i))
      x <- forwardsolve(L_chol, m_iN)
      sigma2_new[i] <- 1/rgamma(1, shape = alpha[i] + N/2,
                                 rate = beta[i] + 0.5*(inner_prod_y[i] - sum(x^2)))
      mi <- backsolve(t(L_chol), x)
      Lambda_new[i, filter] <- mi + sqrt(sigma2_new[i]) * backsolve(t(L_chol), rnorm(length(filter)))
    } else {
      sigma2_new[i] <- 1/rgamma(1, shape = alpha[i] + N/2, rate = beta[i] + 0.5*inner_prod_y[i])
    }
  }
  return(list(Lambda_new = Lambda_new, sigma2_new = sigma2_new))
}