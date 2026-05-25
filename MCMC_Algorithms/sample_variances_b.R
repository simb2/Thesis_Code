#' Sample idiosyncratic variances (PLT model)
#'
#' Draws \eqn{\sigma^2_k \sim G^{-1}((\nu + N)/2,\; (\nu s^2 + \mathrm{SSR}_k)/2)} for all k simultaneously,
#' where \eqn{\mathrm{SSR}_k = \sum_t (y_{kt} - \lambda_k' f_t)^2}.
#'
#' @param nu Degrees-of-freedom hyperparameter \eqn{\nu}.
#' @param s2 Scale hyperparameter \eqn{s^2}.
#' @param data \eqn{v \times N} data matrix \eqn{Y}.
#' @param factors \eqn{q \times N} factor matrix \eqn{F}.
#' @param Lambda \eqn{v \times q} loading matrix \eqn{\Lambda}.
#' @return Length-v vector of draws \eqn{\sigma^2_1, \ldots, \sigma^2_v}.
sample_variances <- function(nu, s2, data, factors, Lambda) {
  N <- dim(data)[2]
  residuals <- data - Lambda %*% factors  # V x N
  di <- rowSums(residuals^2)
  sigma2_new <- 1/rgamma(dim(data)[1], shape = (nu + N)/2, rate = (nu*s2 + di)/2)
  return(sigma2_new)
}
