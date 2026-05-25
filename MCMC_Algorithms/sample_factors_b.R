#' Sample latent factors (PLT model)
#'
#' Draws each column \eqn{f_t \sim N(V_F \Lambda' D^{-1} y_t,\; V_F)} jointly for all t, where
#' \eqn{D = \mathrm{diag}(\sigma^2)} and \eqn{V_F = (I_q + \Lambda' D^{-1} \Lambda)^{-1}} is the posterior covariance.
#'
#' @param Lambda \eqn{v \times q} loading matrix \eqn{\Lambda}.
#' @param sigma2 Length-v vector of idiosyncratic variances \eqn{\sigma^2}.
#' @param data \eqn{v \times N} data matrix \eqn{Y}.
#' @param q Number of factors \eqn{q}.
#' @return \eqn{q \times N} matrix of factor draws \eqn{F}.
sample_factors <- function(Lambda, sigma2, data, q) {
  N <- ncol(data)
  Vf <- solve(diag(q) + crossprod(Lambda / sqrt(sigma2)))
  mean_mat <- Vf %*% crossprod(Lambda, data / sigma2)  # q x N
  factors <- mean_mat + t(chol(Vf)) %*% matrix(rnorm(q * N), q, N)
  return(factors)
}
