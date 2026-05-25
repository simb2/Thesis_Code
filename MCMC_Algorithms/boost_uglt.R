#' Parameter expansion boost step (UGLT)
#'
#' Improves mixing via GIG-based parameter expansion. Samples auxiliary scales
#' \eqn{\psi_j \sim \mathrm{GIG}(s_j/2 - 3/2 - N/2,\; \chi_j,\; \eta_j)} for each column j,
#' where \eqn{s_j} is the number of non-zero loadings in column j, then
#' back-transforms \eqn{(\Lambda, F)} to an equivalent parameterisation.
#'
#' @param Lambda \eqn{v \times q} loading matrix \eqn{\Lambda}.
#' @param factors \eqn{q \times N} factor matrix \eqn{F}.
#' @param sigma2 Length-v vector of idiosyncratic variances \eqn{\sigma^2}.
#' @param theta Length-q column shrinkage vector \eqn{\theta}.
#' @return List with \code{Lambda_new} (\eqn{v \times q}) and \code{factors_new} (\eqn{q \times N}).
boost_uglt <- function(Lambda, factors, sigma2, theta) {
  # First we sample phi:
  q <- ncol(Lambda)
  V <- nrow(Lambda)
  N <- ncol(factors)
  psi <- diag(1 / rgamma(q, 1.5, rate = 1.5))
  
  factors_psi <- psi^0.5 %*% factors
  Lambda_psi <- Lambda %*% solve(psi)^(0.5)
  
  psi_new_vec <- numeric(q)
  
  
  s <- colSums(Lambda != 0)
  
  for (j in 1:q) {
    psi_new_vec[j] = GIGrvg::rgig(1, lambda = s[j]/2 - 1.5 - N/2,
                                  chi = 3 + sum(factors_psi[j, ]^2),
                                  psi = sum(Lambda_psi[, j]^2 / sigma2) / theta[j])
    if (is.na(psi_new_vec[j])) browser()
  }
  
  # un transofring:
  psi_new <- diag(psi_new_vec)
  factors_new <- solve(psi_new)^(0.5) %*% factors_psi
  Lambda_new <- Lambda_psi %*% psi_new^(0.5)
  return(list(Lambda_new = Lambda_new, factors_new = factors_new))
}