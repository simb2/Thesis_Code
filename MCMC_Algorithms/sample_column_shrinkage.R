#' Sample column shrinkage parameters
#'
#' Draws \eqn{\theta_j \sim G^{-1}(a_\theta + d_j/2,\; b_\theta + r_j/2)} for each j, where
#' \eqn{d_j = \sum_k \delta_{kj}} counts active loadings and
#' \eqn{r_j = \sum_k \delta_{kj} \lambda_{kj}^2 / \sigma^2_k} is the precision-weighted sum of squared loadings.
#'
#' @param theta.shape Scalar shape hyperparameter \eqn{a_\theta}.
#' @param theta.rate Scalar rate hyperparameter \eqn{b_\theta}.
#' @param Lambda \eqn{v \times q} loading matrix \eqn{\Lambda}.
#' @param sigma2 Length-v vector of idiosyncratic variances \eqn{\sigma^2_1, \ldots, \sigma^2_v}.
#' @param delta \eqn{v \times q} binary sparsity matrix \eqn{\delta}.
#' @return Length-q vector of draws \eqn{\theta_1, \ldots, \theta_q}.
sample_column_shrinkage <- function(theta.shape, theta.rate, Lambda, sigma2, delta) {
  column.sizes <- colSums(delta)
  rate.update <- colSums(delta * (Lambda^2 / sigma2))
  theta.new <- 1/rgamma(ncol(delta), theta.shape + 0.5*column.sizes, theta.rate + 0.5*rate.update)
  return(theta.new)
}