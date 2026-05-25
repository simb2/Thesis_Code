#' Sample loadings (PLT model)
#'
#' Draws rows of \eqn{\Lambda} from their Normal conditional posteriors under the PLT
#' lower-triangular structure. Row h draws \eqn{\lambda_{h,1:h}} using factors 1:h only;
#' rows \eqn{h > q} draw all q entries. Prior: \eqn{\lambda_{kj} \sim N(0, c_0)}.
#'
#' @param data \eqn{v \times N} data matrix \eqn{Y}.
#' @param sigma2 Length-v vector of idiosyncratic variances \eqn{\sigma^2}.
#' @param Factors \eqn{q \times N} factor matrix \eqn{F}.
#' @param c_0 Scalar prior variance for loading entries.
#' @return \eqn{v \times q} loading matrix \eqn{\Lambda}.
sample_loadings <- function(data, sigma2, Factors, c_0) {
  V <- dim(data)[1]
  q <- dim(Factors)[1]
  N <- dim(data)[2]

  Lambda_new <- matrix(NA, V, q)

  # First q rows: lower-triangular structure, row h uses only factors 1:h
  for (h in 1:q) {
    F_h <- Factors[1:h, , drop = FALSE]  # h x N
    post_cov <- solve((1/c_0) * diag(h) + (1/sigma2[h]) * tcrossprod(F_h))
    post_mean <- post_cov %*% ((1/sigma2[h]) * F_h %*% data[h, ])
    sample_vals <- MASS::mvrnorm(1, post_mean, post_cov)
    Lambda_new[h, ] <- c(sample_vals, rep(0, q - h))
  }

  # Rows q+1 to V: all q factors active; tcrossprod(Factors) is the same for all
  FF <- tcrossprod(Factors)  # q x q, precomputed once
  for (h in (q+1):V) {
    post_cov <- solve((1/c_0) * diag(q) + (1/sigma2[h]) * FF)
    post_mean <- post_cov %*% ((1/sigma2[h]) * Factors %*% data[h, ])
    Lambda_new[h, ] <- MASS::mvrnorm(1, post_mean, post_cov)
  }

  return(Lambda_new)
}
