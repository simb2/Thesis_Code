#' Sample slab probabilities
#'
#' Draws \eqn{\tau_j \sim \mathrm{Beta}(a_\tau + d_j - 1,\; b_\tau + v - l_j - d_j + 1)} for each j,
#' where \eqn{d_j = \sum_k \delta_{kj}} counts active loadings and \eqn{l_j} is the pivot row index.
#'
#' @param hyperparams List with elements \code{aH} (\eqn{a_\tau}) and \code{bH} (\eqn{b_\tau}).
#' @param delta \eqn{v \times q} binary sparsity matrix \eqn{\delta}.
#' @param pivots Length-q integer vector of pivot row indices \eqn{l_1, \ldots, l_q}.
#' @return Length-q vector of draws \eqn{\tau_1, \ldots, \tau_q}.
sample_tau <- function(hyperparams, delta, pivots) {
  colsums.delta <- colSums(delta)
  new.tau <- numeric(ncol(delta))

  for (j in 1:dim(delta)[2]) {
    if (hyperparams$bH + nrow(delta) - pivots[j] -  colsums.delta[j] + 1 <  0) {
      browser()
    }
    
    new.tau[j] <- rbeta(1, hyperparams$aH + colsums.delta[j] - 1,
                        hyperparams$bH + nrow(delta) - pivots[j] -  colsums.delta[j] + 1)
    
    if (is.na(new.tau[j])) {
      browser()
    }
  }
  return(new.tau)
}