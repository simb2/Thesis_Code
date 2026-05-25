#' Remove spurious factors from posterior draws
#'
#' A factor column j is spurious when only its pivot row is active
#' (\eqn{\sum_k \delta_{kj} = 1}). Removes each spurious column from \eqn{\Lambda}, \eqn{F},
#' \eqn{\delta}, \eqn{\theta}, and \eqn{\tau} across all draws, and recomputes pivot indices accordingly.
#'
#' @param pivot_draws List of length-q pivot vectors, one per draw.
#' @param factor_draws List of \eqn{q \times N} factor matrices, one per draw.
#' @param lambda_draws List of \eqn{v \times q} loading matrices, one per draw.
#' @param delta_draws List of \eqn{v \times q} sparsity matrices, one per draw.
#' @param theta_draws List of length-q \eqn{\theta} vectors, one per draw.
#' @param sigma_draws List of length-v \eqn{\sigma^2} vectors, one per draw.
#' @param tau_draws List of length-q \eqn{\tau} vectors, one per draw.
#' @param n_draws Number of MCMC draws.
#' @return List with the same fields, spurious columns removed from each draw.
filter_factors <- function(pivot_draws, factor_draws, lambda_draws, delta_draws,
                           theta_draws, sigma_draws, tau_draws, n_draws) {
  # First we identify the spurious columns
  for (i in 1:n_draws) {
    spur_col <- colSums(delta_draws[[i]]) == 1
    if (length(spur_col) > 0) {
      # we merge the spurious column into the idiosyncratic variance:
      col_sp <- which(spur_col)
      k = 0
      for (col in col_sp) {
        # adding the variance component to 'sigma'
        # first we need to find the row
        col = col - k
        row <- which(delta_draws[[i]][, col] == 1)[1]
        sigma_draws[[i]][row] <- sigma_draws[[i]][row]
        factor_draws[[i]] <- factor_draws[[i]][-col, ]
        lambda_draws[[i]] <- lambda_draws[[i]][ , -col]
        delta_draws[[i]] <- delta_draws[[i]][, -col]
        if (is.null(dim(delta_draws[[i]]))) {
          pivot_draws[[i]] <- which(delta_draws[[i]] == 1)[1]
        } else {
          pivot_draws[[i]] <- apply(delta_draws[[i]], 2, function(j) which(j != 0)[1])
        }
        theta_draws[[i]] <- theta_draws[[i]][-col]
        tau_draws[[i]] <- tau_draws[[i]][-col]
        k = k + 1
      }
    }
  }
  return(list(sigma_draws = sigma_draws, factor_draws = factor_draws, 
              lambda_draws = lambda_draws, delta_draws = delta_draws,
              tau_draws = tau_draws, theta_draws = theta_draws, pivot_draws = pivot_draws))
}