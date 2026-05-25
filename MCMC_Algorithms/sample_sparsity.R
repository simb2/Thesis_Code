#' Sample the sparsity matrix
#'
#' Draws \eqn{\delta} using the vectorized multi-move sampler of Frühwirth-Schnatter,
#' Hosszejni & Lopes (2025), Algorithm D.1/D.2. For each column j, rows below
#' the pivot \eqn{l_j} are split into zero rows (Step I-a, fully vectorized) and
#' non-zero rows (Step I-b, one Cholesky per row). Indicators are then updated
#' jointly via MH accept/reject (Step I-c).
#'
#' @param y \eqn{v \times N} data matrix \eqn{Y}.
#' @param factors \eqn{q \times N} factor matrix \eqn{F}.
#' @param tau Length-q slab probability vector \eqn{\tau_1, \ldots, \tau_q}.
#' @param theta Length-q column shrinkage vector \eqn{\theta_1, \ldots, \theta_q}.
#' @param delta \eqn{v \times q} binary sparsity matrix \eqn{\delta} (current state).
#' @param alpha Length-v shape parameters \eqn{\alpha_1, \ldots, \alpha_v} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param beta Length-v rate parameters \eqn{\beta_1, \ldots, \beta_v} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param inner_prod_y Length-v vector of precomputed row inner products \eqn{\sum_t y_{kt}^2}.
#' @return List with \code{delta_new} (\eqn{v \times q} updated sparsity matrix).

sample_sparsity <- function(y, factors, tau, theta, delta, alpha, beta, inner_prod_y) {
  V <- nrow(y)
  N <- ncol(y)
  q <- nrow(factors) # factors is q x N

  delta_new <- delta
  inner_prod_f <- rowSums(factors^2) # q-vector: sum_t f_{jt}^2
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  log_prior_odds <- log(tau) - log(1 - tau)
  all_cols <- seq_len(q)

  for (j in seq_len(q)) {
    if (pivots[j] >= V) next
    I_j <- seq.int(pivots[j] + 1L, V)
    n_j <- length(I_j)

    # classifying rows: 
    # "zero": no active loading outside column j (I-a);
    # "nonzero" = at least one other active loading  (I-b).
    other_active <- rowSums(delta_new[I_j, -j, drop = FALSE])
    zero_mask <- other_active == 0L

    O_ij <- numeric(n_j)
    # Posterior shape is alpha[i] + N/2
    posterior_shape <- alpha[I_j] + 0.5 * N

    if (any(zero_mask)) {
      zr <- I_j[zero_mask]
      B_iT <- theta[j] / (1 + theta[j] * inner_prod_f[j])
      D_0 <- -0.5 * log1p(theta[j] * inner_prod_f[j])

      m_iT <- drop(y[zr, , drop = FALSE] %*% factors[j, ]) # sum_t f_{jt} y_{it}
      posterior_rate_null <- beta[zr] + 0.5 * inner_prod_y[zr] # null model scale
      posterior_rate1 <- beta[zr] + 0.5 * (inner_prod_y[zr] - m_iT^2 * B_iT)

      O_ij[zero_mask] <- posterior_shape[zero_mask] * log(posterior_rate_null / posterior_rate1) + D_0
    }

    for (k in which(!zero_mask)) {
      i <- I_j[k]

      # Active set with delta[i,j] forced to 1; j placed last so the
      other_cols <- all_cols[-j][delta_new[i, -j] == 1L]
      ord_cols <- c(other_cols, j)
      q_i <- length(ord_cols)

      X_i <- t(factors[ord_cols, , drop = FALSE]) # N x q_i
      P_i <- crossprod(X_i) + diag(1 / theta[ord_cols], nrow = q_i)
      m_i <- drop(crossprod(X_i, y[i, ]))

      # Single Cholesky P_i = L_i L_i'  (lower triangular)
      L_i <- t(chol(P_i))
      x_i <- forwardsolve(L_i, m_i) # L_i x = m_i

      # posterior scale with delta[i,j]=1  (eq. C.14 / C.9)
      posterior_rate1_i <- beta[i] + 0.5 * (inner_prod_y[i] - sum(x_i^2))

      # set for delta[i,j]=0 via last element of x  (eq. D.10)
      x_star <- x_i[q_i]
      C_0_i <- posterior_rate1_i + 0.5 * x_star^2

      # D_ij  (eq. D.11)
      D_ij_i <- -log(L_i[q_i, q_i]) - 0.5 * log(theta[j])

      O_ij[k] <- posterior_shape[k] * log(C_0_i / posterior_rate1_i) + D_ij_i
    }

    O_post <- O_ij + log_prior_odds[j]
    log_u <- log(runif(n_j))
    cur <- delta_new[I_j, j]

    to_1 <- (cur == 0L) & (log_u <= O_post) # propose 0 -> 1
    to_0 <- (cur == 1L) & (log_u <= -O_post) # propose 1 -> 0

    delta_new[I_j[to_1], j] <- 1L
    delta_new[I_j[to_0], j] <- 0L
  }

  list(delta_new = delta_new)
}
