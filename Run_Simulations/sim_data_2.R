sim_data_2 <- function(N, V, q) {
  Phi <- diag(q)
  Sigma <- diag(V)
  Lambda <- matrix(0, nrow = V, ncol = 2)
  for (i in seq_len(V)) {
    for (j in seq_len(min(i, 2))) {
      Lambda[i, j] <- if (i == j) abs(rnorm(1, sd = sqrt(0.5))) else rnorm(1, sd = sqrt(0.5))
    }
  }
  factors <- t(MASS::mvrnorm(N, rep(0, 2), Phi))
  data <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), Sigma))
  list(data = data, factors = factors, Lambda = Lambda, Sigma = Sigma, N = N, V = V)
}
