sim_data_2 <- function(N, V, q) {
  Phi <- diag(q)
  Sigma <- diag(V)
  Lambda <- matrix(0, nrow = V, ncol = 2)
  for (i in 1:(V)) {
    for (j in 1:min(i, 2)) {
      if (i == j) {
        # Diagonal entries: positive values
        Lambda[i, j] <- abs(rnorm(1, sd = sqrt(0.5)))
      } else {
        # Off-diagonal entries: can be any value in specified range
        Lambda[i, j] <- rnorm(1, sd = sqrt(0.5))
      }
    }
  }
  data <- matrix(data = NA, nrow = dim(Lambda)[1], ncol = N)
  factors <- matrix(data = NA, nrow = dim(Lambda)[2], ncol = N)
  for (i in 1:N) {
    f <- mvrnorm(1, mu = rep(0, dim(Lambda)[2]), Sigma = Phi)
    factors[, i] <- f
    data[, i] <- mvrnorm(1, mu = Lambda %*% f, Sigma = Sigma)
  }
  return(list(data = data, factors = factors, Lambda = Lambda, Sigma = Sigma, N = N, V = N))
}
