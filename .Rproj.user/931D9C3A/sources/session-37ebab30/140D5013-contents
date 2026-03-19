# sim_data3

sim_data_3 <- function(N, Lambda, Sigma, factors) {
  data <- matrix(data = NA, nrow = dim(Lambda)[1], ncol = N)
  for (i in 1:N) {
    f <- factors[, i]
    mvrnorm(1, mu = Lambda %*% factors[, i], Sigma = Sigma)
    data[, i] <- mvrnorm(1, mu = Lambda %*% factors[, i], Sigma = Sigma)
  }
  return(data)
}