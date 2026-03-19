# sample loadings for the non sparse plt model

sample_factors <- function(Lambda, sigma2, data, q) {
  SigmaInv <- diag(sigma2^(-1))
  factors <- matrix(NA, nrow = q, ncol = dim(data)[2])
  
  Vf <- solve(diag(q) + t(Lambda) %*% SigmaInv %*% Lambda)
  
  for(i in 1:dim(data)[2]) {
    mean_vec <- Vf %*% t(Lambda) %*% SigmaInv %*% data[, i]
    factors[, i] <- MASS::mvrnorm(1, mu = mean_vec, Sigma = Vf)
  }
  return(factors)
}