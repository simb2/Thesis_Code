sample_column_shrinkage <- function(theta.shape, theta.rate, Lambda, sigma2, delta) {
  theta.new <- numeric(ncol(delta))
  column.sizes <- colSums(delta)
  rate.update <- colSums(delta*(diag(1/sigma2)%*%Lambda^2))
  for (j in 1:dim(Lambda)[2]) {
    theta.new[j] <- 1/rgamma(1, theta.shape[j] + 0.5*column.sizes[j], theta.rate[j] + 0.5*rate.update[j])
  }
  return(theta.new)
}