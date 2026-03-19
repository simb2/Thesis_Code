# sample loadings for the non sparse plt model

sample_loadings <- function(data, sigma2, Factors, c_0) {
  V <- dim(data)[1]  # number of variables
  q <- dim(Factors)[1]  # number of factors
  N <- dim(data)[2]  # number of observations
  
  Lambda_new <- matrix(NA, V, q)
  
  # iterating the rows of Lambda_new, for h in 1:q first (sampling the first q rows)
  for (h in 1:q) {
    post_cov <- (1/c_0) * diag(h)  # Prior covariance for first h elements
    s <- rep(0, h)  
    
    for (i in 1:N) {
      # updating covariance - use first h factors only
      f_h <- Factors[1:h, i]
      post_cov <- post_cov + (1/sigma2[h]) * (f_h %*% t(f_h))
      s <- s + (1/sigma2[h]) * f_h * data[h, i] 
    }
    
    post_cov <- solve(post_cov)
    post_mean <- post_cov %*% s
    
    # Sample and enforce positive diagonal constraint
    sample_vals <- MASS::mvrnorm(1, post_mean, post_cov)
    
    # first h entries should be lower triangular.
    Lambda_new[h, ] <- c(sample_vals, rep(0, q - h)) 
  }
  
  # iterating the rest of the q+1 to V rows of Lambda_new
  for (h in (q+1):V) {
    post_cov <- (1/c_0) * diag(q)  # Prior covariance for all q factors
    s <- rep(0, q)
    
    for (i in 1:N) {
      f <- Factors[, i]
      post_cov <- post_cov + (1/sigma2[h]) * (f %*% t(f))
      s <- s + (1/sigma2[h]) * f * data[h, i]
    }
    
    post_cov <- solve(post_cov)
    post_mean <- post_cov %*% s
    
    Lambda_new[h, ] <- MASS::mvrnorm(1, post_mean, post_cov)
  }
  
  return(Lambda_new)
}