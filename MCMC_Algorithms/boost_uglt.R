# Boost uglt code
boost_uglt <- function(Lambda, factors, sigma2, theta) {
  # First we sample phi:
  q <- ncol(Lambda)
  V <- nrow(Lambda)
  N <- ncol(factors)
  psi <- diag(1 / rgamma(q, 1.5, rate = 1.5))
  
  factors_psi <- psi^0.5 %*% factors
  Lambda_psi <- Lambda %*% solve(psi)^(0.5)
  
  psi_new_vec <- numeric(q)
  
  
  s <- colSums(Lambda != 0)
  
  for (j in 1:q) {
    psi_new_vec[j] = GIGrvg::rgig(1, lambda = s[j]/2 - 1.5 - N/2,
                                  chi = 3 + t(factors_psi[j, ]) %*% factors_psi[j, ], 
                                  psi = theta[j]^(-1)*t(Lambda_psi[, j]) %*% solve(diag(sigma2)) %*% Lambda_psi[, j])
    if (is.na(psi_new_vec[j])) browser()
  }
  
  # un transofring:
  psi_new <- diag(psi_new_vec)
  factors_new <- solve(psi_new)^(0.5) %*% factors_psi
  Lambda_new <- Lambda_psi %*% psi_new^(0.5)
  return(list(Lambda_new = Lambda_new, factors_new = factors_new))
}