boost <- function(Lambda, factors, c_0) {
  # First we sample phi:
  q <- ncol(Lambda)
  V <- nrow(Lambda)
  N <- ncol(factors)
  psi <- diag(1 / rgamma(q, 1.5, rate = 1.5))
  
  factors_psi <- psi^0.5 %*% factors
  Lambda_psi <- Lambda %*% solve(psi)^(0.5)
  
  psi_new_vec <- numeric(q)
  
  for (j in 1:q) {
    psi_new_vec[j] = GIGrvg::rgig(1, lambda = (V - j + 1)/2 - 1.5 - N/2,
                                  chi = 3 + t(factors_psi[j, ]) %*% factors_psi[j, ], psi = c_0^(-1) * t(Lambda_psi[, j]) %*% Lambda_psi[, j])
    if (is.na(psi_new_vec[j])) browser()
  }
  
  # un transofring:
  psi_new <- diag(psi_new_vec)
  factors_new <- solve(psi_new)^(0.5) %*% factors_psi
  Lambda_new <- Lambda_psi %*% psi_new^(0.5)
  return(list(Lambda_new = Lambda_new, factors_new = factors_new))
}
