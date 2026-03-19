library(tidyverse)
source(here('MCMC_Algorthims' ,'compute_log_likelihood_ratio.R'))
#' Samples the sparsity matrix
#'
#'
#' @param delta The current sparsity matrix (V x q)
#' @param y the data matrix (V x N)
#' @param factors a matrix containing the factors (q x N)
#' @param tau a vector containing spike and slab shrinkage (q x 1)
#' @param theta a vector containing the column wise L2 shrinkage (q x 1)
#' @param alpha a vector containing the shape parameter for the priors on the sigma's (V x 1)
#' @param beta a vector for the scale parameters for the prior on the sigma's (V x 1)


sample_sparsity <- function(y, factors, tau, theta, delta, alpha, beta) {
  accepted = 0
  V <- dim(y)[1]
  q <- dim(factors)[1]
  delta_new <- delta # this is to ensure that I can mutate this
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  for (j in 1:q) {
    if (pivots[j] != V) {
      for (i in (pivots[j] + 1):V) {
        O_ij <- compute_log_likelihood_ratio(delta_new, i, j, factors, y, alpha, beta, theta)
        PO <- O_ij + log(tau[j]) - log(1 - tau[j])
        if (is.na(PO)) browser()
        u <- runif(1)
        if (log(u) < PO) {
          delta_new[i,j] <- 1
          accepted = accepted + 1
        } else {
          delta_new[i,j] <- 0
        }
      }
    }
  }
  return(list(delta_new = delta_new, accepted = accepted))
}
