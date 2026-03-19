# Code for running the non-sparse plt model
library(MASS)
library(cli)
library(here)
source(here('MCMC_Algorithms', 'sample_factors_b.R'))
source(here('MCMC_Algorithms', 'sample_loadings_b.R'))
source(here('MCMC_Algorithms', 'sample_variances_b.R'))
source(here('MCMC_Algorithms', 'boost.R'))
run_mcmc <- function(Lambda_0, sigma2_0, n_runs = 1000, data, q, nu, s2, c_0, thin = 1, burn = 1) {
  cli_progress_bar("Sampling from Posterior", total = n_runs)
  W <- array(data = NA, dim = c(n_runs, q, dim(data)[2]))
  L <- array(data = NA, dim = c(n_runs, dim(data)[1], q))
  L[1, , ] <- Lambda_0
  S <- array(data = NA, dim = c(n_runs, dim(data)[1]))
  S[1, ] <- sigma2_0
  T_stat <- numeric(n_runs)
  for (c in 1:n_runs) {
    if (c > 1) {
      W[c, , ] <- sample_factors(L[c-1, , ], S[c-1,], data, q)
      L[c, , ] <- sample_loadings(data, S[c-1, ], W[c-1, , ], c_0)
    } else {
      W[c, , ] <- sample_factors(L[1, , ], S[1,], data, q)
      L[c, , ] <- sample_loadings(data, S[1, ], W[1, , ], c_0)
    }
    S[c, ] <- sample_variances(nu, s2, data, W[c, , ], L[c, , ])
    
    boost <- boost(Lambda = L[c, , ], factors = W[c, , ], c_0)
    L[c, , ] <- boost$Lambda_new
    W[c, , ] <- boost$factors_new
    T_stat[c] <- sum(diag(((L[c, , ])  %*% t(L[c, , ]))) + solve(diag(S[c, ])))
    cli_progress_update()
  }
  for (l in 1:n_runs) {
    for (i in 1:dim(data)[2]) {
      W[l, , i] <- diag(sign(diag(L[l, , ]))) %*% W[l, , i]
    }
    L[l , , ] <- L[l , , ] %*% diag(sign(diag(L[l, , ])))
  }
  thin_burn <- seq(burn, n_runs, by = thin)
  return(list(Factors = W[thin_burn, , ], Lambda = L[thin_burn, , ], 
              Variances = S[thin_burn, ], T_stat = T_stat[thin_burn]))
}


# test matrix