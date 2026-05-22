---
title: "profiling_code"
format: html
editor: visual
---

## Installing libraries, sourcing other files 

```{r}
library(profvis)
library(here)
library(tidyverse)
library(reshape2)
library(MASS)
library(scoringRules) # CRPS
```


Generating data:

```{r}
library(here)

source(here("MCMC_Algorithms", "run_mcmc_UGLT.R"))
source(here("Run_Simulations" ,"sim_data_3.R"))

# Prior specifications and Globals ----------------------------------------
set.seed(8)
n_subj <- c(100)
n_vars <- c(50)
n_factors <- 10

# Helper Functions --------------------------------------------------------

# First a function to simulate the data

sim_data_UGLT <- function(N, V, q) {
  # priors for theta and sigma respecively
  # True value for sigma
  sigma2 <- rep(1, V)
  Sigma <- diag(sigma2)
  theta <- rep(1, q) # L2 shrinkage
  tau <- (1:q) / (2*q + 1) # L1 shrinkage / spike probability
  hyperparams <- list(aH = 2, bH = 2) # priors for tau
  pivots <- sample(1:(V-3), q)
  pivots <- pivots[order(pivots)]
  pivots <- 1:q # for simulating a plt structure. 
  delta <- matrix(0, nrow = V, ncol = q)
  Phi <- diag(q)
  for (j in 1:q) {
    for (i in pivots[j]:V) {
      if (i == pivots[j]) {
        delta[i, j] <- 1
      } else {
        if (runif(1) > tau[j]) {
          delta[i, j] <- 1
        }
      }
    }
  }
  # now sample Lambda according to the prior.
  
  Lambda <- matrix(0, nrow = V, ncol = q)
  # First computing
  for (i in 1:V) {
    # First we make sure the row is non zero
    if (sum(delta[i, ]) != 0) {
      filter <- which(delta[i, ] == 1)
      theta_a <- theta[filter]
      L_0 <- diag(length(theta_a)) * theta_a
      Lambda[i, filter] <- mvrnorm(1, mu = rep(0, sum(delta[i, ])), Sigma = L_0 * sigma2[i])
    }
  }
  Lambda <- Lambda%*%diag(diag(sign(Lambda[pivots, ])))
  data <- matrix(data = NA, nrow = dim(Lambda)[1], ncol = N)
  factors <- matrix(data = NA, nrow = dim(Lambda)[2], ncol = N)
  for (i in 1:N) {
    f <- MASS::mvrnorm(1, mu = rep(0, dim(Lambda)[2]), Sigma = Phi)
    factors[, i] <- f
    data[, i] <- MASS::mvrnorm(1, mu = Lambda %*% f, Sigma = Sigma)
  }
  return(list(
    data = data, factors = factors, Lambda = Lambda, Sigma = Sigma, N = N, V = V,
    delta = delta, pivots = pivots, theta = theta)
  )
}

get_starting_vals <- function(samp) {
  
  data <- samp$data
  V <- samp$V
  
  return(list(
    N = samp$N, q = nrow(samp$factors) + 3, n_runs = 6000, alpha = rep(1.5, V), beta = rep(1.5, V), 
    theta.shape = rep(1.5, V), theta.rate = rep(1.5, V), hyperparams = list(aH = 2, bH = 2), data = data, thin = 2, burn = 1000
  ))
}

```



```{r}

##############################################################

# Here is where I actually run the simulation

settings <- tidyr::crossing(n_subj, n_vars, n_factors)
names(settings) <- c("N", "V", "q")
samples <- purrr::pmap(settings, sim_data_UGLT)
starting_vals <- purrr::map(samples, get_starting_vals)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_UGLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_ests <- purrr::map(post_samps, "estimates")
trace_acf <- purrr::map(post_draws, make_trace_acf_plot)

```