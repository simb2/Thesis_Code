library(coda)
par(mfrow = c(1, 2))

set.seed(42)
n <- 5000

target_density <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)

mh_sampler <- function(n_iter, init, proposal_sd) {
  chain <- numeric(n_iter)
  chain[1] <- init
  accepted <- 0
  for (i in 2:n_iter) {
    proposal <- rnorm(1, mean = chain[i - 1], sd = proposal_sd)
    log_alpha <- target_density(proposal) - target_density(chain[i - 1])
    if (log(runif(1)) < log_alpha) {
      chain[i] <- proposal
      accepted  <- accepted + 1
    } else {
      chain[i] <- chain[i - 1]
    }
  }
  cat("Acceptance rate:", round(accepted / n_iter, 3), "\n")
  chain
}

good_mh <- mh_sampler(n_iter = n, init = 0, proposal_sd = 2)
acf(good_mh,
    main = "",
    col  = "steelblue",
    lwd  = 2)

bad_mh <- mh_sampler(n_iter = n, init = 0, proposal_sd = 0.1)
acf(bad_mh,
    main = "",
    col  = "firebrick",
    lwd  = 2)

