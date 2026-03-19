sample_tau <- function(hyperparams, delta, pivots) {
  colsums.delta <- colSums(delta)
  new.tau <- numeric(ncol(delta))
  # Recompute the pivots
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  
  for (j in 1:dim(delta)[2]) {
    colsums.delta[j]
    if (hyperparams$bH + nrow(delta) - pivots[j] -  colsums.delta[j] + 1 <  0) {
      browser()
    }
    
    new.tau[j] <- rbeta(1, hyperparams$aH + colsums.delta[j] - 1,
                        hyperparams$bH + nrow(delta) - pivots[j] -  colsums.delta[j] + 1)
    
    if (is.na(new.tau[j])) {
      browser()
    }
  }
  return(new.tau)
}