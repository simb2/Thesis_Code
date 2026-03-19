# Here, we filter out any "un-nessassary factors, and kind of transform them by not only droping the columns, 
# but also making sure we add the spurious factor to the variance component


filter_factors <- function(pivot_draws, factor_draws, lambda_draws, delta_draws, 
                           theta_draws, sigma_draws, tau_draws, n_draws) {
  # First we identify the spurious columns
  for (i in 1:n_draws) {
    spur_col <- colSums(delta_draws[[i]]) == 1
    if (length(spur_col) > 0) {
      # we merge the spurious column into the idiosyncratic variance:
      col_sp <- which(spur_col)
      k = 0
      for (col in col_sp) {
        # adding the variance component to 'sigma'
        # first we need to find the row
        col = col - k
        row <- which(delta_draws[[i]][, col] == 1)[1]
        sigma_draws[[i]][row] <- sigma_draws[[i]][row]
        factor_draws[[i]] <- factor_draws[[i]][-col, ]
        lambda_draws[[i]] <- lambda_draws[[i]][ , -col]
        delta_draws[[i]] <- delta_draws[[i]][, -col]
        if (is.null(dim(delta_draws[[i]]))) {
          print(delta_draws[[i]])
          
          pivot_draws[[i]] <- which(delta_draws[[i]] == 1)[1]
        } else {
          pivot_draws[[i]] <- apply(delta_draws[[i]], 2, function(j) which(j != 0)[1])
        }
        theta_draws[[i]] <- theta_draws[[i]][-col]
        tau_draws[[i]] <- tau_draws[[i]][-col]
        k = k + 1
      }
    }
  }
  return(list(sigma_draws = sigma_draws, factor_draws = factor_draws, 
              lambda_draws = lambda_draws, delta_draws = delta_draws,
              tau_draws = tau_draws, theta_draws = theta_draws, pivot_draws = pivot_draws))
}