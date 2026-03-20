source(here('MCMC_Algorithms', 'compute_log_likelihood_ratio.R'))
#' Samples pivots
#'
#' @param delta The sparsity matrix (V x q)
#' @param pivots Current pivot indexes l = (l1, ..., lq)
#' @param factors Current factor draws (q x N matrix)
#' @param y the data matrix (V x N matrix)
#' @param hyperparams List with hyperparameters (aH, bH for 1PB prior)
#' @param move_probs List with probabilities (pshift, pswitch, pa)
#' @param alpha The scale parameters for the IG prior on the observation variances
#' @param beta the rate parameters for the IG prior on the observation variances (TODO maybe make this a list)
#' @return List with updated delta and pivots


update_pivots <- function(delta, pivots, factors, y, theta, alpha, beta,
                          hyperparams, move_probs = list(pshift = 0.3, pswitch = 0.3,
                                                         pa = 0.5) ){
  # Specifies the probability of each 'move type', since move_probs isn't stricly a probability measure.
  probs_type_move <- list(pshift = move_probs$pshift, pswitch = move_probs$pswitch,
                          pad = 1 - move_probs$pshift - move_probs$pswitch)
  # randomize the order in which you draw the pivots
  j_rand <- sample(1:dim(delta)[2])
  
  delta_new <- delta
  
  for (j in j_rand) {
    u <- sample(1:3, 1,  prob = probs_type_move)
    
    if (u == 1) {
      result <- shift_pivot_move(delta, pivots, factors, y, theta, alpha, beta, hyperparams, j)
    } else if (u == 2) {
      result <- switch_pivots_move(delta = delta, pivots = pivots, factors = factors, y = y, 
                                   theta = theta, alpha = alpha, beta = beta, hyperparams = hyperparams, j)
    } else {
      result <- add_delete_pivot_move(delta = delta, pivots = pivots, factors = factors, y = y, theta = theta, alpha = alpha, beta = beta, 
                                      hyperparams = hyperparams, pa = move_probs$pa, j = j)
    }
    # Update the pivots and delta
    pivots <- apply(result$delta, 2, function(col) which(col != 0)[1])
    delta <- result$delta
    
    stopifnot(length(unique(result$pivots)) == length(result$pivots))
    
  }
  
  return(list(pivots = pivots, delta = delta))
}

#' Shift pivot move
#' 
#' @param j the current column index being samples
#' @return a list containing the updated pivots and corresponding sparsity matrix

shift_pivot_move <- function(delta, pivots, factors, y, theta, alpha, beta, hyperparams, j) {
  V <- nrow(delta)
  q <- ncol(delta)
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  lj <- pivots[j]
  
  # if the pivot is equal to V, the column will be spurious.
  if (lj == V) {
    return(list(delta = delta, pivots = pivots))
  }
  
  # Find first non-zero row below lj
  below_pivot <- which(delta[(lj+1):V, j] == 1)[1]
  if (length(below_pivot) == 0) {
    return(list(delta = delta, pivots = pivots))
  }
  l_star <- below_pivot[1] + lj
  if (is.na(l_star)) {
    return(list(delta = delta, pivots = pivots))
  }
  # If l_star <= 2, don't shift (because you won't have any room to do so)
  if (l_star <= 2) {
    return(list(delta = delta, pivots = pivots))
  }
  
  # Available positions for new pivot
  available_positions <- setdiff(1:(l_star-1), pivots)
  if (length(available_positions) == 0) {
    return(list(delta = delta, pivots = pivots))
  }
  
  # Propose new pivot position (sampling uniformly)
  if (length(available_positions) == 1) {
    l_new <- available_positions
  } else {
    l_new <- sample(available_positions, 1)
  }
  if (l_new %in% pivots[-j]) {
    warning("l_new is already a pivot! This should not happen.")
    return(list(delta = delta, pivots = pivots))
  }
  # Calculate log likelihood ratios
  O_new <- compute_log_likelihood_ratio(delta, l_new, j, factors, y, alpha, beta, theta)
  O_old <- compute_log_likelihood_ratio(delta, lj, j, factors, y, alpha, beta, theta)
  # Calculate prior ratio for shift move
  dj <- sum(delta[, j])  # number of non-zero elements in column j
  if (hyperparams$bH + V - l_new - dj + 1 < 0 || hyperparams$bH + V - lj - dj + 1 < 0) {
    browser()
  }
  R_shift <- lbeta(hyperparams$aH + dj - 1, hyperparams$bH + V - l_new - dj + 1) -
    lbeta(hyperparams$aH + dj - 1, hyperparams$bH + V - lj - dj + 1)
  
  # Acceptance probability
  log_alpha <- O_new - O_old + R_shift
  alpha_shift <- exp(log_alpha)
  # Accept or reject
  u = runif(1)
  if (u < alpha_shift) {
    # Update sparsity matrix
    new_delta <- delta
    new_delta[l_new, j] <- 1
    new_delta[lj, j] <- 0
    new_pivots <- apply(new_delta, 2, function(col) which(col != 0)[1])
    
    return(list(delta = new_delta, pivots = new_pivots))
  } else {
    return(list(delta = delta, pivots = pivots))
  }
}



#' Switch pivots move
#'
#' @return List with updated delta and pivots
switch_pivots_move <- function(delta, pivots, factors, y, theta, alpha, beta, hyperparams, j) {
  V <- nrow(delta)
  q <- ncol(delta)
  
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  
  if (dim(delta)[2] <= 1) {
    return(list(delta = delta, pivots = pivots))
  }
  
  lj <- pivots[j]
  
  # Choose another non-zero column
  other_cols <- setdiff(1:q, union(c(j), which(colSums(delta) == 0)))
  
  # If there isn't another column, don't sample
  if (length(other_cols) == 0) {
    return(list(delta = delta, pivots = pivots))
  }
  
  # Sample uniformly over the available columns.
  l <- sample(other_cols, 1)
  l_piv <- pivots[l]
  
  # Find rows between pivots that differ
  min_pivot <- min(lj, l_piv)
  max_pivot <- max(lj, l_piv)
  rows_between <- min_pivot:max_pivot
  
  # Find rows where indicators differ
  S_j_l <- intersect(rows_between, which(delta[, j] != delta[, l]))
  
  if (length(S_j_l) == 0) {
    return(list(delta = delta, pivots = pivots))
  }
  
  # Calculate likelihood ratios for affected rows (i.e over the set wherein the indicators are different)
  delta_cond <- delta
  log_likelihood_sum <- 0
  for (i in S_j_l) {
    # given l
    delta_cond[i, l] <- 0
    O_ij_given_l_piv <- compute_log_likelihood_ratio(delta_cond, i, j, factors, y, alpha, beta, theta)
    delta_cond[i,l] <- delta[i,l]
    # given j
    delta_cond[i, j] <- 0
    O_i_ell_given_j <- compute_log_likelihood_ratio(delta_cond, i, l, factors, y, alpha, beta, theta)
    delta_cond[i, j] <- delta[i, j]
    
    if (delta[i, j] == 0) {
      log_likelihood_sum <- log_likelihood_sum + O_ij_given_l_piv - O_i_ell_given_j
    } else {
      log_likelihood_sum <- log_likelihood_sum + O_i_ell_given_j - O_ij_given_l_piv
    }
  }
  
  # Calculate prior ratio for switch move
  delta_new <- delta
  delta_new[S_j_l, j] <- delta[S_j_l, l]
  delta_new[S_j_l, l] <- delta[S_j_l, j]
  dj_new <- sum(delta_new[, j])
  d_l_piv_new <- sum(delta_new[, l])
  dj <- sum(delta[, j])
  d_l_piv <- sum(delta[, l])
  # After switching, column counts remain the same, only pivot positions change
  if (hyperparams$aH + dj_new - 1 < 0 || hyperparams$bH + V - l_piv - dj_new + 1 < 0) {
    browser()
  }
  if (hyperparams$aH + d_l_piv_new - 1 < 0 || hyperparams$bH + V - l_piv - dj_new + 1 < 0) {
    browser()
  } else if (hyperparams$bH + V - lj - d_l_piv_new + 1 < 0) {
    browser()
  } else if (hyperparams$bH + V - lj - dj + 1 < 0 ) {
    browser()
  } else if (hyperparams$bH + V - l_piv - d_l_piv + 1 < 0) {
    browser()
  }
  R_switch <- lbeta(hyperparams$aH + dj_new - 1, hyperparams$bH + V - l_piv - dj_new + 1) +
    lbeta(hyperparams$aH + d_l_piv_new - 1, hyperparams$bH + V - lj - d_l_piv_new + 1) -
    (lbeta(hyperparams$aH + dj - 1, hyperparams$bH + V - lj - dj + 1) +
       lbeta(hyperparams$aH + d_l_piv - 1, hyperparams$bH + V - l_piv - d_l_piv + 1))
  
  # Acceptance probability
  log_alpha <- log_likelihood_sum + R_switch
  
  alpha_switch <- exp(log_alpha)
  # Accept or reject
  if (runif(1) < alpha_switch) {
    new_pivots <- apply(delta_new, 2, function(col) which(col != 0)[1])
    return(list(delta = delta_new, pivots = new_pivots, accepted = TRUE))
  } else {
    return(list(delta = delta, pivots = pivots, accepted = FALSE))
  }
}

#' Add/delete move
#'
#' @param pa the (tuning) probability of adding a pivot
#' @return List with updated delta and pivots

add_delete_pivot_move <- function(delta, pivots, factors, y, theta, alpha, beta,
                                  hyperparams, pa, j) {
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  
  V <- nrow(delta)
  q <- ncol(delta)
  lj <- pivots[j]
  
  if (lj == V) {
    return(list(delta = delta, pivots = pivots))
  }
  # Check if add move is possible
  if (lj <= 1) {
    available_for_add <- NULL
  } else {
    available_for_add <- setdiff(1:max(lj-1, 1), pivots)
  }
  can_add <- length(available_for_add) > 0
  
  # Check if delete move is possible
  below_pivot <- which(delta[(lj+1):V, j] == 1)
  can_delete <- length(below_pivot) > 0
  
  if (can_delete) {
    l_star <- below_pivot[1] + lj
    delta_new <- delta
    delta_new[lj, j] <- 0
    piv_new <- apply(delta_new, 2, function(col) which(col != 0)[1])
    condition2 <- (length(unique(piv_new)) == length(piv_new))
    can_delete <- condition2
  }
  
  # If I can't add or can't delete
  if (!can_add && !can_delete) {
    return(list(delta = delta, pivots = pivots))
  }
  
  # Decide on add vs delete
  p_add_delta <- NA
  if (can_add && can_delete) {
    p_add_delta <- pa
  } else if (can_add) {
    p_add_delta <- 1
  } else {
    p_add_delta <- 0
  }
  
  if (runif(1) < p_add_delta) {
    # Add move
    return(add_pivot_move(delta = delta, pivots = pivots, j = j, factors = factors, y = y, 
                          theta = theta, alpha = alpha, beta = beta,
                          hyperparams = hyperparams, p_add_delta = p_add_delta,
                          available_positions = available_for_add, pa = pa))
  } else {
    # Delete move
    return(delete_pivot_move(delta, pivots, j, l_star, factors, y, theta, alpha, beta,
                             hyperparams, p_add_delta, pa))
  }
}


#' Add pivot move
#' @param p_add_delta the precomputed probability of adding a pivot (given pa)
#' @param available_positions the sample space (a list of row indeces available)
#' @return List with updated delta and pivots

add_pivot_move <- function(delta, pivots, j, factors, y, theta, alpha, beta,
                           hyperparams, available_positions, p_add_delta, pa) {
  V <- nrow(delta)
  q <- ncol(delta)
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  
  lj <- pivots[j]
  
  # Propose a new position uniformly
  if(length(available_positions) == 1) {
    l_new <- available_positions
  } else {
    l_new <- sample(available_positions, 1)
  }
  pivots_new <- pivots
  pivots_new[j] <- l_new
  
  # Calculate the 'reversibale' part here
  if (l_new == 1) {
    available_positions_new <- NULL
  } else{
    available_positions_new <- setdiff(1:(l_new - 1), pivots_new[-j])
  }
  
  p_add_new <- NA
  if (length(available_positions_new) == 0) {
    p_add_new <- 0
  } else {
    p_add_new <- pa
  }
  
  proposal_ratio <- ((1 - p_add_new) * length(available_positions)) / (p_add_delta)
  
  # Calculate likelihood ratio
  delta_new <- delta
  delta_new[l_new, j] <- 1
  O_new <- compute_log_likelihood_ratio(delta_new, l_new, j, factors, y, alpha, beta, theta)
  
  # Calculate prior ratio for add move
  dj <- sum(delta[, j])
  R_add <- lbeta(hyperparams$aH + dj, hyperparams$bH + V - l_new - dj) -
    lbeta(hyperparams$aH + dj - 1, hyperparams$bH + V - lj - dj + 1)
  
  # Note: p_add_new will not be 1, was can delete this pivot (being the previous pivot row incolumn j)
  
  proposal_ratio <- ((1 - p_add_new)* length(available_positions)) / (p_add_delta)
  
  # Acceptance probability
  log_alpha <- O_new + R_add + log(proposal_ratio)
  alpha_add <- exp(log_alpha)
  
  # Accept or reject
  if (runif(1) < alpha_add) {
    new_delta <- delta
    new_delta[l_new, j] <- 1
    pivots_new <- apply(new_delta, 2, function(col) which(col != 0)[1])
    return(list(delta = new_delta, pivots = pivots_new))
  } else {
    pivots <- apply(delta, 2, function(col) which(col != 0)[1])
    
    return(list(delta = delta, pivots = pivots))
  }
}

#' Delete pivot move

delete_pivot_move <- function(delta, pivots, j, l_star, factors, y, theta, alpha, beta,
                              hyperparams, p_add_delta, pa) {
  V <- nrow(delta)
  pivots <- apply(delta, 2, function(col) which(col != 0)[1])
  lj <- pivots[j]
  
  # Calculate likelihood ratio (negative of add move)
  O_delete <- - compute_log_likelihood_ratio(delta, lj, j, factors, y, alpha, beta, theta)
  # Calculate prior ratio for delete move (equation G.5)
  dj <- sum(delta[, j])
  R_delete <- lbeta(hyperparams$aH + dj - 2, hyperparams$bH + V - l_star - dj + 2) /
    lbeta(hyperparams$aH + dj - 1, hyperparams$bH + V - lj - dj + 1)
  
  # computing the 'reversabilization' part
  available  <- setdiff(1:(l_star-1), pivots[-j])
  p_add_new <- ifelse(length(available) > 0, pa, 0)
  proposal_ratio <- (p_add_new) / ((1 - p_add_delta)* length(available))
  
  # Acceptance probability
  log_alpha <- O_delete + R_delete + log(proposal_ratio)
  alpha_delete <- min(1, exp(log_alpha))
  # Accept or reject
  if (runif(1) < alpha_delete) {
    new_delta <- delta
    new_delta[lj, j] <- 0
    
    new_pivots <- apply(new_delta, 2, function(col) which(col != 0)[1])
    
    return(list(delta = new_delta, pivots = new_pivots))
  } else {
    pivots <- apply(delta, 2, function(col) which(col != 0)[1])
    return(list(delta = delta, pivots = pivots))
  }
}



