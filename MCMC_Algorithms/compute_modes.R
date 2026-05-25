#' Compute the posterior mode of \eqn{\delta}
#'
#' Finds the most frequently visited sparsity pattern across a set of draws
#' with the same number of factors r.
#'
#' @param x Tibble of draws (filtered to a fixed r), with a \code{delta_test} list-column.
#' @return \eqn{v \times r} matrix of the modal sparsity pattern.
post_mode_delta <- function(x) {
  delta_test <- array(NA, dim = c(length(x$delta_test), 
                                  nrow(x$delta_test[[1]]), # V
                                  ncol(x$delta_test[[1]])))
  
  for (i in seq_along(x$delta_test)) {
    delta_test[i, , ] <- x$delta_test[[i]]
  }
  
  counted_sparsities <- as_tibble(delta_test) |>
    group_by_all() |>
    summarise(n = n())
  post_mode_delta <- as.numeric(
    counted_sparsities[which.max(counted_sparsities$n), ]
  )
  matrix(post_mode_delta[1:(nrow(x$delta_test[[1]]) * ncol(x$delta_test[[1]]))], 
         nrow = nrow(x$delta_test[[1]]))
}

#' Compute the posterior mean of \eqn{\Lambda} under the modal \eqn{\delta}
#'
#' Averages loading matrices over draws where the sparsity pattern matches
#' \code{post_mode_delta} exactly.
#'
#' @param x Tibble of draws with \code{Lambda_test} and \code{delta_test} list-columns.
#' @param post_mode_delta \eqn{v \times r} modal sparsity matrix from \code{post_mode_delta()}.
#' @return \eqn{v \times r} posterior mean loading matrix.
post_est_lambda <- function(x, post_mode_delta) {
  matches <- which(vapply(x$delta_test, function(d) all(d == post_mode_delta), logical(1)))
  if (length(matches) == 0)
    return(matrix(0, nrow(post_mode_delta), ncol(post_mode_delta)))
  matched <- x$Lambda_test[matches]
  Reduce("+", matched) / length(matched)
}

#' Compute posterior means of \eqn{\sigma^2}, \eqn{\tau}, \eqn{\theta}, and \eqn{F}
#'
#' Averages \eqn{\sigma^2}, \eqn{\tau}, \eqn{\theta}, and \eqn{F} across all draws in \code{x}.
#'
#' @param x Tibble of draws with list-columns \code{sigma_test}, \code{tau_test},
#'   \code{theta_test}, and \code{W}.
#' @return List with \code{sigma2_mean}, \code{tau_mean}, \code{theta_mean}, \code{factors_est}.
compute_post_means <- function(x) {
  sigma2_mean <- numeric(length = length(x$sigma_test[[1]]))
  tau_mean <- theta_mean <- numeric(length = length(x$tau_test[[1]]))
  factors_est <- matrix(0, nrow = nrow(x$W[[1]]), ncol = ncol(x$W[[1]]))
  for (i in seq_along(x$Lambda_test)) {
    sigma2_mean <- sigma2_mean + x$sigma_test[[i]]
    tau_mean <- tau_mean + x$tau_test[[i]]
    theta_mean <- theta_mean + x$theta_test[[i]]
    factors_est <- factors_est + x$W[[i]]
  }
  sigma2_mean <- (sigma2_mean)/length(x$Lambda_test)
  tau_mean <- (tau_mean)/length(x$Lambda_test)
  theta_mean <- (theta_mean)/length(x$Lambda_test)
  factors_est <- (factors_est)/length(x$Lambda_test)
  
  
  list(sigma2_mean = sigma2_mean, tau_mean = tau_mean, theta_mean = theta_mean, 
       factors_est = factors_est)
}
