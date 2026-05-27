// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]

double compute_log_likelihood_ratio <- function(umat delta, uword l_new, uword j, mat factors, mat y, 
  vec alpha, vec beta, vec theta) {
  log_lik <- 0.5*(log(det(L_iN)) - log(det(L_0))) - 0.5*N*log(2*pi) +
    alpha[i]*log(beta[i]) - lgamma(alpha[i])+ lgamma(N/2 + alpha[i]) -
    (N/2 + alpha[i])*log(beta[i] * 0.5*(t(y[i, ])  %*% (y[i, ]) -
    t(M_i) %*% solve(L_iN) %*% M_i))
  
      log_lik_null <- 0.5*(log(det(L_iN)) - log(det(L_0))) - 0.5*N*log(2*pi) +
      alpha[i]*log(beta[i]) - lgamma(alpha[i])+ lgamma(N/2 + alpha[i]) -
      (N/2 + alpha[i])*log(beta[i] * 0.5*(t(y[i, ])  %*% y[i, ] -
                                            t(M_i) %*% solve(L_iN) %*% M_i))