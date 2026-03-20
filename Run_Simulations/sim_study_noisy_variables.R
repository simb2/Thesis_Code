# Sim study analysis
library(tidyverse)
library(reshape2)
library(here)
library(MASS)
library(cli)
library(scoringRules) # CRPS
setwd(here())


source(here("Run_Simulations", "sim_data_2.R"))
source(here("MCMC_Algorithms", "run_mcmc_PLT.R"))  
source(here("Run_Simulations", "sim_data_3.R"))
source(here("Run_Simulations", 'sim_study_noisy_variables_helpers.R'))

n_subj <- c(200, 400)
n_vars <- c(10, 20)
n_factors <- 2
set.seed(8)


settings <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors)
samples <- pmap(settings, sim_data_L)
data <- lapply(samples, function(x) x$data)



# ------------------------- Running the spurious factor model -------------- #
starting_values <- map(data, get_starting_vals)
post_samps <- map(starting_values, ~ do.call(run_mcmc, .x))
posterior_estimates <- map(post_samps, compute_posterior_estimates)

gof_inputs <- settings |> mutate(
  Lambda = map(samples, "Lambda"),
  Sigma = map(samples, "Sigma"),
  factors = map(samples, "factors"),
  Lambda_est = map(posterior_estimates, "Lambda_estimates"),
  sigma_est = map(posterior_estimates, "sigma_est"),
  factors_est = map(posterior_estimates, "factors_est")
)

mse_vals <- pmap(gof_inputs, compute_goodness_of_fit)
gof_inputs$mse <- mse_vals

crps_inputs <- tibble(
  true_val = gof_inputs$Lambda,
  lambda_samps1 = map(post_samps, "Lambda")
)

crps <- pmap(crps_inputs, avg_crps)
gof_inputs$crps <- unlist(crps)


# Making plots

est_heat_map <- purrr::map(gof_inputs$Lambda_est, make_heat_map)
true_heat_map <- purrr::map(gof_inputs$Lambda, make_heat_map)
trace_acf <- purrr::map(post_samps, make_trace_acf_plot)

# Creating the overall data frame to save everything in
contr_noise_inputs <- tibble(
  L = map(posterior_estimates, "Lambda_estimates"),
  y = map(samples, "data")
)

noise <- pmap(contr_noise_inputs, contributed_variance_noise)

tru_contr_noise_inputs <- tibble(
  L = map(samples, "Lambda"), 
  y = map(samples, "data")
)

tru_noise <- pmap(tru_contr_noise_inputs, contributed_variance_noise)

sim_results <- gof_inputs |> mutate(
  est_heat_map = est_heat_map,
  true_heat_map = true_heat_map,
  trace_acf = trace_acf,
  post_samps = post_samps,
  noise = noise, 
  tru_noise = tru_noise
)

gof_inputs$noise = noise

saveRDS(gof_inputs, "analysis_summary_no_nsp2.rds")
saveRDS(sim_results, "results_no_nsp2.rds")


########################################
########################################
########################################

get_starting_vals <- function(samp) {
  
  data <- samp$data
  V <- samp$V
  
  return(list(
    N = samp$N, q = nrow(samp$factors), n_runs = 2000, alpha = rep(1.5, V), beta = rep(1.5, V), 
    theta.shape = rep(1.5, V), theta.rate = rep(1.5, V), hyperparams = list(aH = 2, bH = 2), data = data, thin = 2, ident = TRUE,
    burn = 500
  ))
  
}


##############################################
##############################################
##############################################


source(here("MCMC_Algorithms" ,'run_mcmc_UGLT_Static.R'))
settings <- tidyr::crossing(n_subj, n_vars, n_factors)
settings_uglt <- settings
samples_uglt <- samples
names(settings_uglt) <- c("N", "V", "q")

starting_vals <- purrr::map(samples_uglt, get_starting_vals)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_UGLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_estimates <- purrr::map(post_samps, "estimates")

make_trace_acf_plot <- function(post_draws) {
  T_stat_df <- tibble(
    Value = post_draws$T_stat,
    Iteration = 1:length(post_draws$T_stat)
  )
  acf <- acf(post_draws$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  return(list(acf = acf, trace = trace))
}

trace_acf <- map(post_draws, make_trace_acf_plot)

# Computing CRPS 
estimates <- list()
for (i in 1:4) {
  estimates[[i]] <- post_estimates[[i]][[1]]
}




gof_inputs <- settings_uglt |> mutate(
  Lambda = map(samples_uglt, "Lambda"),
  Sigma = map(samples_uglt, "Sigma"),
  factors = map(samples_uglt, "factors"),
  Lambda_est = map(estimates, "lambda_est"),
  sigma_est = map(estimates, "sigma2_mean"),
  factors_est = map(estimates, "factors_est")
)


mse_vals <- pmap(gof_inputs, compute_goodness_of_fit)
mse_vals

gof_inputs$mse <- mse_vals

crps_inputs <- tibble(
  true_val = gof_inputs$Lambda,
  lambda_samps1 = map(post_draws, "Lambda_test")
)

avg_crps_list <- function(lambda_samps1, true_val) {
  # compute the avg crps
  
  V <- nrow(lambda_samps1[[1]])
  q <- ncol(lambda_samps1[[1]])
  lambda_samps <- array(data = NA, dim = c(length(lambda_samps1), V, q))
  for (i in seq_along(lambda_samps1)) {
    lambda_samps[i, , ] <- lambda_samps1[[i]]
  }
  
  crps <- 0
  for (i in 1:V) {
    for (j in 1:q) {
      crps <- crps + crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
    }
  }
  return(crps / (V * q))
}
# Creating the overall data frame to save everything in


crps <- pmap(crps_inputs, avg_crps_list)


contr_noise_inputs <- tibble(
  L = map(estimates, "lambda_est"),
  y = map(samples_uglt, "data")
)
noise <- pmap(contr_noise_inputs, contributed_variance_noise)

gof_inputs$crps <- unlist(crps)

gof_inputs <- gof_inputs %>% mutate(
  noise = noise
)
# Making plots

est_heat_map <- map(gof_inputs$Lambda_est, make_heat_map)
true_heat_map <- map(gof_inputs$Lambda, make_heat_map)
# trace_acf <- map(post_samps, make_trace_acf_plot)

# Creating the overall data frame to save everything in

sim_results <- gof_inputs |> mutate(
  post_samps = post_samps
)

saveRDS(sim_results, "results_sp2_UGLT.rds")
saveRDS(gof_inputs, "analysis_summary_sp2_UGLT.rds")

#################################################################
#################################################################
#################################################################



source(here("MCMC_Algorithms" ,'run_mcmc_sparse_PLT_Static.R'))
settings <- tidyr::crossing(n_subj, n_vars, n_factors)
names(settings) <- c("N", "V", "q")
settings_sparse_plt <- settings
samples_sparse_plt <- samples
starting_vals <- purrr::map(samples_sparse_plt, get_starting_vals)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_sparse_PLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_estimates <- purrr::map(post_samps, "estimates")

make_trace_acf_plot <- function(post_draws) {
  T_stat_df <- tibble(
    Value = post_draws$T_stat,
    Iteration = 1:length(post_draws$T_stat)
  )
  acf <- acf(post_draws$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  return(list(acf = acf, trace = trace))
}

trace_acf <- map(post_draws, make_trace_acf_plot)
trace_acf
# Computing CRPS 
estimates <- list()
for (i in 1:4) {
  estimates[[i]] <- post_estimates[[i]][[1]]
}




gof_inputs <- settings_sparse_plt |> mutate(
  Lambda = map(samples_sparse_plt, "Lambda"),
  Sigma = map(samples_sparse_plt, "Sigma"),
  factors = map(samples_sparse_plt, "factors"),
  Lambda_est = map(estimates, "lambda_est"),
  sigma_est = map(estimates, "sigma2_mean"),
  factors_est = map(estimates, "factors_est")
)


mse_vals <- pmap(gof_inputs, compute_goodness_of_fit)
mse_vals

gof_inputs$mse <- mse_vals

crps_inputs <- tibble(
  true_val = gof_inputs$Lambda,
  lambda_samps1 = map(post_draws, "Lambda_test")
)

avg_crps_list <- function(lambda_samps1, true_val) {
  # compute the avg crps
  
  V <- nrow(lambda_samps1[[1]])
  q <- ncol(lambda_samps1[[1]])
  lambda_samps <- array(data = NA, dim = c(length(lambda_samps1), V, q))
  for (i in seq_along(lambda_samps1)) {
    lambda_samps[i, , ] <- lambda_samps1[[i]]
  }
  
  crps <- 0
  for (i in 1:V) {
    for (j in 1:q) {
      crps <- crps + crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
    }
  }
  return(crps / (V * q))
}
# Creating the overall data frame to save everything in


crps <- pmap(crps_inputs, avg_crps_list)


contr_noise_inputs <- tibble(
  L = map(estimates, "lambda_est"),
  y = map(samples_sparse_plt, "data")
)
noise <- pmap(contr_noise_inputs, contributed_variance_noise)

gof_inputs$crps <- unlist(crps)

gof_inputs <- gof_inputs %>% mutate(
  noise = noise
)
# Making plots

est_heat_map <- map(gof_inputs$Lambda_est, make_heat_map)
true_heat_map <- map(gof_inputs$Lambda, make_heat_map)
# trace_acf <- map(post_samps, make_trace_acf_plot)

# Creating the overall data frame to save everything in

sim_results <- gof_inputs |> mutate(
  post_samps = post_samps
)

saveRDS(sim_results, "results_sp2_Sparse_PLT.rds")
saveRDS(gof_inputs, "analysis_summary_sp2_Sparse_PLT.rds")
