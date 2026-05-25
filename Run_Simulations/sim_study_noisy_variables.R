library(tidyverse)
library(MASS)
library(cli)
library(scoringRules)
library(here)
setwd(here())

source(here("Run_Simulations", "sim_data_3.R"))
source(here("Run_Simulations", "sim_helpers.R"))
source(here("Run_Simulations", "sim_study_noisy_variables_helpers.R"))
source(here("MCMC_Algorithms", "run_mcmc_PLT.R"))

n_subj <- c(200, 400)
n_vars <- c(10, 20)
n_factors <- 2
set.seed(8)

settings <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors)
samples <- pmap(settings, sim_data_L)
data_list <- lapply(samples, function(x) x$data)

# ---- Non-sparse PLT model ------------------------------------------------
starting_values <- map(data_list, get_starting_vals)
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
gof_inputs$mse <- pmap(gof_inputs, compute_goodness_of_fit)
gof_inputs$crps <- unlist(pmap(
  tibble(lambda_samps = map(post_samps, "Lambda"), true_val = gof_inputs$Lambda),
  avg_crps_array
))

noise <- pmap(
  tibble(L = map(posterior_estimates, "Lambda_estimates"), y = map(samples, "data")),
  contributed_variance_noise
)
tru_noise <- pmap(
  tibble(L = map(samples, "Lambda"), y = map(samples, "data")),
  contributed_variance_noise
)

sim_results <- gof_inputs |> mutate(
  est_heat_map = purrr::map(Lambda_est, make_heat_map),
  true_heat_map = purrr::map(Lambda, make_heat_map),
  trace_acf = purrr::map(post_samps, make_trace_acf_plot),
  post_samps = post_samps,
  noise = noise,
  tru_noise = tru_noise
)
gof_inputs$noise <- noise

saveRDS(gof_inputs, "analysis_summary_no_nsp2.rds")
saveRDS(sim_results, "results_no_nsp2.rds")

# ---- UGLT model ----------------------------------------------------------
source(here("MCMC_Algorithms", "run_mcmc_UGLT_Static.R"))
starting_vals <- purrr::map(samples, get_starting_vals_uglt, n_runs = 5000, thin = 2, burn = 500)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_UGLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_estimates <- purrr::map(post_samps, "estimates")
estimates <- purrr::map(post_estimates, 1)
trace_acf <- purrr::map(post_draws, make_trace_acf_plot)

gof_inputs <- settings |> mutate(
  Lambda = map(samples, "Lambda"),
  Sigma = map(samples, "Sigma"),
  factors = map(samples, "factors"),
  Lambda_est = map(estimates, "lambda_est"),
  sigma_est = map(estimates, "sigma2_mean"),
  factors_est = map(estimates, "factors_est")
)
gof_inputs$mse <- pmap(gof_inputs, compute_goodness_of_fit)
gof_inputs$crps <- unlist(pmap(
  tibble(lambda_samps = map(post_draws, "Lambda_test"), true_val = gof_inputs$Lambda),
  avg_crps
))
gof_inputs$noise <- pmap(
  tibble(L = map(estimates, "lambda_est"), y = map(samples, "data")),
  contributed_variance_noise
)

sim_results <- gof_inputs |> mutate(
  est_heat_map  = map(Lambda_est, make_heat_map),
  true_heat_map = map(Lambda, make_heat_map),
  trace_acf = trace_acf,
  post_samps = post_samps
)
saveRDS(sim_results, "results_sp2_UGLT.rds")
saveRDS(gof_inputs, "analysis_summary_sp2_UGLT.rds")

# ---- Sparse PLT model ----------------------------------------------------
source(here("MCMC_Algorithms", "run_mcmc_sparse_PLT_Static.R"))
starting_vals <- purrr::map(samples, get_starting_vals_uglt, n_runs = 5000, thin = 2, burn = 500)
post_samps <- purrr::map(starting_vals, ~ do.call(run_mcmc_sparse_PLT, .x))
post_draws <- purrr::map(post_samps, "draws")
post_estimates <- purrr::map(post_samps, "estimates")
estimates <- purrr::map(post_estimates, 1)
trace_acf <- purrr::map(post_draws, make_trace_acf_plot)

gof_inputs <- settings |> mutate(
  Lambda = map(samples, "Lambda"),
  Sigma = map(samples, "Sigma"),
  factors = map(samples, "factors"),
  Lambda_est = map(estimates, "lambda_est"),
  sigma_est = map(estimates, "sigma2_mean"),
  factors_est = map(estimates, "factors_est")
)
gof_inputs$mse <- pmap(gof_inputs, compute_goodness_of_fit)
gof_inputs$crps <- unlist(pmap(
  tibble(lambda_samps = map(post_draws, "Lambda_test"), true_val = gof_inputs$Lambda),
  avg_crps
))
gof_inputs$noise <- pmap(
  tibble(L = map(estimates, "lambda_est"), y = map(samples, "data")),
  contributed_variance_noise
)

sim_results <- gof_inputs |> mutate(
  est_heat_map = map(Lambda_est, make_heat_map),
  true_heat_map = map(Lambda, make_heat_map),
  trace_acf = trace_acf,
  post_samps = post_samps
)
saveRDS(sim_results, "results_sp2_Sparse_PLT.rds")
saveRDS(gof_inputs, "analysis_summary_sp2_Sparse_PLT.rds")
