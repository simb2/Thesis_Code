library(tidyverse)
library(MASS)
library(scoringRules)
library(here)

source(here("Run_Simulations", "sim_data_3.R"))
source(here("Run_Simulations", "sim_helpers.R"))

set.seed(20)
n_subj <- 2400
n_vars <- c(20, 30, 40)
n_vars <- c(20)
n_factors <- 3

settings <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors)
samples <- pmap(settings, sim_data_UGLT)
starting_vals <- purrr::map(samples, get_starting_vals_uglt)

# ---- UGLT ----------------------------------------------------------------
source(here("MCMC_Algorithms", "run_mcmc_UGLT_Static.R"))
post_samps <- profvis::profvis(purrr::map(starting_vals, ~ do.call(run_mcmc_UGLT, .x)))
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

tru_noise <- pmap(
  tibble(L = map(samples, "Lambda"), y = map(samples, "data")),
  contributed_variance_noise
)

sim_results <- gof_inputs |> mutate(
  est_heat_map = map(Lambda_est, make_heat_map),
  true_heat_map = map(Lambda, make_heat_map),
  trace_acf = trace_acf,
  post_samps = post_samps,
  tru_noise = tru_noise
)

saveRDS(sim_results, "results_uglt_setting_UGLT.rds")
saveRDS(gof_inputs, "analysis_uglt_setting_UGLT.rds")

# ---- Sparse PLT (same data) -----------------------------------------------
source(here("MCMC_Algorithms", "run_mcmc_sparse_PLT_Static.R"))
post_samps <- profvis::profvis(purrr::map(starting_vals, ~ do.call(run_mcmc_sparse_PLT, .x)))
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

saveRDS(sim_results, "results_uglt_setting_plt.rds")
saveRDS(gof_inputs, "analysis__uglt_setting_plt.rds")
