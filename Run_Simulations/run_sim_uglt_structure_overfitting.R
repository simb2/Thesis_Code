library(tidyverse)
library(MASS)
library(scoringRules)
library(here)

source(here("Run_Simulations", "sim_data_3.R"))
source(here("Run_Simulations", "sim_helpers.R"))
source(here("MCMC_Algorithms", "run_mcmc_UGLT.R"))
source(here("MCMC_Algorithms", "run_mcmc_sparse_PLT.R"))

set.seed(8)
n_subj <- 2000
n_vars <- 30
n_factors <- c(4, 10)

# ---- Data: UGLT structure (pivots sampled uniformly from {1,...,v-3}) -----
settings_uglt <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors)
samples_uglt <- purrr::pmap(settings_uglt, sim_data_UGLT, plt_structure = FALSE)
starting_vals_uglt <- purrr::map(samples_uglt, get_starting_vals_uglt, n_runs = 6000, q_overfit = 3)

post_samps_uglt_uglt <- purrr::map(starting_vals_uglt, ~ do.call(run_mcmc_UGLT, .x))
post_draws_uglt_uglt <- purrr::map(post_samps_uglt_uglt, "draws")

post_samps_plt_uglt <- purrr::map(starting_vals_uglt, ~ do.call(run_mcmc_sparse_PLT, .x))
post_draws_plt_uglt <- purrr::map(post_samps_plt_uglt, "draws")

saveRDS(post_draws_uglt_uglt, "overfit_uglt_structure_post_draws_UGLT.rds")
saveRDS(post_draws_plt_uglt, "overfit_uglt_structure_post_draws_sparse_PLT.rds")
saveRDS(samples_uglt, "overfit_uglt_structure_init_val.rds")

# ---- Data: PLT structure (pivots fixed at l_j = j) -----------------------
settings_plt <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors)
samples_plt <- purrr::pmap(settings_plt, sim_data_UGLT, plt_structure = TRUE)
starting_vals_plt <- purrr::map(samples_plt, get_starting_vals_uglt, n_runs = 6000, q_overfit = 3)

post_samps_uglt_plt <- purrr::map(starting_vals_plt, ~ do.call(run_mcmc_UGLT, .x))
post_draws_uglt_plt <- purrr::map(post_samps_uglt_plt, "draws")

post_samps_plt_plt <- purrr::map(starting_vals_plt, ~ do.call(run_mcmc_sparse_PLT, .x))
post_draws_plt_plt <- purrr::map(post_samps_plt_plt, "draws")

saveRDS(post_draws_uglt_plt, "overfit_plt_structure_post_draws_UGLT.rds")
saveRDS(post_draws_plt_plt, "overfit_plt_structure_post_draws_sparse_PLT.rds")
saveRDS(samples_plt, "overfit_plt_structure_init_val.rds")
