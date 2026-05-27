# convergence diagnostics: simulation study 3 (q = 10)
library(tidyverse)
library(here)

set.seed(1)

# run_sim_uglt_structure_overfitting.R saves one combined file per (data structure, model).
# Each file is a list of 2: [[1]] = q=4, [[2]] = q=10.
# Extract q=10 results (index 2) for this script.
uglt_draws_UGLT   <- readRDS(here("overfit_uglt_structure_post_draws_UGLT.rds"))
uglt_draws_sparse <- readRDS(here("overfit_uglt_structure_post_draws_sparse_PLT.rds"))
plt_draws_UGLT    <- readRDS(here("overfit_plt_structure_post_draws_UGLT.rds"))
plt_draws_sparse  <- readRDS(here("overfit_plt_structure_post_draws_sparse_PLT.rds"))

overfit_post_draws_10_UGLT                   <- uglt_draws_UGLT[2]
overfit_post_draws_plt_10_sparse             <- uglt_draws_sparse[2]
overfit_plt_setting_10_post_draws_UGLT       <- plt_draws_UGLT[2]
overfit_plt_setting_10_post_draws_plt_sparse <- plt_draws_sparse[2]

# ---- UGLT structure -----------------------------------------------------------
geweke_results <- purrr::map(overfit_post_draws_plt_10_sparse, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_results <- purrr::map(overfit_post_draws_plt_10_sparse, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengths <- purrr::map(overfit_post_draws_plt_10_sparse, function(x) {
  length(x$T_stat)
})

geweke_resultsUGLT <- purrr::map(overfit_post_draws_10_UGLT, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_resultsUGLT <- purrr::map(overfit_post_draws_10_UGLT, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengthsUGLT <- purrr::map(overfit_post_draws_10_UGLT, function(x) {
  length(x$T_stat)
})

df_plt <- tibble(
  mod = "UGLT setting, Sparse PLT fit",
  samps = unlist(lengths),
  ESS = unlist(ESS_results),
  Geweke = unlist(geweke_results)
) |> mutate(
  ESS_Ratio = ESS / samps
)
df_uglt <- tibble(
  mod = "UGLT Setting, UGLT fit",
  samps = unlist(lengthsUGLT),
  ESS = unlist(ESS_resultsUGLT),
  Geweke = unlist(geweke_resultsUGLT)
) |> mutate(
  ESS_Ratio = ESS / samps
)
df_plt
df_uglt

# ---- PLT structure ------------------------------------------------------------
geweke_results_plt_plt_sparse <- purrr::map(overfit_plt_setting_10_post_draws_plt_sparse, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_results_plt_plt_sparse <- purrr::map(overfit_plt_setting_10_post_draws_plt_sparse, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengths_plt_plt_sparse <- purrr::map(overfit_plt_setting_10_post_draws_plt_sparse, function(x) {
  length(x$T_stat)
})
geweke_results_plt_uglt <- purrr::map(overfit_plt_setting_10_post_draws_UGLT, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_results_plt_uglt <- purrr::map(overfit_plt_setting_10_post_draws_UGLT, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengths_plt_uglt <- purrr::map(overfit_plt_setting_10_post_draws_UGLT, function(x) {
  length(x$T_stat)
})

df_plt_plt_sparse <- tibble(
  mod = "PLT setting, sparse plt fit",
  samps = unlist(lengths_plt_plt_sparse),
  ESS = unlist(ESS_results_plt_plt_sparse),
  Geweke = unlist(geweke_results_plt_plt_sparse)
) |> mutate(
  ESS_Ratio = ESS / samps
)

df_plt_uglt <- tibble(
  mod = "PLT setting, UGLT fit",
  samps = unlist(lengths_plt_uglt),
  ESS = unlist(ESS_results_plt_uglt),
  Geweke = unlist(geweke_results_plt_uglt)
) |> mutate(
  ESS_Ratio = ESS / samps
)

df_plt
df_uglt
df_plt_plt_sparse
df_plt_uglt
