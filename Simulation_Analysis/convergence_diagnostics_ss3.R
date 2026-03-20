# convergence_diagnostics ss2

overfit_post_draws_UGLT <- readRDS(here("overfit_post_draws_UGLT.rds"))
overfit_post_draws_plt_sparse <- readRDS(here("overfit_post_draws_plt_sparse.rds"))
overfit_plt_setting_post_draws_plt_sparse <- readRDS(here("overfit_plt_setting_post_draws_plt_sparse.rds"))
overfit_plt_setting_post_draws_UGLT <- readRDS(here("overfit_plt_setting_post_draws_UGLT.rds"))
geweke_results <- purrr::map(overfit_post_draws_plt_sparse, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_results <- purrr::map(overfit_post_draws_plt_sparse, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengths <- purrr::map(overfit_post_draws_plt_sparse, function(x) {
  length(x$T_stat)
})

geweke_resultsUGLT <- purrr::map(overfit_post_draws_UGLT, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_resultsUGLT <- purrr::map(overfit_post_draws_UGLT, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengthsUGLT <- purrr::map(overfit_post_draws_UGLT, function(x) {
  length(x$T_stat)
})

library(tidyverse)

df_plt <- tibble(
  mod = "UGLT setting, Sparse PLT fit",
  samps = unlist(lengths),
  ESS = unlist(ESS_results), 
  Geweke = unlist(geweke_results)
) |> mutate(
  ESS_Ratio = ESS/samps
)
df_uglt <- tibble(
  mod = "UGLT Setting, UGLT fit", 
  samps = unlist(lengthsUGLT),
  ESS = unlist(ESS_resultsUGLT), 
  Geweke = unlist(geweke_resultsUGLT)
) |> mutate(
  ESS_Ratio = ESS/samps
)


geweke_results_plt_plt_sparse <- purrr::map(overfit_plt_setting_post_draws_plt_sparse, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_results_plt_plt_sparse <- purrr::map(overfit_plt_setting_post_draws_plt_sparse, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengths_plt_plt_sparse <- purrr::map(overfit_plt_setting_post_draws_plt_sparse, function(x) {
  length(x$T_stat)
})
geweke_results_plt_uglt <- purrr::map(overfit_plt_setting_post_draws_UGLT, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$T_stat)
})
ESS_results_plt_uglt <- purrr::map(overfit_plt_setting_post_draws_UGLT, function(x) {
  LaplacesDemon::ESS(x$T_stat)
})
lengths_plt_uglt <- purrr::map(overfit_plt_setting_post_draws_UGLT, function(x) {
  length(x$T_stat)
})


df_plt_plt_sparse <- tibble(
  mod = "PLT setting, sparse plt fit", 
  samps = unlist(lengths_plt_plt_sparse),
  ESS = unlist(ESS_results_plt_plt_sparse), 
  Geweke = unlist(geweke_results_plt_plt_sparse)
) |> mutate(
  ESS_Ratio = ESS/samps
)


df_plt_uglt <- tibble(
  mod = "PLT setting, UGLT fit", 
  samps = unlist(lengths_plt_uglt),
  ESS = unlist(ESS_results_plt_uglt), 
  Geweke = unlist(geweke_results_plt_uglt)
)|> mutate(
  ESS_Ratio = ESS/samps
)


df_plt
df_uglt
df_plt_plt_sparse
df_plt_uglt
