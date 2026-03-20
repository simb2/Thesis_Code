# convergence_diagnostics ss2

results_uglt_setting_plt <- readRDS(here('results_uglt_setting_plt.rds'))
results_uglt_setting_UGLT <- readRDS(here("results_uglt_setting_UGLT.rds"))


geweke_results <- purrr::map(results_uglt_setting_plt$post_samps, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$draws$T_stat)
})
ESS_results <- purrr::map(results_uglt_setting_plt$post_samps, function(x) {
  LaplacesDemon::ESS(x$draws$T_stat)
})
lengths <- purrr::map(results_uglt_setting_plt$post_samps, function(x) {
  length(x$draws$T_stat)
})

geweke_results_UGLT <- purrr::map(results_uglt_setting_UGLT$post_samps, function(x){
  LaplacesDemon::Geweke.Diagnostic(x$draws$T_stat)
})

ESS_results_UGLT <- purrr::map(results_uglt_setting_UGLT$post_samps, function(x) {
  LaplacesDemon::ESS(x$draws$T_stat)
})

lengths_UGLT <- purrr::map(results_uglt_setting_UGLT$post_samps, function(x) {
  length(x$draws$T_stat)
})

library(tidyverse)

df_plt <- tibble(
  mod = "Sparse PLT", 
  V = results_uglt_setting_plt$V, 
  N = results_uglt_setting_plt$N, 
  samps = unlist(lengths),
  ESS = unlist(ESS_results), 
  Geweke = unlist(geweke_results)
) |> 
  arrange(V) |> 
  mutate(
    ESS_ratio = ESS/samps
  )
df_uglt <- tibble(
  mod = "Sparse PLT", 
  V = results_uglt_setting_UGLT$V, 
  N = results_uglt_setting_UGLT$N, 
  samps = unlist(lengths_UGLT),
  ESS = unlist(ESS_results_UGLT), 
  Geweke = unlist(geweke_results_UGLT)
) |> 
  arrange(V) |> 
  mutate(
    ESS_ratio = ESS/samps
  )

df_plt
df_uglt

