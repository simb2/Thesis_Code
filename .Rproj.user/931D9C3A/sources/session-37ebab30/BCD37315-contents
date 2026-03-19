df_builder <- function(obj, func) {
  get_diag <- function(obj, func) {
    diagnost <- purrr::map(obj$post_samps, function(x){
      if (!is.null(x$draws)) {
        func(x$draws$T_stat)
      } else {
        func(x$T_stat)
      }
    }) 
    return(diagnost)
  }
  geweke_results <- get_diag(obj, LaplacesDemon::Geweke.Diagnostic)
  ESS_results <- get_diag(obj, LaplacesDemon::ESS)
  lengths <- get_diag(obj, length)
  tibble(
    V = obj$V, 
    N = obj$N, 
    samps = unlist(lengths),
    ESS = unlist(ESS_results), 
    Geweke = unlist(geweke_results)
  ) |> 
    arrange(V) |> mutate(
      ESS_ratio = ESS/samps
    )
}

library(tidyverse)

df_plt_no_noise <- df_builder(results_no_nsp2) %>% 
  mutate(
    mod = 'plt no noise'
  )
df_plt_noise <- df_builder(results_nsp2) %>% 
  mutate(
    mod = 'plt noise'
  )
df_sparse_plt_no_noise <- df_builder(results_no_sp2_Sparse_PLT ) %>% 
  mutate(
    mod = 'sparse plt no noise'
  )
df_sparse_plt_noise <- df_builder(results_sp2_Sparse_PLT) %>% 
  mutate(
    mod = 'sparse plt noise'
  )

df_uglt_no_noise <- df_builder(results_no_sp2_UGLT) %>% 
  mutate(
    mod = 'uglt no noise'
  )
df_uglt_noise <- df_builder(results_sp2_UGLT) %>% 
  mutate(
    mod = 'uglt noise'
  )

df_plt_no_noise
df_plt_noise
df_sparse_plt_no_noise
df_sparse_plt_noise
df_uglt_noise
df_uglt_no_noise
