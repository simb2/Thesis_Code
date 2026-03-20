# performance measures 2
colours <- c('#332288', '#117733', '#44AA99', '#88CCEE',
             '#DDCC77', '#CC6677', '#AA4499', '#882255')
# Model results
library(tidyverse)
library(MASS)
results_uglt_setting_plt <- readRDS(here("results_uglt_setting_plt.rds"))
results_uglt_setting_UGLT <- readRDS(here("results_uglt_setting_UGLT.rds"))
source("Run_Simulations" ,'sim_data_3.R')

compute_goodness_of_fit <- function(N, V, q, Lambda, Sigma, factors, Lambda_est, sigma_est, factors_est) {
  runs <- 20
  mse_vec <- numeric(runs)
  
  for (i in 1:runs) {
    samp <- sim_data_3(N, Lambda, Sigma, factors = factors)
    samp2 <- sim_data_3(N, Lambda_est, Sigma = diag(sigma_est), factors = factors_est)
    mse_vec[i] <- mean((samp - samp2)^2)
  }
  return(list(mean(mse_vec), mse_vec))
}

compute_mse <- function(Lambda, Sigma, Lambda_est, sigma_est){
  true_cov <- Lambda %*% t(Lambda) + Sigma
  est_cov <- Lambda_est %*% t(Lambda_est) + diag(sigma_est)
  mse <- mean((true_cov - est_cov)^2)
  return(mse)
}

run_compute_mse <- function(x) {
  inputs <- tibble(
    Lambda = x$Lambda,
    Sigma = x$Sigma,
    Lambda_est = x$Lambda_est,
    sigma_est = x$sigma_est
  )
  return(pmap(inputs, compute_mse))
}

run_gof <- function(x) {
  gof_inputs <- tibble(
    N = x$N, 
    V = x$V, 
    q = x$q,
    Lambda = x$Lambda,
    Sigma = x$Sigma,
    factors = x$factors,
    Lambda_est = x$Lambda_est,
    sigma_est = x$sigma_est,
    factors_est = x$factors_est
  )
  return(pmap(gof_inputs, compute_goodness_of_fit))
}

results_uglt_setting_plt$gof <- run_gof(results_uglt_setting_plt)
results_uglt_setting_UGLT$gof <- run_gof(results_uglt_setting_UGLT)


results_uglt_setting_plt$mse <- run_compute_mse(results_uglt_setting_plt)
results_uglt_setting_UGLT$mse <- run_compute_mse(results_uglt_setting_UGLT)


GOF1 <- unlist(purrr::map(results_uglt_setting_plt$gof, function(x) {
  x[[1]]
}))
GOF2 <- unlist(purrr::map(results_uglt_setting_UGLT$gof, function(x){
  x[[1]]
}))
MSE1 <- unlist(purrr::map(results_uglt_setting_plt$mse, function(x) {
  x[[1]]
}))

MSE2 <- unlist(purrr::map(results_uglt_setting_UGLT$mse, function(x) {
  x[[1]]
}))

tru_noise <- unlist(results_uglt_setting_UGLT$tru_noise)

df_plt <- tibble(
  mod = "Sparse PLT", 
  V = results_uglt_setting_plt$V, 
  N = results_uglt_setting_plt$N, 
  MSE = MSE1, 
  gof = GOF1,
  crps = unlist(results_uglt_setting_plt$crps), 
  noise = abs(unlist(results_uglt_setting_plt$noise) - tru_noise)
) |> 
  arrange(V)

df_uglt <- tibble(
  mod = "UGLT", 
  V = results_uglt_setting_UGLT$V, 
  N = results_uglt_setting_UGLT$N, 
  MSE = MSE2, 
  gof = GOF2,
  crps = results_uglt_setting_UGLT$crps, 
  noise = abs(unlist(results_uglt_setting_UGLT$noise) - tru_noise) 
) |> 
  arrange(V)

df_plt
df_uglt


df_baseline <- bind_rows(
  df_uglt,
  df_plt
) %>%
  pivot_longer(cols = c(MSE, crps, noise, gof), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "crps" = "CRPS",
                         "MSE" = "MSE", 
                         "noise" = "Noise Ratio", 
                         "gof" = "Goodness of Fit"
  ))

p1 <- ggplot(df_baseline, aes(x = V, y = value, fill = mod)) +
  geom_col(position = position_dodge(width = 3.1), width = 3) + 
  facet_grid(metric ~ ., scales = "free_y") +
  theme_minimal() +
  labs(
    fill = "Model",
    x = "Number of Variables (V)",
    y = ""
  ) + 
  theme(legend.position = "bottom") +
  scale_fill_manual(values = colours) +
  scale_x_continuous(breaks = c(20, 30, 40))
p1
