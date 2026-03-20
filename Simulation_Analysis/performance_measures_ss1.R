# Model results
colours <- c('#332288', '#117733', '#44AA99', '#88CCEE',
             '#DDCC77', '#CC6677', '#AA4499', '#882255')
options(ggplot2.discrete.colour= colours)


# I'm gonna have to read over this part again very very carefully!!!!
library(tidyverse)
library(here)
library(MASS)
results_nsp2 <- readRDS(here('results_nsp2.rds'))
results_no_nsp2 <- readRDS(here('results_no_nsp2.rds'))
results_no_sp2_Sparse_PLT  <- readRDS(here('results_no_sp2_Sparse_PLT.rds'))
results_no_sp2_UGLT  <- readRDS(here('results_no_sp2_UGLT.rds'))
results_sp2_Sparse_PLT <- readRDS(here('results_sp2_Sparse_PLT.rds'))
results_sp2_UGLT <- readRDS(here('results_sp2_UGLT.rds'))

source(here("Run_Simulations", 'sim_study_noisy_variables_helpers.R'))
source(here("Run_Simulations", 'sim_data_3.R'))
# creating a function to recompute mse. 

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


results_nsp2$mse <- run_compute_mse(results_nsp2)
results_no_nsp2$mse <- run_compute_mse(results_no_nsp2)
results_sp2_Sparse_PLT$mse <- run_compute_mse(results_sp2_Sparse_PLT)
results_no_sp2_Sparse_PLT $mse <- run_compute_mse(results_no_sp2_Sparse_PLT )
results_sp2_UGLT$mse <- run_compute_mse(results_sp2_UGLT)
results_no_sp2_UGLT $mse <- run_compute_mse(results_no_sp2_UGLT )


results_nsp2$gof <- run_gof(results_nsp2)
results_no_nsp2$gof <- run_gof(results_no_nsp2)
results_sp2_Sparse_PLT$gof <- run_gof(results_sp2_Sparse_PLT)
results_no_sp2_Sparse_PLT $gof <- run_gof(results_no_sp2_Sparse_PLT )
results_sp2_UGLT$gof <- run_gof(results_sp2_UGLT)
results_no_sp2_UGLT $gof <- run_gof(results_no_sp2_UGLT )

GOF1 <- unlist(purrr::map(results_nsp2$gof, function(x) {
  x[[1]]
}))

GOF2 <- unlist(purrr::map(results_no_nsp2$gof, function(x) {
  x[[1]]
}))

GOF3 <- unlist(purrr::map(results_sp2_Sparse_PLT$gof, function(x) {
  x[[1]]
}))

GOF4 <- unlist(purrr::map(results_no_sp2_Sparse_PLT $gof, function(x) {
  x[[1]]
}))

GOF5 <- unlist(purrr::map(results_sp2_UGLT$gof, function(x) {
  x[[1]]
}))

GOF6 <- unlist(purrr::map(results_no_sp2_UGLT $gof, function(x) {
  x[[1]]
}))

extract_mse <- function(results) unlist(purrr::map(results$mse, ~ .x[[1]]))

MSE1 <- extract_mse(results_nsp2)
MSE2 <- extract_mse(results_no_nsp2)
MSE3 <- extract_mse(results_sp2_Sparse_PLT)
MSE4 <- extract_mse(results_no_sp2_Sparse_PLT )
MSE5 <- extract_mse(results_sp2_UGLT)   # was MSE4 (duplicate) in your code
MSE6 <- extract_mse(results_no_sp2_UGLT )


crps <- unlist(results_nsp2$crps)
noise <- unlist(results_nsp2$noise)
tru_noise <- unlist(results_nsp2$tru_noise)

# First two (special cases)
df_plt_noise <- tibble(
  mod = "Non-sparse PLT",
  V = results_nsp2$V,
  N = results_nsp2$N,
  gof = GOF1,
  mse = MSE1, 
  crps = crps,
  noise = abs(noise - tru_noise),
  tru_noise = tru_noise
) |> arrange(V)

df_plt_no_noise <- tibble(
  mod = "Non-sparse PLT",
  V = results_no_nsp2$V,
  N = results_no_nsp2$N,
  gof = GOF2,
  mse = MSE2, 
  crps = unlist(results_no_nsp2$crps),
  noise = abs(unlist(results_no_nsp2$noise) - unlist(results_no_nsp2$tru_noise)),
  tru_noise = unlist(results_no_nsp2$tru_noise)
) |> arrange(V)

# Last four (uniform structure)
make_df <- function(mod_name, results, gof) {
  tibble(
    mod = mod_name,
    V = results$V,
    N = results$N,
    gof = gof,
    crps = results$crps
  ) |> arrange(V)
}

df_sparse_plt_noise <- make_df("Sparse PLT", results_sp2_Sparse_PLT, GOF3) |> mutate(
  mse = MSE3,
  noise = abs(unlist(results_nsp2$tru_noise)[-1] - unlist(results_sp2_Sparse_PLT$noise))
)
df_sparse_plt_no_noise <- make_df("Sparse PLT", results_no_sp2_Sparse_PLT , GOF4) |> mutate(
  mse = MSE4,
  noise = abs(unlist(results_no_nsp2$tru_noise) - unlist(results_no_sp2_Sparse_PLT$noise))
)
df_uglt_noise <- make_df("UGLT", results_sp2_UGLT, GOF5) |> mutate(
  mse = MSE5,
  noise = abs(unlist(results_nsp2$tru_noise)[-1] - unlist(results_sp2_UGLT$noise))
)
df_uglt_no_noise <- make_df("UGLT", results_no_sp2_UGLT , GOF6) |> mutate(
  mse = MSE6,
  noise = abs(unlist(results_no_nsp2$tru_noise) - unlist(results_no_sp2_UGLT$noise))
)



df_plt_noise
df_plt_no_noise
df_sparse_plt_noise
df_sparse_plt_no_noise
df_uglt_noise
df_uglt_no_noise

df_all <- bind_rows(
  df_plt_noise %>% mutate(noisy = "Noisy"),
  df_plt_no_noise %>% mutate(noisy = "No Noise"),
  df_sparse_plt_noise %>% mutate(noisy = "Noisy"),
  df_sparse_plt_no_noise %>% mutate(noisy = "No Noise"),
  df_uglt_noise %>% mutate(noisy = "Noisy"),
  df_uglt_no_noise %>% mutate(noisy = "No Noise")
) %>%
  pivot_longer(cols = c(gof, crps, noise), names_to = "metric", values_to = "value")

ggplot(df_all, aes(x = N, y = value, color = noisy, shape = mod)) +
  geom_line() +
  geom_point(size = 2, position = position_dodge(width = 1), alpha = 0.7) +
  facet_grid(metric ~ V, scales = "free_y", labeller = labeller(V = label_both)) +
  theme_minimal() +
  labs(
    color = "Noise Conditio",
    shape = "Model",
    x = "Number of Observations (N)",
    y = ""
  ) +
  scale_color_manual(values = colours)


df_baseline <- bind_rows(
  df_plt_no_noise,
  df_sparse_plt_no_noise,
  df_uglt_no_noise
) %>%
  pivot_longer(cols = c(mse, crps, noise, gof), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "crps" = "CRPS",
                         "mse" = "MSE", 
                         "gof" = "Goodness of Fit",
                         "noise" = "Noise Ratio"
  ))

p1 = ggplot(df_baseline, aes(x = factor(N), y = value, fill = mod)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(metric ~ V, scales = "free_y", labeller = labeller(V = label_both)) +
  theme_minimal() +
  labs(
    fill = "Model",
    x = "Number of Observations (N)",
    y = ""
  ) + theme(legend.position = "none") +
  scale_fill_manual(values = colours)



df_noise <- bind_rows(
  df_plt_noise,
  df_sparse_plt_noise,
  df_uglt_noise
) %>%
  pivot_longer(cols = c(mse, gof, crps, noise), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "crps" = "CRPS",
                         "mse" = "MSE", 
                         "noise" = "Noise Ratio", 
                         "gof" = "Goodness of Fit"
  ))

p2 = ggplot(df_noise, aes(x = factor(N), y = value, fill = mod)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(metric ~ V, scales = "free_y", labeller = labeller(V = label_both)) +
  theme_minimal() +
  labs(
    fill = "Model",
    x = "Number of Observations (N)",
    y = ""
  ) + theme(legend.position = "none") +
  scale_fill_manual(values = colours)



library(patchwork)
p1
p2
p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#####################################################################
# Heat map:


library(ggplot2)
library(reshape2)

# Your lambda matrix (rows = i, cols = j)
lambda <- results_sp2_UGLT$Lambda[[3]]

# Convert to long format
df <- reshape2::melt(lambda)
df2 <- reshape2::melt(results_sp2_UGLT$Lambda_est[[3]])
colnames(df) <- c("i", "j", "value")
colnames(df2) <- c("i", "j", "value")
p3 <- ggplot(df, aes(x=j, y=i, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(
    low="blue", mid="white", high="red",
    midpoint=0
  ) +
  scale_y_reverse() +   # flip y so row 1 is at top
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_x_continuous(breaks = 1:2)



p4 <- ggplot(df2, aes(x=j, y=i, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(
    low="blue", mid="white", high="red",
    midpoint=0
  ) +
  scale_y_reverse() +   # flip y so row 1 is at top
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_x_continuous(breaks = 1:2)

p3 + p4
