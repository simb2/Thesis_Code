# Make trace plots: simulation study 3 (q = 10, PLT structure)
library(tidyverse)
library(here)

set.seed(1)

# run_sim_uglt_structure_overfitting.R saves one combined file per (data structure, model).
# Each file is a list of 2: [[1]] = q=4, [[2]] = q=10.
# Extract q=10 results (index 2) for these plots.
overfit_post_draws_10_UGLT                   <- readRDS(here("overfit_uglt_structure_post_draws_UGLT.rds"))[2]
overfit_plt_setting_10_post_draws_UGLT       <- readRDS(here("overfit_plt_structure_post_draws_UGLT.rds"))[2]
overfit_post_draws_plt_10_sparse             <- readRDS(here("overfit_uglt_structure_post_draws_sparse_PLT.rds"))[2]
overfit_plt_setting_10_post_draws_plt_sparse <- readRDS(here("overfit_plt_structure_post_draws_sparse_PLT.rds"))[2]

make_trace_acf_plot <- function(post_samp) {
  T_stat_df <- tibble(
    Value     = post_samp$T_stat,
    Iteration = seq_along(post_samp$T_stat)
  )
  acf   <- acf(post_samp$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  list(acf = acf, trace = trace)
}

plot(acf(overfit_plt_setting_10_post_draws_plt_sparse[[1]]$T_stat), main = '')

T_stat_df <- tibble(
  Value     = overfit_plt_setting_10_post_draws_plt_sparse[[1]]$T_stat,
  Iteration = seq_along(overfit_plt_setting_10_post_draws_plt_sparse[[1]]$T_stat)
)

trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
  geom_line(linewidth = 0.4, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Iteration", y = "Value")

trace
