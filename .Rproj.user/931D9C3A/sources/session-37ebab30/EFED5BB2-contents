# make trace plot
overfit_post_draws_10_UGLT <- readRDS("~/Projects/Thesis/masters_thesis/BayesianFactorMod/overfit_post_draws_10_UGLT.rds")
overfit_plt_setting_10_post_draws_UGLT <- readRDS("~/Projects/Thesis/masters_thesis/BayesianFactorMod/overfit_plt_setting_10_post_draws_UGLT.rds")
overfit_post_draws_plt_10_sparse <- readRDS("~/Projects/Thesis/masters_thesis/BayesianFactorMod/overfit_post_draws_plt_10_sparse.rds")
overfit_plt_setting_10_post_draws_plt_sparse <- readRDS("~/Projects/Thesis/masters_thesis/BayesianFactorMod/overfit_plt_setting_10_post_draws_plt_sparse.rds")

make_trace_acf_plot <- function(post_samp) {
  T_stat_df <- tibble(
    Value = post_samp$T_stat,
    Iteration = 1:length(post_samp$T_stat)
  )
  acf <- acf(post_samp$T_stat)
  trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Iteration", y = "Value")
  return(list(acf = acf, trace = trace))
}

plot(acf(overfit_plt_setting_10_post_draws_plt_sparse[[1]]$T_stat), main = '')


T_stat_df <- tibble(
  Value = overfit_plt_setting_10_post_draws_plt_sparse[[1]]$T_stat,
  Iteration = seq_along(overfit_plt_setting_10_post_draws_plt_sparse[[1]]$T_stat)
)

trace <- ggplot(T_stat_df, mapping = aes(x = Iteration, y = Value)) +
  geom_line(linewidth = 0.4, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Iteration", y = "Value")


trace