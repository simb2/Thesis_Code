# sim study 3 performance measures
library(tidyverse)
library(here)

set.seed(1)

# run_sim_uglt_structure_overfitting.R saves one combined file per (data structure, model)
# Each file is a list of 2: [[1]] = q=4, [[2]] = q=10.
uglt_draws_UGLT   <- readRDS(here("overfit_uglt_structure_post_draws_UGLT.rds"))
uglt_draws_sparse <- readRDS(here("overfit_uglt_structure_post_draws_sparse_PLT.rds"))
plt_draws_UGLT    <- readRDS(here("overfit_plt_structure_post_draws_UGLT.rds"))
plt_draws_sparse  <- readRDS(here("overfit_plt_structure_post_draws_sparse_PLT.rds"))

# ---- q = 4 (index 1) ----------------------------------------------------------
overfit_post_draws_UGLT                   <- uglt_draws_UGLT[1]
overfit_post_draws_plt_sparse             <- uglt_draws_sparse[1]
overfit_plt_setting_post_draws_UGLT       <- plt_draws_UGLT[1]
overfit_plt_setting_post_draws_plt_sparse <- plt_draws_sparse[1]

table(overfit_post_draws_UGLT[[1]]$r)
table(overfit_post_draws_plt_sparse[[1]]$r)
table(overfit_plt_setting_post_draws_plt_sparse[[1]]$r)
table(overfit_plt_setting_post_draws_UGLT[[1]]$r)

# ---- q = 10 (index 2) ---------------------------------------------------------
overfit_post_draws_10_UGLT                   <- uglt_draws_UGLT[2]
overfit_post_draws_plt_10_sparse             <- uglt_draws_sparse[2]
overfit_plt_setting_10_post_draws_UGLT       <- plt_draws_UGLT[2]
overfit_plt_setting_10_post_draws_plt_sparse <- plt_draws_sparse[2]

table(overfit_post_draws_10_UGLT[[1]]$r)
table(overfit_post_draws_plt_10_sparse[[1]]$r)
table(overfit_plt_setting_10_post_draws_UGLT[[1]]$r)
table(overfit_plt_setting_10_post_draws_plt_sparse[[1]]$r)
