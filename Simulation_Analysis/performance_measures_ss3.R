# sim study 3 performance measures
library(tidyverse)
library(here)
overfit_post_draws_UGLT <- readRDS(here("overfit_post_draws_UGLT.rds"))
overfit_post_draws_plt_sparse <- readRDS(here("overfit_post_draws_plt_sparse.rds"))
table(overfit_post_draws_UGLT[[1]]$r)
table(overfit_post_draws_plt_sparse[[1]]$r)


overfit_plt_setting_post_draws_plt_sparse <- readRDS(here("overfit_plt_setting_post_draws_plt_sparse.rds"))
overfit_plt_setting_post_draws_UGLT <- readRDS(here("overfit_plt_setting_post_draws_UGLT.rds"))
table(overfit_plt_setting_post_draws_plt_sparse[[1]]$r)
table(overfit_plt_setting_post_draws_UGLT[[1]]$r)


overfit_post_draws_10_UGLT <- readRDS(here("overfit_post_draws_10_UGLT.rds"))
overfit_post_draws_plt_10_sparse <- readRDS(here("overfit_post_draws_plt_10_sparse.rds"))
table(overfit_post_draws_10_UGLT[[1]]$r)
table(overfit_post_draws_plt_10_sparse[[1]]$r)


overfit_plt_setting_10_post_draws_UGLT <- readRDS(here("overfit_plt_setting_10_post_draws_UGLT.rds"))
overfit_plt_setting_10_post_draws_plt_sparse <- readRDS(here("overfit_plt_setting_10_post_draws_plt_sparse.rds"))
table(overfit_plt_setting_10_post_draws_UGLT[[1]]$r)
table(overfit_plt_setting_10_post_draws_plt_sparse[[1]]$r)
