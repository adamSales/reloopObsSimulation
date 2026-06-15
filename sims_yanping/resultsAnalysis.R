library(dplyr)
library(tidyr)

# read data
setwd("/Users/ypei/Research/obsPLOOP/results")
yhatRF <- readRDS("yhatRF_sim_results_all.rds")
yhatX1RF <- readRDS("yhatX1RF_sim_results_all.rds")

yhatRF$spec <- "Yhat RF"
yhatX1RF$spec <- "Yhat + X1 RF"

sim_results <- bind_rows(
  yhatRF %>% mutate(spec = "Yhat RF"),
  yhatX1RF %>% mutate(spec = "Yhat + X1 RF")
)

table(sim_results$spec)
table(sim_results$scenario)
range(sim_results$sim)
length(unique(sim_results$sim))

# est
est_cols <- c(
  "simple_difference",
  "rebar_x1",
  "rebar_y0hat",
  "reloop_rf_po",
  "reloop_rf_v12",
  "reloop_rf_interp"
)

sim_long <- sim_results %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "estimator",
    values_to = "estimate"
  ) %>%
  mutate(error = estimate - true_tau)

mc_summary <- sim_long %>%
  group_by(spec, scenario, estimator) %>%
  summarise(
    n = n(),
    mean_estimate = mean(estimate),
    bias = mean(error),
    sd = sd(estimate),
    rmse = sqrt(mean(error^2)),
    mae = mean(abs(error)),
    mcse_bias = sd(error) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(scenario, spec, rmse)

mc_summary
write.csv(mc_summary, "mc_summary.csv", row.names = FALSE)
saveRDS(mc_summary, "mc_summary.rds")

# which one best under each condition
best_by_rmse <- mc_summary %>%
  group_by(spec, scenario) %>%
  slice_min(rmse, n = 6) %>%
  ungroup()

best_by_rmse
# print(best_by_rmse, n = Inf)

best_by_bias <- mc_summary %>%
  group_by(spec, scenario) %>%
  slice_min(abs(bias), n = 6) %>%
  ungroup()

print(best_by_bias, n = Inf)


# plot
library(ggplot2)

ggplot(mc_summary, aes(x = estimator, y = bias, color = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_wrap(~ scenario, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Bias",
    color = "Specification",
    title = "Monte Carlo Bias by Scenario"
  )
mc_summary %>% filter(!is.finite(bias))


ggplot(mc_summary, aes(x = estimator, y = rmse, color = spec)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_wrap(~ scenario, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "RMSE",
    color = "Specification",
    title = "Monte Carlo RMSE by Scenario"
  )

# dignostics
diag_summary <- sim_results %>%
  group_by(spec, scenario) %>%
  summarise(
    avg_treatment_rate = mean(treatment_rate_full),
    avg_n_pairs = mean(n_pairs),
    avg_abs_DeltaX1 = mean(mean_abs_DeltaX1),
    avg_abs_DeltaX2 = mean(mean_abs_DeltaX2),
    avg_abs_DeltaX3 = mean(mean_abs_DeltaX3),
    avg_abs_DeltaY0_hat = mean(mean_abs_DeltaY0_hat),
    avg_cor_Y0hat_Y0 = mean(cor_DeltaY0hat_DeltaY0),
    avg_cor_X1_Y0 = mean(cor_DeltaX1_DeltaY0),
    avg_cor_X2_Y0 = mean(cor_DeltaX2_DeltaY0),
    avg_cor_X3_Y0 = mean(cor_DeltaX3_DeltaY0),
    prop_unit1_treated = mean(prop_unit1_treated),
    .groups = "drop"
  )

diag_summary
write.csv(diag_summary, "diagnostic_summary.csv", row.names = FALSE)
