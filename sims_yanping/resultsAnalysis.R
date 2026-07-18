library(dplyr)
library(tidyr)

setwd("/Users/ypei/Research/obsPLOOP/results")
sim_results<- readRDS("yhatRF_Sim_results_all.rds")

# estimator columns
est_cols <- c(
  "simple_difference",
  # "rebar_x1",
  "rebar_y0hat",
  "reloop_rf_po",
  "reloop_rf_v12",
  "reloop_rf_interp",
  "cf_loop_p05",
  "cf_ipw_rf",
  "cf_aipw_rf"
)

# long format
sim_long <- sim_results %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "estimator",
    values_to = "estimate"
  ) %>%
  mutate(error = estimate - true_tau)

# Monte Carlo summary
mc_summary <- sim_long %>%
  group_by(scenario, functional_form, matching_quality, estimator) %>%
  summarise(
    n = sum(is.finite(estimate)),
    mean_estimate = mean(estimate, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae = mean(abs(error), na.rm = TRUE),
    mcse_bias = sd(error, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(scenario, functional_form, matching_quality, rmse)

saveRDS(mc_summary, "Sim_mc_summary.rds")
write.csv(mc_summary, "Sim_mc_summary.csv", row.names = FALSE)

# best estimators by RMSE and absolute bias
best_by_rmse <- mc_summary %>%
  filter(is.finite(rmse)) %>%
  group_by(scenario, functional_form, matching_quality) %>%
  slice_min(rmse, n = 3) %>%
  ungroup()

best_by_bias <- mc_summary %>%
  filter(is.finite(bias)) %>%
  group_by(scenario, functional_form, matching_quality) %>%
  slice_min(abs(bias), n = 9) %>%
  ungroup()

saveRDS(best_by_rmse, "Sim_best_by_rmse.rds")
saveRDS(best_by_bias, "Sim_best_by_bias.rds")
write.csv(best_by_rmse, "Sim_best_by_rmse.csv", row.names = FALSE)
write.csv(best_by_bias, "Sim_best_by_bias.csv", row.names = FALSE)

# diagnostics
diag_summary <- sim_results %>%
  group_by(scenario, functional_form, matching_quality) %>%
  summarise(
    n_sims = n_distinct(sim),
    avg_treatment_rate = mean(treatment_rate_full),
    avg_n_pairs = mean(n_pairs),
    avg_mean_DeltaX1 = mean(mean_DeltaX1),
    avg_mean_DeltaX2 = mean(mean_DeltaX2),
    avg_mean_DeltaX3 = mean(mean_DeltaX3),
    avg_abs_DeltaX1 = mean(mean_abs_DeltaX1),
    avg_abs_DeltaX2 = mean(mean_abs_DeltaX2),
    avg_abs_DeltaX3 = mean(mean_abs_DeltaX3),
    avg_abs_DeltaY0_hat = mean(mean_abs_DeltaY0_hat),
    avg_cor_Y0hat_Y0 = mean(cor_DeltaY0hat_DeltaY0, na.rm = TRUE),
    avg_cor_X1_Y0 = mean(cor_DeltaX1_DeltaY0, na.rm = TRUE),
    avg_cor_X2_Y0 = mean(cor_DeltaX2_DeltaY0, na.rm = TRUE),
    avg_cor_X3_Y0 = mean(cor_DeltaX3_DeltaY0, na.rm = TRUE),
    avg_e_hat_cf = mean(mean_e_hat_cf, na.rm = TRUE),
    avg_sd_e_hat_cf = mean(sd_e_hat_cf, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(diag_summary, "Sim_diagnostic_summary.rds")
write.csv(diag_summary, "Sim_diagnostic_summary.csv", row.names = FALSE)

print(mc_summary, n = Inf)
print(best_by_rmse, n = Inf)
print(best_by_bias, n = Inf)
# print(diag_summary, n = Inf)
print(diag_summary, n = Inf, width = Inf)



# `cf_ipw_rf` and `cf_aipw_rf` perform poorly:
## `cf_ipw_rf` and `cf_aipw_rf` have large negative bias in many settings, especially `cf_ipw_rf`.
## This suggests that using the full-data RF propensity score as a post-matching IPW weight may over-correct the estimates or target a different population than the matched sample.
## This is useful to know. IPW/AIPW is not automatically well-suited to this post-matching paired setting, especially after the matched sample has already been forced to have a 1:1 treated-control structure.
## The RF propensity score still reflects the original observational treatment assignment mechanism, while the final analysis is conducted only on the matched sample.

# `cf_loop_p05` often performs well:
## This is also reasonable. In the matched pairs, the treated-control ratio is fixed at 1:1.
## Therefore, a LOOP-style estimator with fixed treatment probability `p = 0.5` is more closely aligned with the matched design.
## In other words, `cf_loop_p05` respects the post-matching structure better than the IPW estimators that use a variable estimated propensity score.

# `reloop_rf_po` is usually more stable than `reloop_rf_v12` and `reloop_rf_interp`:
## Across many settings, `reloop_rf_po` ranks near the top. This suggests that the individuals-as-units RF imputation strategy works better in this simulation design.
## By contrast, the pair-as-units version, `reloop_rf_v12`, and the interpolated version, `reloop_rf_interp`,
## appear less stable or more biased in several settings. This may mean that, in this data-generating process, 
## learning individual-level potential outcomes is easier than learning pair-level contrasts directly.


library(dplyr)
library(tidyr)
library(ggplot2)

sim_results <- readRDS("yhatX1RF_Sim_results_all.rds")
mc_summary <- readRDS("Sim_mc_summary.rds")
diag_summary <- readRDS("Sim_diagnostic_summary.rds")

# common setup
est_order <- c(
  "simple_difference",
  # "rebar_x1",
  "rebar_y0hat",
  "reloop_rf_po",
  "reloop_rf_v12",
  "reloop_rf_interp",
  "cf_loop_p05",
  "cf_ipw_rf",
  "cf_aipw_rf"
)

est_labels <- c(
  simple_difference = "Simple diff",
  # rebar_x1 = "Rebar X1",
  rebar_y0hat = "Rebar Yhat",
  reloop_rf_po = "ReLOOP RF PO",
  reloop_rf_v12 = "ReLOOP RF V12",
  reloop_rf_interp = "ReLOOP RF interp",
  cf_loop_p05 = "CF LOOP p=0.5",
  cf_ipw_rf = "CF IPW RF",
  cf_aipw_rf = "CF AIPW RF"
)

# q1: Diagnostics balance plot
# Does good matching actually reduce residual covariate imbalance compared to bad matching?
p_balance <- diag_summary %>%
  select(
    scenario, functional_form, matching_quality,
    avg_mean_DeltaX1, avg_mean_DeltaX2, avg_mean_DeltaX3
  ) %>%
  pivot_longer(
    cols = starts_with("avg_mean_Delta"),
    names_to = "covariate",
    values_to = "imbalance"
  ) %>%
  mutate(
    covariate = recode(
      covariate,
      avg_mean_DeltaX1 = "Delta X1",
      avg_mean_DeltaX2 = "Delta X2",
      avg_mean_DeltaX3 = "Delta X3"
    ),
    matching_quality = factor(matching_quality, levels = c("bad", "good")),
    functional_form = factor(functional_form, levels = c("linear", "nonlinear"))
  ) %>%
  ggplot(aes(x = matching_quality, y = imbalance, fill = matching_quality)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(scenario + functional_form ~ covariate, scales = "free_y") +
  scale_fill_manual(values = c(bad = "#C44E52", good = "#4C72B0")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Average signed imbalance",
    fill = "Matching quality",
    title = "Covariate Imbalance After Matching"
  )

p_balance
ggsave("Sim_plot_1_balance_diagnostics.png", p_balance, width = 12, height = 8, dpi = 300)

# q2: Signed bias heatmap
# Do different estimators exhibit positive or negative bias under different scenarios, functional forms, and matching qualities? 
# What are the overall bias patterns?
p_bias_heatmap <- mc_summary %>%
  mutate(
    estimator = factor(estimator, levels = est_order, labels = est_labels[est_order]),
    setting = paste(functional_form, matching_quality, sep = " / "),
    setting = factor(
      setting,
      levels = c("linear / bad", "linear / good", "nonlinear / bad", "nonlinear / good")
    )
  ) %>%
  ggplot(aes(x = estimator, y = setting, fill = bias)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(bias, 2)), size = 3) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Bias",
    title = "Signed Bias by Estimator and Setting"
  )

p_bias_heatmap
ggsave("Sim_plot_2_signed_bias_heatmap.png", p_bias_heatmap, width = 13, height = 8, dpi = 300)

# q3: RMSE heatmap
# Which estimator achieves the best bias–variance tradeoff overall?
p_rmse_heatmap <- mc_summary %>%
  mutate(
    estimator = factor(estimator, levels = est_order, labels = est_labels[est_order]),
    setting = paste(functional_form, matching_quality, sep = " / "),
    setting = factor(
      setting,
      levels = c("linear / bad", "linear / good", "nonlinear / bad", "nonlinear / good")
    )
  ) %>%
  ggplot(aes(x = estimator, y = setting, fill = rmse)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(rmse, 2)), size = 3) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "RMSE",
    title = "RMSE by Estimator and Setting"
  )

p_rmse_heatmap
ggsave("Sim_plot_3_rmse_heatmap.png", p_rmse_heatmap, width = 13, height = 8, dpi = 300)

# q4:Top 3 RMSE table
top3_rmse <- mc_summary %>%
  filter(is.finite(rmse)) %>%
  group_by(scenario, functional_form, matching_quality) %>%
  slice_min(rmse, n = 3) %>%
  ungroup() %>%
  arrange(scenario, functional_form, matching_quality, rmse) %>%
  mutate(
    estimator_label = recode(estimator, !!!est_labels)
  ) %>%
  select(
    scenario, functional_form, matching_quality,
    estimator = estimator_label,
    bias, sd, rmse, mae
  )

print(top3_rmse, n = Inf)
write.csv(top3_rmse, "Sim_table_3_top3_by_rmse.csv", row.names = FALSE)

# q5: Focus on IPW / AIPW instability
# Do CF-IPW-RF and CF-AIPW-RF exhibit unstable performance in the post-matching sample?
# How do they compare with CF-Loop-P05 and ReLoop-RF-PO?
focus_estimators <- c(
  "cf_ipw_rf",
  "cf_aipw_rf",
  "cf_loop_p05",
  "reloop_rf_po"
)

p_focus <- mc_summary %>%
  filter(estimator %in% focus_estimators) %>%
  mutate(
    estimator = factor(estimator, levels = focus_estimators, labels = est_labels[focus_estimators]),
    setting = paste(functional_form, matching_quality, sep = " / "),
    setting = factor(
      setting,
      levels = c("linear / bad", "linear / good", "nonlinear / bad", "nonlinear / good")
    )
  ) %>%
  ggplot(aes(x = setting, y = bias, color = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ scenario, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    x = NULL,
    y = "Bias",
    color = "Estimator",
    title = "Focused Comparison: IPW/AIPW vs LOOP-style and ReLOOP"
  )

p_focus
ggsave("Sim_plot_4_focus_ipw_aipw.png", p_focus, width = 11, height = 8, dpi = 300)

#####
compact_table <- mc_summary %>%
  filter(estimator %in% c(
    "simple_difference",
    "rebar_y0hat",
    "reloop_rf_po",
    "cf_loop_p05",
    "cf_ipw_rf",
    "cf_aipw_rf"
  )) %>%
  mutate(
    estimator = recode(estimator, !!!est_labels)
  ) %>%
  select(
    scenario, functional_form, matching_quality,
    estimator, bias, rmse
  ) %>%
  arrange(scenario, functional_form, matching_quality, estimator)

print(compact_table, n = Inf)
write.csv(compact_table, "Sim_table_4_compact_key_estimators.csv", row.names = FALSE)


# 1. Good matching does reduce imbalance, but not completely
## The balance plot shows that good matching substantially reduces residual imbalance in `X1`, especially compared with bad matching.
## This is exactly what we hoped to see, because bad matching only matches on `X3`, while good matching also includes `X1` and, in the nonlinear setting, the nonlinear assignment terms.
## However, good matching does not eliminate imbalance completely. In the nonlinear settings, there is still noticeable residual imbalance in `X1` and `X3`. So nonlinear settings remain harder even after good matching.

# 2. Nonlinear settings are much harder
## The signed bias and RMSE heatmaps show that estimator performance is much worse in nonlinear settings, especially in the aligned scenario.
## For example, in the aligned nonlinear settings, the simple difference and `rebar_yhat` still have large positive bias.
## This makes sense because treatment assignment and outcomes share nonlinear structure involving `X1`, so residual imbalance after matching translates into large baseline outcome imbalance.

# 3. `cf_loop_p05` is the most stable overall
## Across many settings, `cf_loop_p05` has small bias and low RMSE. This is especially clear in the RMSE heatmap.
## This is reasonable because the matched sample has a fixed 1:1 treated-control structure.
## A LOOP-style estimator with fixed `p = 0.5` is therefore more aligned with the post-matching design than estimators using a variable estimated propensity score.

# 4. `reloop_rf_po` also performs well, especially compared with other ReLOOP variants
## `reloop_rf_po` is generally more stable than `reloop_rf_v12` and `reloop_rf_interp`.
## This suggests that, in this simulation, the individuals-as-units RF imputation strategy works better than directly modeling pair-level contrasts.
## So among the P-LOOP/ReLOOP estimators, `reloop_rf_po` seems to be the strongest version.

# 5. IPW-style estimators are unstable in this post-matching setting
## `cf_ipw_rf` has large negative bias and high RMSE in almost every setting. `cf_aipw_rf` improves on pure IPW but is still often negatively biased and less stable than `cf_loop_p05` or `reloop_rf_po`.
## This suggests that using a full-data RF propensity score as a post-matching weighting score may be mismatched with the matched-pair target sample.
## After matching, the sample already has a forced 1:1 treated-control structure, while the RF propensity still reflects the original observational treatment assignment mechanism.

# Bottom-line conclusion
## The main pattern is: Nonlinear settings and bad matching are harder. Fixed-`p = 0.5` LOOP-style adjustment and `reloop_rf_po` are the most stable estimators.
## IPW-style estimators using RF propensity scores are unstable in this post-matching paired setting, likely because the propensity weighting does not align well with the matched-sample design.
