library(dRCT)
library(MatchIt)
library(randomForest)

simulate_observational_data <- function(
    N_units = 1000,
    scenario = c("aligned", "unrelated", "opposite"),
    alpha_x3 = 0.5,
    beta_x1 = 0.5,
    beta_x3 = 1.0,
    treat_intercept = -2.0,
    gamma = 0.5,
    gamma_x3 = 1.0,
    tau = 0,
    sigma_x = 1,
    sigma_y = 1,
    seed = NULL
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  
  X3 <- rnorm(N_units)
  X1 <- alpha_x3 * X3 + rnorm(N_units, sd = sigma_x)
  X2 <- rnorm(N_units)
  
  p_treat <- plogis(treat_intercept + beta_x3 * X3 + beta_x1 * X1)
  Z <- rbinom(N_units, 1, p_treat)
  
  if (scenario == "aligned") {
    Y0 <- gamma_x3 * X3 + gamma * X1 + rnorm(N_units, sd = sigma_y)
  }
  
  if (scenario == "unrelated") {
    Y0 <- gamma_x3 * X3 + gamma * X2 + rnorm(N_units, sd = sigma_y)
  }
  
  if (scenario == "opposite") {
    Y0 <- gamma_x3 * X3 - gamma * X1 + rnorm(N_units, sd = sigma_y)
  }
  
  Y1 <- Y0 + tau
  Y <- ifelse(Z == 1, Y1, Y0)
  
  data.frame(
    id = 1:N_units,
    scenario = scenario,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    p_treat = p_treat,
    Z = Z,
    Y0 = Y0,
    Y1 = Y1,
    Y = Y
  )
}

match_pairs_on_x3 <- function(dat) {
  m.out <- matchit(
    Z ~ X3,
    data = dat,
    method = "nearest",
    distance = "glm",
    ratio = 1,
    replace = FALSE
  )
  
  dat_matched <- as.data.frame(match.data(m.out))
  dat_matched$pair <- as.numeric(as.factor(dat_matched$subclass))
  
  dat_long <- do.call(
    rbind,
    lapply(
      split(dat_matched, dat_matched$pair),
      function(x) {
        x <- as.data.frame(x)
        x <- x[sample(seq_len(nrow(x))), ]
        x$unit <- seq_len(nrow(x))
        x
      }
    )
  )
  
  dat_long <- as.data.frame(dat_long)
  rownames(dat_long) <- NULL
  dat_long <- dat_long[order(dat_long$pair, dat_long$unit), ]
  
  matched_ids <- dat_long$id
  dat_remnant <- dat[!(dat$id %in% matched_ids), ]
  
  list(
    long = dat_long,
    remnant = dat_remnant,
    matchit = m.out
  )
}

add_remnant_y0_predictions <- function(
    dat_long,
    dat_remnant,
    ntree = 500,
    mtry = NULL,
    nodesize = 5
) {
  remnant_controls <- dat_remnant[dat_remnant$Z == 0, ]
  
  if (nrow(remnant_controls) < 10) {
    stop("Too few remnant controls to train the Y0 prediction model.")
  }
  
  if (is.null(mtry)) {
  mtry <- max(1, floor(p / 3))
  }
  
  y0_model <- randomForest(
    Y ~ X1 + X2 + X3,
    data = remnant_controls,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize
  )
  
  dat_long$Y0_hat_remnant <- as.numeric(
    predict(y0_model, newdata = dat_long)
  )
  
  list(
    long = dat_long,
    y0_model = y0_model
  )
}

make_pair_data <- function(dat_long) {
  pair_sizes <- table(dat_long$pair)
  
  if (any(pair_sizes != 2)) {
    stop("Each matched pair must contain exactly two units.")
  }
  
  treated_counts <- tapply(dat_long$Z, dat_long$pair, sum)
  
  if (any(treated_counts != 1)) {
    stop("Each matched pair must contain exactly one treated and one control unit.")
  }
  
  treated <- dat_long[dat_long$Z == 1, ]
  control <- dat_long[dat_long$Z == 0, ]
  
  treated <- treated[order(treated$pair), ]
  control <- control[order(control$pair), ]
  
  if (!identical(treated$pair, control$pair)) {
    stop("Treated and control units are not aligned by pair.")
  }
  
  dat_pair <- data.frame(
    pair = treated$pair,
    scenario = treated$scenario,
    id_treat = treated$id,
    id_ctrl = control$id,
    X1_treat = treated$X1,
    X1_ctrl = control$X1,
    X2_treat = treated$X2,
    X2_ctrl = control$X2,
    X3_treat = treated$X3,
    X3_ctrl = control$X3,
    Y0_hat_treat = treated$Y0_hat_remnant,
    Y0_hat_ctrl = control$Y0_hat_remnant,
    Y0_treat_unit = treated$Y0,
    Y0_ctrl_unit = control$Y0,
    Y1_treat_unit = treated$Y1,
    Y1_ctrl_unit = control$Y1,
    Y_treat = treated$Y,
    Y_ctrl = control$Y
  )
  
  dat_pair$DeltaX1 <- dat_pair$X1_treat - dat_pair$X1_ctrl
  dat_pair$DeltaX2 <- dat_pair$X2_treat - dat_pair$X2_ctrl
  dat_pair$DeltaX3 <- dat_pair$X3_treat - dat_pair$X3_ctrl
  
  dat_pair$DeltaY_obs <- dat_pair$Y_treat - dat_pair$Y_ctrl
  dat_pair$DeltaY0 <- dat_pair$Y0_treat_unit - dat_pair$Y0_ctrl_unit
  dat_pair$DeltaY0_hat <- dat_pair$Y0_hat_treat - dat_pair$Y0_hat_ctrl
  dat_pair$d_true <- dat_pair$DeltaY0
  
  dat_pair
}

simulate_matched_observational_data <- function(
    N_units = 1000,
    scenario = c("aligned", "unrelated", "opposite"),
    alpha_x3 = 0.5,
    beta_x1 = 0.5,
    beta_x3 = 1.0,
    treat_intercept = -2.0,
    gamma = 0.5,
    gamma_x3 = 1.0,
    tau = 0,
    sigma_x = 1,
    sigma_y = 1,
    rf_ntree = 500,
    rf_mtry = NULL,
    rf_nodesize = 5,
    seed = NULL
) {
  scenario <- match.arg(scenario)
  
  dat_full <- simulate_observational_data(
    N_units = N_units,
    scenario = scenario,
    alpha_x3 = alpha_x3,
    beta_x1 = beta_x1,
    beta_x3 = beta_x3,
    treat_intercept = treat_intercept,
    gamma = gamma,
    gamma_x3 = gamma_x3,
    tau = tau,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    seed = seed
  )
  
  matched <- match_pairs_on_x3(dat_full)
  
  pred <- add_remnant_y0_predictions(
    dat_long = matched$long,
    dat_remnant = matched$remnant,
    ntree = rf_ntree,
    mtry = rf_mtry,
    nodesize = rf_nodesize
  )
  
  dat_long <- pred$long
  dat_pair <- make_pair_data(dat_long)
  
  list(
    full = dat_full,
    long = dat_long,
    pair = dat_pair,
    remnant = matched$remnant,
    matchit = matched$matchit,
    y0_model = pred$y0_model
  )
}

get_tauhat <- function(x) {
  as.numeric(x["tauhat"])
}

estimate_all <- function(sim, weighted_imp = FALSE) {
  pair <- sim$pair
  long <- sim$long
  
  tau_simple <- mean(pair$DeltaY_obs)
  tau_rebar_x1 <- mean(pair$DeltaY_obs - pair$DeltaX1)
  tau_rebar_y0hat <- mean(pair$DeltaY_obs - pair$DeltaY0_hat)
  
  Z_covariates <- cbind(
    Y0_hat_remnant = long$Y0_hat_remnant,
    X1 = long$X1
  )
  
  reloop_rf_po <- p_loop(
    Y = long$Y,
    Tr = long$Z,
    Z = Z_covariates,
    P = long$pair,
    pred = p_rf_po,
    weighted_imp = weighted_imp
  )
  
  reloop_rf_v12 <- p_loop(
    Y = long$Y,
    Tr = long$Z,
    Z = Z_covariates,
    P = long$pair,
    pred = p_rf_v12,
    weighted_imp = weighted_imp
  )
  
  reloop_rf_interp <- p_loop(
    Y = long$Y,
    Tr = long$Z,
    Z = Z_covariates,
    P = long$pair,
    pred = p_rf_interp,
    weighted_imp = weighted_imp
  )
  
  data.frame(
    scenario = unique(pair$scenario),
    simple_difference = tau_simple,
    rebar_x1 = tau_rebar_x1,
    rebar_y0hat = tau_rebar_y0hat,
    reloop_rf_po = get_tauhat(reloop_rf_po),
    reloop_rf_v12 = get_tauhat(reloop_rf_v12),
    reloop_rf_interp = get_tauhat(reloop_rf_interp),
    true_tau = mean(
      c(
        pair$Y1_treat_unit - pair$Y0_treat_unit,
        pair$Y1_ctrl_unit - pair$Y0_ctrl_unit
      )
    ),
    mean_DeltaY0 = mean(pair$DeltaY0),
    mean_DeltaY0_hat = mean(pair$DeltaY0_hat),
    mean_abs_DeltaX1 = mean(abs(pair$DeltaX1)),
    mean_abs_DeltaX2 = mean(abs(pair$DeltaX2)),
    mean_abs_DeltaX3 = mean(abs(pair$DeltaX3)),
    mean_abs_DeltaY0_hat = mean(abs(pair$DeltaY0_hat)),
    cor_DeltaY0hat_DeltaY0 = cor(pair$DeltaY0_hat, pair$DeltaY0),
    cor_DeltaX1_DeltaY0 = cor(pair$DeltaX1, pair$DeltaY0),
    cor_DeltaX2_DeltaY0 = cor(pair$DeltaX2, pair$DeltaY0),
    cor_DeltaX3_DeltaY0 = cor(pair$DeltaX3, pair$DeltaY0),
    prop_unit1_treated = mean(long$Z[long$unit == 1]),
    n_treated_full = sum(sim$full$Z == 1),
    n_control_full = sum(sim$full$Z == 0),
    treatment_rate_full = mean(sim$full$Z),
    n_matched = nrow(long),
    n_pairs = nrow(pair),
    n_remnant = nrow(sim$remnant),
    n_remnant_controls = sum(sim$remnant$Z == 0)
  )
}

one_run <- function(
    i,
    N_units = 1000,
    treat_intercept = -2.0,
    rf_ntree = 500,
    weighted_imp = FALSE
) {
  dat_aligned <- simulate_matched_observational_data(
    N_units = N_units,
    scenario = "aligned",
    treat_intercept = treat_intercept,
    rf_ntree = rf_ntree,
    seed = 10000 + i
  )
  
  dat_unrelated <- simulate_matched_observational_data(
    N_units = N_units,
    scenario = "unrelated",
    treat_intercept = treat_intercept,
    rf_ntree = rf_ntree,
    seed = 20000 + i
  )
  
  dat_opposite <- simulate_matched_observational_data(
    N_units = N_units,
    scenario = "opposite",
    treat_intercept = treat_intercept,
    rf_ntree = rf_ntree,
    seed = 30000 + i
  )
  
  results <- rbind(
    estimate_all(dat_aligned, weighted_imp = weighted_imp),
    estimate_all(dat_unrelated, weighted_imp = weighted_imp),
    estimate_all(dat_opposite, weighted_imp = weighted_imp)
  )
  
  results$sim <- i
  results
}

args <- commandArgs(trailingOnly = TRUE)

task_id <- if (length(args) >= 1) as.integer(args[1]) else 1
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1

if (is.na(task_id) || is.na(n_tasks)) {
  stop("task_id and n_tasks must be integers.")
}

set.seed(123)

n_sims <- 500
N_units <- 1000
treat_intercept <- -2.0
rf_ntree <- 500
weighted_imp <- FALSE

if (task_id < 1 || task_id > n_tasks) {
  stop("task_id must be between 1 and n_tasks.")
}

all_sims <- seq_len(n_sims)
chunk_id <- cut(all_sims, breaks = n_tasks, labels = FALSE)
sim_chunks <- split(all_sims, chunk_id)
sims_this_task <- sim_chunks[[as.character(task_id)]]

if (length(sims_this_task) == 0) {
  stop("No simulations assigned to this task. Use fewer array tasks.")
}

cat("Task:", task_id, "of", n_tasks, "\n")
cat("Simulations:", paste(sims_this_task, collapse = ", "), "\n")

sim_results <- lapply(
  sims_this_task,
  one_run,
  N_units = N_units,
  treat_intercept = treat_intercept,
  rf_ntree = rf_ntree,
  weighted_imp = weighted_imp
)

sim_results <- do.call(rbind, sim_results)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

out_file <- sprintf("results/sim_results_task_%03d.rds", task_id)
saveRDS(sim_results, file = out_file)

cat("Saved:", out_file, "\n")
print(sim_results)