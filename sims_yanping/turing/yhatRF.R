library(dRCT)
library(MatchIt)
library(randomForest)
library(dplyr)
library(tidyr)

# -----------------------------
# Helpers
# -----------------------------

clip01 <- function(x, eps = 0.02) {
  pmin(pmax(x, eps), 1 - eps)
}

get_tauhat <- function(x) {
  as.numeric(x["tauhat"])
}

safe_tauhat <- function(expr) {
  out <- tryCatch(expr, error = function(e) NA_real_)
  if (inherits(out, "loopEst")) return(get_tauhat(out))
  if (is.numeric(out) && length(out) == 1) return(out)
  NA_real_
}

make_pair_folds <- function(pair_id, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pairs <- unique(pair_id)
  pairs <- sample(pairs)
  fold_id_by_pair <- rep(seq_len(K), length.out = length(pairs))
  names(fold_id_by_pair) <- pairs
  as.integer(fold_id_by_pair[as.character(pair_id)])
}

rf_reg_predict <- function(train, test, outcome = "Y",
                           covars = c("X1", "X2", "X3"),
                           ntree = 300,
                           nodesize = 5) {
  if (nrow(train) < 10) {
    return(rep(mean(train[[outcome]], na.rm = TRUE), nrow(test)))
  }

  fit <- randomForest(
    x = train[, covars, drop = FALSE],
    y = train[[outcome]],
    ntree = ntree,
    nodesize = nodesize
  )

  as.numeric(predict(fit, newdata = test[, covars, drop = FALSE]))
}

rf_prop_predict <- function(train, test,
                            covars = c("X1", "X2", "X3"),
                            ntree = 300,
                            nodesize = 5) {
  if (length(unique(train$Z)) < 2 || nrow(train) < 20) {
    return(rep(mean(train$Z), nrow(test)))
  }

  z_factor <- factor(train$Z, levels = c(0, 1))

  fit <- randomForest(
    x = train[, covars, drop = FALSE],
    y = z_factor,
    ntree = ntree,
    nodesize = nodesize
  )

  pred <- predict(fit, newdata = test[, covars, drop = FALSE], type = "prob")
  clip01(as.numeric(pred[, "1"]))
}

# -----------------------------
# Data generating process
# -----------------------------

simulate_observational_data <- function(
    N_units = 1000,
    scenario = c("aligned", "unrelated", "opposite"),
    functional_form = c("linear", "nonlinear"),
    treat_intercept = -1.5,
    tau = 0,
    sigma_x = 1,
    sigma_y = 1,
    seed = NULL
) {
  scenario <- match.arg(scenario)
  functional_form <- match.arg(functional_form)

  if (!is.null(seed)) set.seed(seed)

  X3 <- rnorm(N_units)
  X1 <- 0.5 * X3 + rnorm(N_units, sd = sigma_x)
  X2 <- rnorm(N_units, sd = sigma_x)

  if (functional_form == "linear") {
    linpred <- treat_intercept + 1.0 * X3 + 0.5 * X1

    if (scenario == "aligned") {
      Y0 <- 1.0 * X3 + 0.5 * X1 + rnorm(N_units, sd = sigma_y)
    }

    if (scenario == "unrelated") {
      Y0 <- 1.0 * X3 + 0.5 * X2 + rnorm(N_units, sd = sigma_y)
    }

    if (scenario == "opposite") {
      Y0 <- 1.0 * X3 - 0.5 * X1 + rnorm(N_units, sd = sigma_y)
    }
  }

  if (functional_form == "nonlinear") {
    linpred <- treat_intercept +
      0.9 * X3 +
      0.5 * X1 +
      0.7 * X1 * X3 +
      0.5 * sin(X1)

    if (scenario == "aligned") {
      Y0 <- 1.0 * X3 +
        0.5 * X1 +
        0.7 * X1 * X3 +
        0.5 * sin(X1) +
        rnorm(N_units, sd = sigma_y)
    }

    if (scenario == "unrelated") {
      Y0 <- 1.0 * X3 +
        0.5 * X2 +
        0.7 * X2 * X3 +
        0.5 * sin(X2) +
        rnorm(N_units, sd = sigma_y)
    }

    if (scenario == "opposite") {
      Y0 <- 1.0 * X3 -
        0.5 * X1 -
        0.7 * X1 * X3 -
        0.5 * sin(X1) +
        rnorm(N_units, sd = sigma_y)
    }
  }

  p_treat <- clip01(plogis(linpred))
  Z <- rbinom(N_units, 1, p_treat)

  Y1 <- Y0 + tau
  Y <- ifelse(Z == 1, Y1, Y0)

  data.frame(
    id = seq_len(N_units),
    scenario = scenario,
    functional_form = functional_form,
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

# -----------------------------
# Matching
# -----------------------------

match_pairs <- function(dat,
                        functional_form = c("linear", "nonlinear"),
                        matching_quality = c("good", "bad")) {
  functional_form <- match.arg(functional_form)
  matching_quality <- match.arg(matching_quality)

  if (matching_quality == "bad") {
    fml <- Z ~ X3
  }

  if (matching_quality == "good" && functional_form == "linear") {
    fml <- Z ~ X1 + X3
  }

  if (matching_quality == "good" && functional_form == "nonlinear") {
    fml <- Z ~ X1 + X3 + I(X1 * X3) + I(sin(X1))
  }

  m.out <- matchit(
    fml,
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

# -----------------------------
# External Y0 prediction from remnant controls
# -----------------------------

add_remnant_y0_predictions <- function(
    dat_long,
    dat_remnant,
    covars = c("X1", "X2", "X3"),
    ntree = 300,
    mtry = NULL,
    nodesize = 5
) {
  remnant_controls <- dat_remnant[dat_remnant$Z == 0, ]

  if (nrow(remnant_controls) < 10) {
    stop("Too few remnant controls to train the Y0 prediction model.")
  }

  if (is.null(mtry)) {
    mtry <- max(1, floor(length(covars) / 3))
  }

  y0_model <- randomForest(
    x = remnant_controls[, covars, drop = FALSE],
    y = remnant_controls$Y,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize
  )

  dat_long$Y0_hat_remnant <- as.numeric(
    predict(y0_model, newdata = dat_long[, covars, drop = FALSE])
  )

  list(
    long = dat_long,
    y0_model = y0_model
  )
}

# -----------------------------
# Pair-level data
# -----------------------------

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
    functional_form = treated$functional_form,
    matching_quality = treated$matching_quality,
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

# -----------------------------
# Cross-fitted LOOP/IPW/AIPW estimators
# -----------------------------

crossfit_estimators <- function(
    long,
    remnant = NULL,
    K = 5,
    covars = c("X1", "X2", "X3"),
    ntree = 300,
    nodesize = 5,
    seed = NULL
) {
  long <- as.data.frame(long)
  long$fold <- make_pair_folds(long$pair, K = K, seed = seed)

  n <- nrow(long)
  mu1_hat <- rep(NA_real_, n)
  mu0_hat <- rep(NA_real_, n)
  e_hat <- rep(NA_real_, n)

  if (!is.null(remnant)) {
    remnant <- as.data.frame(remnant)
  }

  for (k in seq_len(K)) {
    test_idx <- which(long$fold == k)
    train_matched <- long[long$fold != k, ]

    if (!is.null(remnant)) {
      train_all <- bind_rows(train_matched, remnant)
    } else {
      train_all <- train_matched
    }

    test <- long[test_idx, ]

    train_t <- train_all[train_all$Z == 1, ]
    train_c <- train_all[train_all$Z == 0, ]

    mu1_hat[test_idx] <- rf_reg_predict(
      train = train_t,
      test = test,
      outcome = "Y",
      covars = covars,
      ntree = ntree,
      nodesize = nodesize
    )

    mu0_hat[test_idx] <- rf_reg_predict(
      train = train_c,
      test = test,
      outcome = "Y",
      covars = covars,
      ntree = ntree,
      nodesize = nodesize
    )

    e_hat[test_idx] <- rf_prop_predict(
      train = train_all,
      test = test,
      covars = covars,
      ntree = ntree,
      nodesize = nodesize
    )
  }

  e_hat <- clip01(e_hat)

  Y <- long$Y
  Z <- long$Z

  cf_loop_p05 <- mean(
    mu1_hat - mu0_hat +
      Z * (Y - mu1_hat) / 0.5 -
      (1 - Z) * (Y - mu0_hat) / 0.5
  )

  cf_ipw_rf <- mean(
    Z * Y / e_hat -
      (1 - Z) * Y / (1 - e_hat)
  )

  cf_aipw_rf <- mean(
    mu1_hat - mu0_hat +
      Z * (Y - mu1_hat) / e_hat -
      (1 - Z) * (Y - mu0_hat) / (1 - e_hat)
  )

  list(
    cf_loop_p05 = cf_loop_p05,
    cf_ipw_rf = cf_ipw_rf,
    cf_aipw_rf = cf_aipw_rf,
    mean_e_hat = mean(e_hat),
    sd_e_hat = sd(e_hat)
  )
}

# -----------------------------
# One simulated matched dataset
# -----------------------------

simulate_matched_observational_data <- function(
    N_units = 1000,
    scenario = c("aligned", "unrelated", "opposite"),
    functional_form = c("linear", "nonlinear"),
    matching_quality = c("good", "bad"),
    treat_intercept = -1.5,
    tau = 0,
    sigma_x = 1,
    sigma_y = 1,
    rf_ntree = 300,
    rf_mtry = NULL,
    rf_nodesize = 5,
    seed = NULL
) {
  scenario <- match.arg(scenario)
  functional_form <- match.arg(functional_form)
  matching_quality <- match.arg(matching_quality)

  dat_full <- simulate_observational_data(
    N_units = N_units,
    scenario = scenario,
    functional_form = functional_form,
    treat_intercept = treat_intercept,
    tau = tau,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    seed = seed
  )

  matched <- match_pairs(
    dat = dat_full,
    functional_form = functional_form,
    matching_quality = matching_quality
  )

  dat_long <- matched$long
  dat_long$matching_quality <- matching_quality

  pred <- add_remnant_y0_predictions(
    dat_long = dat_long,
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

# -----------------------------
# Estimate all estimators
# -----------------------------

estimate_all <- function(
    sim,
    weighted_imp = FALSE,
    rf_ntree_cf = 300,
    K_cf = 5,
    seed = NULL
) {
  pair <- sim$pair
  long <- sim$long
  remnant <- sim$remnant

  tau_simple <- mean(pair$DeltaY_obs)
  # tau_rebar_x1 <- mean(pair$DeltaY_obs - pair$DeltaX1)
  tau_rebar_y0hat <- mean(pair$DeltaY_obs - pair$DeltaY0_hat)

  Z_covariates <- cbind(
    Y0_hat_remnant = long$Y0_hat_remnant) #,
    # X1 = long$X1
  # )

  reloop_rf_po <- safe_tauhat(
    p_loop(
      Y = long$Y,
      Tr = long$Z,
      Z = Z_covariates,
      P = long$pair,
      pred = p_rf_po,
      weighted_imp = weighted_imp
    )
  )

  reloop_rf_v12 <- safe_tauhat(
    p_loop(
      Y = long$Y,
      Tr = long$Z,
      Z = Z_covariates,
      P = long$pair,
      pred = p_rf_v12,
      weighted_imp = weighted_imp
    )
  )

  reloop_rf_interp <- safe_tauhat(
    p_loop(
      Y = long$Y,
      Tr = long$Z,
      Z = Z_covariates,
      P = long$pair,
      pred = p_rf_interp,
      weighted_imp = weighted_imp
    )
  )

  cf <- crossfit_estimators(
    long = long,
    remnant = remnant,
    K = K_cf,
    covars = c("X1", "X2", "X3"),
    ntree = rf_ntree_cf,
    seed = seed
  )

  data.frame(
    scenario = unique(pair$scenario),
    functional_form = unique(pair$functional_form),
    matching_quality = unique(pair$matching_quality),

    simple_difference = tau_simple,
    # rebar_x1 = tau_rebar_x1,
    rebar_y0hat = tau_rebar_y0hat,
    reloop_rf_po = reloop_rf_po,
    reloop_rf_v12 = reloop_rf_v12,
    reloop_rf_interp = reloop_rf_interp,
    cf_loop_p05 = cf$cf_loop_p05,
    cf_ipw_rf = cf$cf_ipw_rf,
    cf_aipw_rf = cf$cf_aipw_rf,

    true_tau = mean(
      c(
        pair$Y1_treat_unit - pair$Y0_treat_unit,
        pair$Y1_ctrl_unit - pair$Y0_ctrl_unit
      )
    ),

    mean_DeltaY0 = mean(pair$DeltaY0),
    mean_DeltaY0_hat = mean(pair$DeltaY0_hat),
    mean_DeltaX1 = mean(pair$DeltaX1),
    mean_DeltaX2 = mean(pair$DeltaX2),
    mean_DeltaX3 = mean(pair$DeltaX3),
    mean_abs_DeltaX1 = mean(abs(pair$DeltaX1)),
    mean_abs_DeltaX2 = mean(abs(pair$DeltaX2)),
    mean_abs_DeltaX3 = mean(abs(pair$DeltaX3)),
    mean_abs_DeltaY0_hat = mean(abs(pair$DeltaY0_hat)),
    cor_DeltaY0hat_DeltaY0 = cor(pair$DeltaY0_hat, pair$DeltaY0),
    cor_DeltaX1_DeltaY0 = cor(pair$DeltaX1, pair$DeltaY0),
    cor_DeltaX2_DeltaY0 = cor(pair$DeltaX2, pair$DeltaY0),
    cor_DeltaX3_DeltaY0 = cor(pair$DeltaX3, pair$DeltaY0),
    prop_unit1_treated = mean(long$Z[long$unit == 1]),
    mean_e_hat_cf = cf$mean_e_hat,
    sd_e_hat_cf = cf$sd_e_hat,

    n_treated_full = sum(sim$full$Z == 1),
    n_control_full = sum(sim$full$Z == 0),
    treatment_rate_full = mean(sim$full$Z),
    n_matched = nrow(long),
    n_pairs = nrow(pair),
    n_remnant = nrow(remnant),
    n_remnant_controls = sum(remnant$Z == 0)
  )
}

# -----------------------------
# One Monte Carlo replication
# -----------------------------

one_run <- function(
    i,
    N_units = 1000,
    treat_intercept = -1.5,
    rf_ntree = 300,
    rf_ntree_cf = 300,
    weighted_imp = FALSE,
    K_cf = 5
) {
  grid <- expand.grid(
    scenario = c("aligned", "unrelated", "opposite"),
    functional_form = c("linear", "nonlinear"),
    matching_quality = c("good", "bad"),
    stringsAsFactors = FALSE
  )

  results <- lapply(seq_len(nrow(grid)), function(g) {
    scenario_g <- grid$scenario[g]
    functional_g <- grid$functional_form[g]
    matching_g <- grid$matching_quality[g]

    seed_g <- 100000 * i +
      ifelse(scenario_g == "aligned", 10000,
             ifelse(scenario_g == "unrelated", 20000, 30000)) +
      ifelse(functional_g == "linear", 1000, 2000) +
      ifelse(matching_g == "good", 100, 200)

    sim <- simulate_matched_observational_data(
      N_units = N_units,
      scenario = scenario_g,
      functional_form = functional_g,
      matching_quality = matching_g,
      treat_intercept = treat_intercept,
      rf_ntree = rf_ntree,
      seed = seed_g
    )

    est <- estimate_all(
      sim = sim,
      weighted_imp = weighted_imp,
      rf_ntree_cf = rf_ntree_cf,
      K_cf = K_cf,
      seed = seed_g + 1
    )

    est$sim <- i
    est
  })

  do.call(rbind, results)
}

# -----------------------------
# Turing array-task interface
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

task_id <- if (length(args) >= 1) as.integer(args[1]) else 1
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1

if (is.na(task_id) || is.na(n_tasks)) {
  stop("task_id and n_tasks must be integers.")
}

set.seed(123)

n_sims <- 500 #500
N_units <- 1000
treat_intercept <- -1.5
rf_ntree <- 300
rf_ntree_cf <- 300
weighted_imp <- FALSE
K_cf <- 5

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
  rf_ntree_cf = rf_ntree_cf,
  weighted_imp = weighted_imp,
  K_cf = K_cf
)

sim_results <- do.call(rbind, sim_results)

#dir.create("results", showWarnings = FALSE, recursive = TRUE)

out_file <- sprintf("results/sim_results_task_%03d.rds", task_id)
saveRDS(sim_results, file = out_file)

cat("Saved:", out_file, "\n")
print(sim_results)
