#!/usr/bin/env Rscript
# =============================================================================
# STABL Feature Selection with Ridge Regression for SGA10 Prediction
# =============================================================================
#
# Description:
#   Identifies neonatal metabolomic biomarkers predictive of small-for-
#   gestational-age status (10th percentile, SGA10) using a modified STABL
#   (Stability Approach to Regularization Based on subsampling) with
#   knockoff-augmented false discovery control, followed by Ridge regression
#   for classification.
#
#   Knockoffs are constructed via residual permutation: for each feature,
#   the residuals from regressing that feature on all others are permuted
#   and added back to the fitted values. This preserves the inter-feature
#   correlation structure while breaking the marginal association with the
#   response.
#
# Pipeline:
#   1. Data cleaning and mean imputation
#   2. STABL feature selection (bootstrap stability + knockoff FDP control)
#   3. Ridge regression on selected features (10-fold CV)
#   4. Evaluation: ROC-AUC with 95% CI, AUPRC
#
# Usage:
#   Rscript stabl_sga10_analysis.R
#
# Dependencies:
#   glmnet, MASS, pROC, ggplot2, PRROC (optional)
#
# =============================================================================

# --- Configuration -----------------------------------------------------------

CONFIG <- list(
  seed           = 223,
  train_fraction = 0.70,
  n_bootstrap    = 50,
  alpha_lasso    = 1,
  ridge_nfolds   = 10,
  lambda_ratio   = 30,
  target_col     = "sga_who",
  output_plot    = "stabl_ridge_top_metabolites.png",
  plot_width     = 6,
  plot_height    = 4,
  plot_dpi       = 600
)

METABOLITE_COLS <- c(
  "ALA", "ARG", "C02", "C03", "C03DC", "C04", "C04DC", "C05", "C051",
  "C05DC", "C05OH", "C06", "C08", "C081", "C10", "C101", "C12", "C121",
  "C14", "C141", "C142", "C14OH", "C16", "C16OH", "C18", "C181",
  "C181OH", "C182", "C18OH", "CIT", "FC", "GLY", "MET", "ORN", "OXP",
  "PHE", "PRO", "R03D10", "R05D03D", "R14_12", "R16O16", "R3_2",
  "R5_3", "R8_10", "RA_O", "RC_A", "RF_C", "RF_Y", "RL_A", "RO_C",
  "RV_F", "SA", "TYR", "VAL", "XLE"
)

# --- Dependencies ------------------------------------------------------------

required_packages <- c("glmnet", "MASS", "pROC", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# --- Knockoff Generation -----------------------------------------------------

#' Generate knockoff variables via residual permutation
#'
#' For each feature j, regresses j on all remaining features, then permutes
#' the residuals to construct a knockoff that preserves the inter-feature
#' correlation structure while breaking the marginal association with the
#' response. All operations are performed on standardized data, then
#' back-transformed to the original scale.
#'
#' @param x Numeric matrix (n x p) of predictor variables.
#' @return  Numeric matrix (n x p) of knockoff variables on the original scale.
create_knockoffs <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  x_std <- scale(x)

  x_knock <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    other_vars <- setdiff(seq_len(p), j)
    fit <- lm(x_std[, j] ~ x_std[, other_vars])
    perm_idx <- sample(n)
    x_knock[, j] <- fitted(fit) + residuals(fit)[perm_idx]
  }

  # Restore original location and scale
  x_knock <- sweep(x_knock, 2, attr(x_std, "scaled:scale"), "*")
  x_knock <- sweep(x_knock, 2, attr(x_std, "scaled:center"), "+")
  return(x_knock)
}

# --- STABL Feature Selection --------------------------------------------------

#' STABL: Stability selection with knockoff-based FDP control
#'
#' Performs bootstrap stability selection on the augmented design matrix
#' [X, X_knockoff]. The maximum selection frequency across the lasso
#' regularization path is computed for each variable. Knockoff frequencies
#' provide an empirical null distribution to estimate the false discovery
#' proportion (FDP), and the threshold minimizing FDP is used to select
#' features.
#'
#' @param x      Numeric matrix (n x p) of predictors.
#' @param y      Response vector (length n).
#' @param n_boot Number of bootstrap resamples.
#' @param alpha  Elastic net mixing parameter (1 = lasso, 0 = ridge).
#' @param family GLM family ("gaussian" or "binomial").
#' @return List with components:
#'   \item{selected}{Integer vector of selected feature indices (1-indexed).}
#'   \item{max_freq}{Numeric vector of max selection frequencies (length 2p).}
#'   \item{threshold}{Optimal frequency threshold.}
#'   \item{fdp}{Estimated false discovery proportion at threshold.}
run_stabl <- function(x, y, n_boot, alpha = 1, family = "gaussian") {
  n <- nrow(x)
  p <- ncol(x)

  # Construct knockoffs via residual permutation
  x_knock <- create_knockoffs(x)
  x_aug   <- cbind(x, x_knock)

  # Establish shared lambda sequence
  fit_init <- glmnet(x_aug, y, family = family)
  lambdas  <- fit_init$lambda
  lambdas  <- lambdas[max(lambdas) / lambdas < CONFIG$lambda_ratio]
  n_lambda <- length(lambdas)

  # Accumulate selection counts across bootstrap samples
  sel_counts <- matrix(0, nrow = 2 * p, ncol = n_lambda)

  message(sprintf("Running %d bootstrap samples...", n_boot))
  for (b in seq_len(n_boot)) {
    if (b %% 10 == 0) message(sprintf("  Bootstrap %d / %d", b, n_boot))
    boot_idx <- sample(n, replace = TRUE)
    fit_b <- glmnet(
      x_aug[boot_idx, ], y[boot_idx],
      alpha = alpha, family = family, lambda = lambdas
    )
    coef_b <- matrix(
      as.numeric(coef(fit_b, s = lambdas)[-1, ]),
      nrow = 2 * p
    )
    sel_counts <- sel_counts + (coef_b != 0)
  }

  # Maximum selection frequency across the regularization path
  freq     <- sel_counts / n_boot
  max_freq <- apply(freq, 1, max)

  # Knockoff-calibrated FDP estimation
  thresholds <- seq(0, 1, length.out = 100)
  fdp <- vapply(thresholds, function(t) {
    n_knock <- sum(max_freq[(p + 1):(2 * p)] > t) + 1
    n_total <- sum(max_freq > t)
    if (n_total == 0) n_total <- 1
    n_knock / n_total
  }, numeric(1))

  best_idx  <- which.min(fdp)
  threshold <- thresholds[best_idx]
  selected  <- which(max_freq[seq_len(p)] > threshold)

  message(sprintf("Threshold: %.4f | Estimated FDP: %.4f | Features selected: %d",
                  threshold, fdp[best_idx], length(selected)))

  list(
    selected  = selected,
    max_freq  = max_freq,
    threshold = threshold,
    fdp       = fdp[best_idx]
  )
}

# --- Data Cleaning Utilities --------------------------------------------------

#' Remove columns where all values are "NULL"
drop_null_columns <- function(df) {
  df[, !vapply(df, function(col) all(col == "NULL"), logical(1))]
}

#' Remove rows with missing or "NULL" values in the target column
drop_missing_target <- function(df, target) {
  df[!is.na(df[[target]]) & df[[target]] != "NULL", ]
}

#' Replace "NULL" strings with NA across all columns
null_to_na <- function(df) {
  data.frame(lapply(df, function(col) {
    col <- as.character(col)
    col[col == "NULL"] <- NA
    col
  }), stringsAsFactors = FALSE)
}

#' Coerce columns to numeric where possible
coerce_numeric <- function(df) {
  data.frame(lapply(df, function(col) {
    num <- suppressWarnings(as.numeric(col))
    if (all(is.na(num) == is.na(col))) num else col
  }), stringsAsFactors = FALSE)
}

#' Impute NA values with column means (numeric columns only)
impute_mean <- function(df) {
  data.frame(lapply(df, function(col) {
    if (is.numeric(col)) col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  }))
}

#' Full cleaning pipeline
clean_metabolomics <- function(df, target) {
  df |>
    drop_null_columns() |>
    drop_missing_target(target) |>
    null_to_na() |>
    coerce_numeric() |>
    impute_mean()
}

# --- Evaluation Utilities -----------------------------------------------------

#' Compute AUPRC via the trapezoidal rule
#'
#' @param y_true  Binary ground truth labels (0/1).
#' @param y_score Predicted probabilities.
#' @return Scalar AUPRC value.
compute_auprc <- function(y_true, y_score) {
  thresholds <- sort(unique(y_score), decreasing = TRUE)

  pr <- vapply(thresholds, function(t) {
    pred <- as.integer(y_score >= t)
    tp <- sum(pred == 1 & y_true == 1)
    fp <- sum(pred == 1 & y_true == 0)
    fn <- sum(pred == 0 & y_true == 1)
    c(
      precision = ifelse(tp + fp == 0, 0, tp / (tp + fp)),
      recall    = ifelse(tp + fn == 0, 0, tp / (tp + fn))
    )
  }, numeric(2))

  precision <- pr["precision", ]
  recall    <- pr["recall", ]
  ord <- order(recall)
  sum(diff(recall[ord]) * (precision[ord][-1] + precision[ord][-length(ord)]) / 2)
}

# --- Main Analysis ------------------------------------------------------------

main <- function() {

  # 1. Data preparation --------------------------------------------------------
  message("=== DATA PREPARATION ===")

  # NOTE: Replace this line with your data loading step, e.g.:
  #   raw_data <- read.csv("path/to/metabolomics_data.csv")
  raw_data <- metabolite_cleans_for_use_with_R

  df <- subset(raw_data, select = c(METABOLITE_COLS, CONFIG$target_col))
  df <- clean_metabolomics(df, CONFIG$target_col)

  message(sprintf("Clean dataset: %d samples, %d metabolites",
                  nrow(df), ncol(df) - 1))

  # 2. Train/test split --------------------------------------------------------
  set.seed(CONFIG$seed)
  train_idx <- sample(nrow(df), size = floor(CONFIG$train_fraction * nrow(df)))

  x_train <- as.matrix(df[train_idx,  !(names(df) %in% CONFIG$target_col)])
  y_train <- df[[CONFIG$target_col]][train_idx]
  x_test  <- as.matrix(df[-train_idx, !(names(df) %in% CONFIG$target_col)])
  y_test  <- df[[CONFIG$target_col]][-train_idx]

  # Safety: zero-fill any residual NAs in design matrices
  x_train[is.na(x_train)] <- 0
  x_test[is.na(x_test)]   <- 0

  message(sprintf("Train: %d | Test: %d", nrow(x_train), nrow(x_test)))

  # 3. STABL feature selection -------------------------------------------------
  message("\n=== STABL FEATURE SELECTION ===")
  stabl_out <- run_stabl(
    x_train, y_train,
    n_boot = CONFIG$n_bootstrap,
    alpha  = CONFIG$alpha_lasso,
    family = "binomial"
  )

  sel_idx   <- stabl_out$selected
  sel_names <- colnames(x_train)[sel_idx]
  message(sprintf("Selected biomarkers (%d): %s",
                  length(sel_names), paste(sel_names, collapse = ", ")))

  if (length(sel_idx) == 0) {
    stop("STABL selected zero features. Consider adjusting n_boot or lambda_ratio.")
  }

  # 4. Ridge regression on selected features -----------------------------------
  message("\n=== RIDGE REGRESSION ===")
  x_train_sel <- x_train[, sel_idx, drop = FALSE]
  x_test_sel  <- x_test[, sel_idx, drop = FALSE]

  ridge_cv <- cv.glmnet(
    x_train_sel, y_train,
    alpha = 0, family = "binomial", nfolds = CONFIG$ridge_nfolds
  )
  y_prob <- as.numeric(
    predict(ridge_cv, x_test_sel, s = "lambda.min", type = "response")
  )

  # 5. Evaluation --------------------------------------------------------------
  message("\n=== EVALUATION ===")

  # ROC-AUC with 95% CI
  roc_obj <- roc(y_test, y_prob, quiet = TRUE)
  auc_val <- auc(roc_obj)
  auc_ci  <- ci.auc(roc_obj, conf.level = 0.95)
  message(sprintf("AUC: %.4f [95%% CI: %.4f - %.4f]",
                  auc_val, auc_ci[1], auc_ci[3]))

  plot.roc(roc_obj, main = "ROC Curve (STABL + Ridge)")

  # AUPRC
  auprc_val <- compute_auprc(y_test, y_prob)
  message(sprintf("AUPRC (trapezoidal): %.4f", auprc_val))

  if (requireNamespace("PRROC", quietly = TRUE)) {
    pr_obj <- PRROC::pr.curve(
      scores.class0 = y_prob[y_test == 1],
      scores.class1 = y_prob[y_test == 0],
      curve = TRUE
    )
    message(sprintf("AUPRC (PRROC): %.4f", pr_obj$auc.integral))
    plot(pr_obj, main = "Precision-Recall Curve (STABL + Ridge)")
  }

  # 6. Coefficient visualization -----------------------------------------------
  message("\n=== COEFFICIENT PLOT ===")

  coefs <- as.numeric(coef(ridge_cv, s = "lambda.min"))[-1]
  coef_df <- data.frame(
    metabolite      = colnames(x_train_sel),
    coefficient     = coefs,
    abs_coefficient = abs(coefs),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[order(-coef_df$abs_coefficient), ]

  n_top <- min(10, nrow(coef_df))
  top_df <- coef_df[seq_len(n_top), ]
  top_df$direction <- ifelse(top_df$coefficient > 0, "positive", "negative")

  message("Top metabolites by |coefficient|:")
  print(top_df[, c("metabolite", "coefficient")], row.names = FALSE)

  p <- ggplot(top_df, aes(
    x = reorder(metabolite, coefficient),
    y = coefficient,
    fill = direction
  )) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = c(positive = "#0072B2", negative = "#E69F00")) +
    labs(
      title = sprintf("Top %d Metabolites (STABL + Ridge)", n_top),
      x = "Metabolite",
      y = "Ridge Coefficient"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))

  print(p)
  ggsave(CONFIG$output_plot, plot = p,
         width = CONFIG$plot_width, height = CONFIG$plot_height,
         dpi = CONFIG$plot_dpi)
  message(sprintf("Plot saved: %s", CONFIG$output_plot))

  # 7. Summary -----------------------------------------------------------------
  message("\n=== SUMMARY ===")
  message(sprintf("Features selected : %d", length(sel_names)))
  message(sprintf("AUC               : %.4f [%.4f, %.4f]",
                  auc_val, auc_ci[1], auc_ci[3]))
  message(sprintf("AUPRC             : %.4f", auprc_val))
  message(sprintf("Biomarkers        : %s", paste(sel_names, collapse = ", ")))

  invisible(list(
    selected_features = sel_names,
    auc       = as.numeric(auc_val),
    auc_ci    = as.numeric(auc_ci),
    auprc     = auprc_val,
    ridge_cv  = ridge_cv,
    stabl     = stabl_out,
    coef_df   = coef_df
  ))
}

# --- Entry Point --------------------------------------------------------------
results <- main()
