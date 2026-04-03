#!/usr/bin/env Rscript
# =============================================================================
# Elastic Net Logistic Regression for SGA10 Prediction
# =============================================================================
#
# Description:
#   Identifies neonatal metabolomic biomarkers predictive of small-for-
#   gestational-age status (10th percentile, SGA10) using elastic net
#   logistic regression (alpha = 0.5) with 10-fold cross-validated lambda
#   tuning optimized on AUROC.
#
# Pipeline:
#   1. Data cleaning and mean imputation
#   2. Stratified 75/25 train/test split
#   3. Elastic net (alpha = 0.5) with 10-fold CV, AUROC-optimized
#   4. Evaluation: ROC-AUC with 95% CI, AUPRC
#   5. Coefficient export and top-10 visualization
#
# Usage:
#   Rscript elastic_net_sga10_analysis.R
#
# Dependencies:
#   glmnet, caret, pROC, ggplot2, PRROC
#
# =============================================================================

# --- Configuration -----------------------------------------------------------

CONFIG <- list(
  seed_split     = 123,
  seed_model     = 233,
  train_fraction = 0.75,
  alpha          = 0.5,
  cv_folds       = 10,
  target_col     = "sga_who",
  coef_csv       = "elastic_net_sga10_coefficients.csv",
  output_plot    = "elastic_net_sga10_top_metabolites.png",
  plot_width     = 6,
  plot_height    = 4,
  plot_dpi       = 500
)

METABOLITE_COLS <- c(
  "ALA", "ARG", "C02", "C03DC", "C04", "C05", "C051", "C05DC", "C05OH",
  "C06", "C08", "C081", "C10", "C101", "C12", "C121", "C14", "C141",
  "C142", "C14OH", "C16", "C16OH", "C18", "C181", "C181OH", "C182",
  "C18OH", "CIT", "FC", "GLY", "MET", "ORN", "PHE", "R03D10", "R14_12",
  "R16O16", "R3_2", "R5_3", "R8_10", "RA_O", "RC_A", "RF_C", "RF_Y",
  "RL_A", "RO_C", "RV_F", "TYR", "VAL", "XLE", "PRO", "OXP", "SA"
)

# --- Dependencies ------------------------------------------------------------

required_packages <- c("glmnet", "caret", "pROC", "ggplot2", "PRROC")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
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

  # 2. Stratified train/test split ---------------------------------------------
  set.seed(CONFIG$seed_split)
  train_idx <- createDataPartition(
    df[[CONFIG$target_col]], p = CONFIG$train_fraction, list = FALSE
  )

  x_train <- as.matrix(df[train_idx, METABOLITE_COLS])
  x_test  <- as.matrix(df[-train_idx, METABOLITE_COLS])
  y_train <- factor(df[[CONFIG$target_col]][train_idx],
                    levels = c(0, 1), labels = c("No", "Yes"))
  y_test_numeric <- df[[CONFIG$target_col]][-train_idx]

  message(sprintf("Train: %d | Test: %d", nrow(x_train), nrow(x_test)))

  # 3. Elastic net with 10-fold CV, AUROC-optimized ----------------------------
  message("\n=== ELASTIC NET TRAINING ===")

  cv_control <- trainControl(
    method          = "cv",
    number          = CONFIG$cv_folds,
    classProbs      = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    index           = createFolds(y_train, CONFIG$cv_folds),
    verboseIter     = TRUE
  )

  tune_grid <- expand.grid(
    alpha  = CONFIG$alpha,
    lambda = 10^seq(-4, 1, length.out = 100)
  )

  set.seed(CONFIG$seed_model)
  enet_model <- train(
    x       = x_train,
    y       = y_train,
    method  = "glmnet",
    metric  = "ROC",
    tuneGrid = tune_grid,
    trControl = cv_control,
    family  = "binomial"
  )

  message(sprintf("Best lambda: %.6f | Best alpha: %.2f",
                  enet_model$bestTune$lambda, enet_model$bestTune$alpha))

  # 4. Evaluation on test set --------------------------------------------------
  message("\n=== EVALUATION ===")

  test_probs <- predict(enet_model, x_test, type = "prob")[, "Yes"]

  # ROC-AUC with 95% CI
  roc_obj <- roc(y_test_numeric, test_probs, quiet = TRUE)
  auc_val <- auc(roc_obj)
  auc_ci  <- ci.auc(roc_obj, conf.level = 0.95)
  message(sprintf("AUC: %.4f [95%% CI: %.4f - %.4f]",
                  auc_val, auc_ci[1], auc_ci[3]))

  plot.roc(roc_obj, main = "ROC Curve (Elastic Net)")

  # AUPRC
  pos_scores <- test_probs[y_test_numeric == 1]
  neg_scores <- test_probs[y_test_numeric == 0]
  pr_obj <- pr.curve(
    scores.class0 = pos_scores,
    scores.class1 = neg_scores,
    curve = FALSE
  )
  auprc_val <- pr_obj$auc.integral
  message(sprintf("AUPRC: %.4f", auprc_val))

  # 5. Coefficient export ------------------------------------------------------
  message("\n=== COEFFICIENT EXPORT ===")

  all_coefs <- as.matrix(
    coef(enet_model$finalModel, s = enet_model$bestTune$lambda)
  )[, 1]

  coef_export <- data.frame(
    variable    = names(all_coefs),
    coefficient = all_coefs,
    stringsAsFactors = FALSE
  )
  write.csv(coef_export, CONFIG$coef_csv, row.names = FALSE)
  message(sprintf("Coefficients saved: %s", CONFIG$coef_csv))

  # 6. Top-10 coefficient visualization ----------------------------------------
  message("\n=== COEFFICIENT PLOT ===")

  # Remove intercept
  feature_coefs <- all_coefs[-1]
  coef_df <- data.frame(
    metabolite      = names(feature_coefs),
    coefficient     = feature_coefs,
    abs_coefficient = abs(feature_coefs),
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
      title = sprintf("Top %d Metabolites (Elastic Net)", n_top),
      x = "Metabolite",
      y = "Coefficient Value"
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
  message(sprintf("Alpha             : %.2f", CONFIG$alpha))
  message(sprintf("Best lambda       : %.6f", enet_model$bestTune$lambda))
  message(sprintf("AUC               : %.4f [%.4f, %.4f]",
                  auc_val, auc_ci[1], auc_ci[3]))
  message(sprintf("AUPRC             : %.4f", auprc_val))

  n_nonzero <- sum(feature_coefs != 0)
  message(sprintf("Non-zero features : %d / %d", n_nonzero, length(feature_coefs)))

  invisible(list(
    model    = enet_model,
    auc      = as.numeric(auc_val),
    auc_ci   = as.numeric(auc_ci),
    auprc    = auprc_val,
    coef_df  = coef_df
  ))
}

# --- Entry Point --------------------------------------------------------------
results <- main()
