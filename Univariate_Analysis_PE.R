#!/usr/bin/env Rscript
# =============================================================================
# Univariate Metabolomic Analysis for Preeclampsia
# =============================================================================
#
# Description:
#   Performs univariate comparison of neonatal metabolite levels between
#   preeclampsia-exposed and unexposed newborns. Preeclampsia is defined as
#   a composite of superimposed, severe (spree), or mild (mpree) subtypes.
#
# Pipeline:
#   1. Construct composite preeclampsia indicator
#   2. Data cleaning and mean imputation
#   3. Wilcoxon rank-sum test per metabolite with Benjamini-Hochberg correction
#   4. Log2 fold-change estimation
#   5. Volcano plot visualization
#
# Usage:
#   Rscript univariate_preeclampsia_analysis.R
#
# Dependencies:
#   dplyr, EnhancedVolcano
#
# =============================================================================

# --- Configuration -----------------------------------------------------------

CONFIG <- list(
  target_col   = "pree",
  p_cutoff     = 0.1,
  fc_cutoff    = 0.1,
  output_plot  = "volcano_preeclampsia_metabolites.png",
  plot_width   = 8,
  plot_height  = 6,
  plot_dpi     = 600
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

required_packages <- c("dplyr", "EnhancedVolcano")
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

# --- Statistical Utilities ----------------------------------------------------

#' Run Wilcoxon rank-sum tests across all metabolites
#'
#' @param df     Cleaned data frame with metabolite columns and binary target.
#' @param target Name of the binary outcome column.
#' @param metabolites Character vector of metabolite column names.
#' @return Data frame with columns: metabolite, p_value, adjusted_p_value.
run_wilcoxon_tests <- function(df, target, metabolites) {
  p_values <- vapply(metabolites, function(met) {
    if (is.numeric(df[[met]]) && length(unique(df[[target]])) == 2) {
      wilcox.test(df[[met]] ~ df[[target]])$p.value
    } else {
      NA_real_
    }
  }, numeric(1))

  valid <- !is.na(p_values)
  adjusted <- p.adjust(p_values[valid], method = "BH")

  data.frame(
    metabolite       = names(p_values[valid]),
    p_value          = p_values[valid],
    adjusted_p_value = adjusted,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

#' Compute log2 fold-change for each metabolite
#'
#' Calculates log2(mean_case + 1) - log2(mean_control + 1) for each metabolite.
#' The +1 pseudocount prevents log(0).
#'
#' @param df     Cleaned data frame.
#' @param target Name of the binary outcome column (1 = case, 0 = control).
#' @param metabolites Character vector of metabolite column names.
#' @return Named numeric vector of log2 fold-changes.
compute_log2fc <- function(df, target, metabolites) {
  vapply(metabolites, function(met) {
    if (is.numeric(df[[met]]) && length(unique(df[[target]])) == 2) {
      mean_case    <- mean(df[[met]][df[[target]] == 1], na.rm = TRUE)
      mean_control <- mean(df[[met]][df[[target]] == 0], na.rm = TRUE)
      log2(mean_case + 1) - log2(mean_control + 1)
    } else {
      NA_real_
    }
  }, numeric(1))
}

# --- Main Analysis ------------------------------------------------------------

main <- function() {

  # 1. Data preparation --------------------------------------------------------
  message("=== DATA PREPARATION ===")

  # NOTE: Replace this line with your data loading step, e.g.:
  #   raw_data <- read.csv("path/to/metabolomics_data.csv")
  raw_data <- metabolite_cleans_for_use_with_R

  # Construct composite preeclampsia indicator
  raw_data <- raw_data %>%
    mutate(pree = ifelse(superimposed == 1 | spree == 1 | mpree == 1, 1, 0))

  df <- subset(raw_data, select = c(METABOLITE_COLS, CONFIG$target_col))
  df <- clean_metabolomics(df, CONFIG$target_col)

  n_case    <- sum(df[[CONFIG$target_col]] == 1)
  n_control <- sum(df[[CONFIG$target_col]] == 0)
  message(sprintf("Clean dataset: %d samples (%d pree, %d control), %d metabolites",
                  nrow(df), n_case, n_control, ncol(df) - 1))

  # 2. Wilcoxon rank-sum tests with BH correction -----------------------------
  message("\n=== WILCOXON RANK-SUM TESTS ===")
  results <- run_wilcoxon_tests(df, CONFIG$target_col, METABOLITE_COLS)

  # 3. Log2 fold-change -------------------------------------------------------
  log2fc <- compute_log2fc(df, CONFIG$target_col, results$metabolite)
  results$log2_fold_change <- log2fc

  n_sig <- sum(results$adjusted_p_value < CONFIG$p_cutoff)
  message(sprintf("Significant metabolites (BH-adjusted p < %.2f): %d / %d",
                  CONFIG$p_cutoff, n_sig, nrow(results)))
  print(results[order(results$adjusted_p_value), ], row.names = FALSE)

  # 4. Volcano plot ------------------------------------------------------------
  message("\n=== VOLCANO PLOT ===")

  volcano <- EnhancedVolcano(
    results,
    lab        = results$metabolite,
    x          = "log2_fold_change",
    y          = "adjusted_p_value",
    xlab       = bquote(~Log[2] ~ "fold change"),
    ylab       = bquote(~-Log[10] ~ italic(P)),
    pCutoff    = CONFIG$p_cutoff,
    FCcutoff   = CONFIG$fc_cutoff,
    col        = c("gray", "forestgreen", "royalblue", "red"),
    pointSize  = 3.5,
    colAlpha   = 0.8,
    drawConnectors  = TRUE,
    widthConnectors = 0.75,
    colConnectors   = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    xlim    = c(-0.5, 0.5),
    ylim    = c(0, 300),
    hline   = c(-log10(CONFIG$p_cutoff)),
    hlineCol = "black",
    vlineCol = "black",
    caption  = ""
  )

  print(volcano)
  ggsave(CONFIG$output_plot, plot = volcano,
         width = CONFIG$plot_width, height = CONFIG$plot_height,
         dpi = CONFIG$plot_dpi)
  message(sprintf("Plot saved: %s", CONFIG$output_plot))

  # 5. Summary -----------------------------------------------------------------
  message("\n=== SUMMARY ===")
  message(sprintf("Total metabolites tested : %d", nrow(results)))
  message(sprintf("Significant (BH < %.2f)  : %d", CONFIG$p_cutoff, n_sig))

  if (n_sig > 0) {
    sig_mets <- results$metabolite[results$adjusted_p_value < CONFIG$p_cutoff]
    message(sprintf("Significant metabolites  : %s", paste(sig_mets, collapse = ", ")))
  }

  invisible(results)
}

# --- Entry Point --------------------------------------------------------------
results <- main()
