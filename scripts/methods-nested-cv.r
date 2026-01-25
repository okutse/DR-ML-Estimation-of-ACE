# Title: Doubly-robust machine learning estimation of the average causal effect of (selected) DEGs on breast cancer stage and survival outcomes
# Author: Amos Okutse
# Date: Jan 2025


# Comparative Methods:
#(1) DR Penalized Logistic Regression
#(2) DR Random Forest
#(3) DR BART
#(4) DR Logistic regression on categorized DEGs (low, mid, high expression levels based on quantiles)

## Load required packages
required_pkgs <- c("msigdbr", "GSVA", "ranger", "pROC", "ROCR", "ggridges", "viridis", "forcats", "tidyr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(sprintf("Install required packages before running this script: %s", paste(missing_pkgs, collapse = ", ")))
}

library(data.table)    # fast I/O and merging
library(dplyr)         # data manipulation
library(stringr)       # string helpers
library(survival)      # survival objects
library(glmnet)        # high-dimensional regression
library(grf)          # generalized random forests
library(BART)         # Bayesian Additive Regression Trees
library(splines)      # natural splines
library(ggplot2)      # plotting
library(ggrepel)      # labeling points cleanly
library(ggpubr)       # arranging plots
library(cowplot)      # arranging plots
library(parallel)     # parallel computing
library(doParallel)   # parallel backend for foreach
library(foreach)      # parallel loops
library(boot)        # for inverse logit function
library(ggvenn)      # for Venn diagrams)
library(tidyr)       # for data reshaping
library(forcats)     # for factor manipulation

library(pathview)         # for KEGG pathway visualization
library(msigdbr)       # MSigDB gene sets
library(GSVA)          # pathway activity scoring
library(ranger)        # fast random forests
library(pROC)          # AUC computation
library(ROCR)          # for PR-AUC computation

# install BiocManager if you don't have it (uncomment and run only run once)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# install the limma package using BiocManager
# BiocManager::install("limma")
library(limma)       # for DEG analysis
library(matrixStats) # for simple matrix operations
library(VennDiagram)  # for Venn diagrams
library(glmnet)       # high-dimensional regression
library(ggridges)     # for ridge plots
library(viridis)      # for color scales


# source helper files
source("scripts/helper.R")

# -----------------------------------------------------------------------------
# Helper utilities used across the workflow
# -----------------------------------------------------------------------------

#' Compute comprehensive classification metrics
compute_classification_metrics <- function(pred_prob, truth, threshold = 0.5) {
  pred_class <- ifelse(pred_prob >= threshold, 1, 0)
  
  # Confusion matrix components
  tp <- sum(pred_class == 1 & truth == 1)
  tn <- sum(pred_class == 0 & truth == 0)
  fp <- sum(pred_class == 1 & truth == 0)
  fn <- sum(pred_class == 0 & truth == 1)
  
  # Metrics
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  balanced_accuracy <- ((tp / max(tp + fn, 1)) + (tn / max(tn + fp, 1))) / 2
  precision <- tp / max(tp + fp, 1)
  recall <- tp / max(tp + fn, 1)  # sensitivity
  specificity <- tn / max(tn + fp, 1)
  f1 <- 2 * precision * recall / max(precision + recall, 1e-10)
  
  # ROC AUC
  roc_auc <- tryCatch({
    roc_obj <- pROC::roc(response = truth, predictor = pred_prob, quiet = TRUE, direction = "<")
    as.numeric(roc_obj$auc)
  }, error = function(e) NA_real_)
  
 # PR AUC 
  pr_auc <- tryCatch({
    fg <- pred_prob[truth == 1]
    bg <- pred_prob[truth == 0]
    if (length(fg) == 0 || length(bg) == 0) return(NA_real_)
    pr <- ROCR::prediction(pred_prob, truth)
    perf <- ROCR::performance(pr, "aucpr")
    as.numeric(perf@y.values[[1]])
  }, error = function(e) NA_real_)
  
  data.frame(
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    precision = precision,
    recall = recall,
    sensitivity = recall,
    specificity = specificity,
    f1 = f1,
    roc_auc = roc_auc,
    pr_auc = pr_auc,
    stringsAsFactors = FALSE
  )
}

#' Create calibration data for plotting
create_calibration_data <- function(pred_prob, truth, n_bins = 10) {
  bins <- cut(pred_prob, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  cal_df <- data.frame(pred = pred_prob, truth = truth, bin = bins) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean_pred = mean(pred, na.rm = TRUE),
      mean_obs = mean(truth, na.rm = TRUE),
      n = dplyr::n(),
      se = sqrt(mean_obs * (1 - mean_obs) / max(n, 1)),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(mean_pred))
  cal_df
}

#' Compute optimal threshold using Youden's J statistic
compute_youden_threshold <- function(pred_prob, truth) {
  tryCatch({
    roc_obj <- pROC::roc(response = truth, predictor = pred_prob, quiet = TRUE, direction = "<")
    coords <- pROC::coords(roc_obj, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
    if (is.data.frame(coords)) {
      # If multiple optimal thresholds, take the first one
      coords <- coords[1, ]
    }
    list(
      threshold = as.numeric(coords["threshold"]),
      sensitivity = as.numeric(coords["sensitivity"]),
      specificity = as.numeric(coords["specificity"]),
      youden_j = as.numeric(coords["sensitivity"]) + as.numeric(coords["specificity"]) - 1
    )
  }, error = function(e) {
    list(threshold = 0.5, sensitivity = NA_real_, specificity = NA_real_, youden_j = NA_real_)
  })
}

#' Compute optimal threshold from PR curve (F1-optimal)
compute_pr_optimal_threshold <- function(pred_prob, truth) {
  tryCatch({
    pred_obj <- ROCR::prediction(pred_prob, truth)
    perf_pr <- ROCR::performance(pred_obj, "prec", "rec")
    prec <- perf_pr@y.values[[1]]
    rec <- perf_pr@x.values[[1]]
    cutoffs <- pred_obj@cutoffs[[1]]
    
    # Compute F1 for each threshold
    f1 <- 2 * prec * rec / (prec + rec + 1e-10)
    f1[is.na(f1)] <- 0
    
    # Find the threshold with maximum F1
    best_idx <- which.max(f1)
    list(
      threshold = cutoffs[best_idx],
      precision = prec[best_idx],
      recall = rec[best_idx],
      f1 = f1[best_idx]
    )
  }, error = function(e) {
    list(threshold = 0.5, precision = NA_real_, recall = NA_real_, f1 = NA_real_)
  })
}

#' Compute classification metrics at a specific threshold
compute_metrics_at_threshold <- function(pred_prob, truth, threshold, threshold_name = "custom") {
  pred_class <- ifelse(pred_prob >= threshold, 1, 0)
  
  tp <- sum(pred_class == 1 & truth == 1)
  tn <- sum(pred_class == 0 & truth == 0)
  fp <- sum(pred_class == 1 & truth == 0)
  fn <- sum(pred_class == 0 & truth == 1)
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  balanced_accuracy <- ((tp / max(tp + fn, 1)) + (tn / max(tn + fp, 1))) / 2
  precision <- tp / max(tp + fp, 1)
  recall <- tp / max(tp + fn, 1)
  specificity <- tn / max(tn + fp, 1)
  f1 <- 2 * precision * recall / max(precision + recall, 1e-10)
  
  data.frame(
    threshold_type = threshold_name,
    threshold = threshold,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    precision = precision,
    recall = recall,
    sensitivity = recall,
    specificity = specificity,
    f1 = f1,
    tp = tp, tn = tn, fp = fp, fn = fn,
    stringsAsFactors = FALSE
  )
}

#' Create confusion matrix visualization
create_confusion_matrix_plot <- function(pred_prob, truth, threshold, model_name, split_name, threshold_type = "0.5") {
  pred_class <- ifelse(pred_prob >= threshold, 1, 0)
  
  cm <- table(Predicted = factor(pred_class, levels = c(0, 1), labels = c("Early", "Late")),
              Actual = factor(truth, levels = c(0, 1), labels = c("Early", "Late")))
  
  cm_df <- as.data.frame(cm)
  
  # Add percentages
  total <- sum(cm_df$Freq)
  cm_df$Percentage <- round(cm_df$Freq / total * 100, 1)
  cm_df$Label <- sprintf("%d\n(%.1f%%)", cm_df$Freq, cm_df$Percentage)
  
  ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 5, color = "white") +
    scale_fill_gradient(low = "#3C8DAD", high = "#1F4E5F", name = "Count") +
    labs(
      title = sprintf("Confusion Matrix: %s (%s)", model_name, split_name),
      subtitle = sprintf("Threshold: %.3f (%s)", threshold, threshold_type),
      x = "Actual Label",
      y = "Predicted Label"
    ) +
    theme_nature() +
    theme(legend.position = "none") +
    coord_fixed()
}

#' Create ROC curve with threshold points
create_roc_curve_with_thresholds <- function(roc_val, roc_test, youden_val, youden_test, pr_opt_val, pr_opt_test, model_name) {
  # Extract ROC data
  roc_val_df <- data.frame(
    fpr = 1 - roc_val$specificities,
    tpr = roc_val$sensitivities,
    split = "Validation"
  )
  roc_test_df <- data.frame(
    fpr = 1 - roc_test$specificities,
    tpr = roc_test$sensitivities,
    split = "Test"
  )
  roc_data <- dplyr::bind_rows(roc_val_df, roc_test_df)
  
  # Threshold points
  threshold_points <- data.frame(
    fpr = c(1 - youden_val$specificity, 1 - youden_test$specificity,
            1 - compute_spec_at_threshold(roc_val, pr_opt_val$threshold),
            1 - compute_spec_at_threshold(roc_test, pr_opt_test$threshold)),
    tpr = c(youden_val$sensitivity, youden_test$sensitivity,
            compute_sens_at_threshold(roc_val, pr_opt_val$threshold),
            compute_sens_at_threshold(roc_test, pr_opt_test$threshold)),
    threshold_type = c("Youden", "Youden", "PR-Opt", "PR-Opt"),
    split = c("Validation", "Test", "Validation", "Test"),
    threshold = c(youden_val$threshold, youden_test$threshold,
                  pr_opt_val$threshold, pr_opt_test$threshold)
  )
  
  ggplot(roc_data, aes(x = fpr, y = tpr, color = split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(data = threshold_points, aes(shape = threshold_type), size = 4, stroke = 1.5) +
    scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
    scale_shape_manual(values = c("Youden" = 17, "PR-Opt" = 15)) +
    labs(
      title = sprintf("ROC Curve: %s", model_name),
      subtitle = sprintf("Val AUC: %.3f | Test AUC: %.3f", as.numeric(roc_val$auc), as.numeric(roc_test$auc)),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Split",
      shape = "Threshold"
    ) +
    coord_equal() +
    theme_nature() +
    theme(legend.position = "right")
}

#' Helper to compute sensitivity at a threshold
compute_sens_at_threshold <- function(roc_obj, threshold) {
  coords <- pROC::coords(roc_obj, threshold, input = "threshold", ret = "sensitivity")
  if (is.na(coords) || length(coords) == 0) return(NA_real_)
  as.numeric(coords)
}

#' Helper to compute specificity at a threshold
compute_spec_at_threshold <- function(roc_obj, threshold) {
  coords <- pROC::coords(roc_obj, threshold, input = "threshold", ret = "specificity")
  if (is.na(coords) || length(coords) == 0) return(NA_real_)
  as.numeric(coords)
}

#' Create PR curve with threshold points
create_pr_curve_with_thresholds <- function(pred_val, truth_val, pred_test, truth_test, youden_val, youden_test, pr_opt_val, pr_opt_test, model_name) {
  # Compute PR curves
  pr_val <- ROCR::prediction(pred_val, truth_val)
  pr_test <- ROCR::prediction(pred_test, truth_test)
  
  perf_val <- ROCR::performance(pr_val, "prec", "rec")
  perf_test <- ROCR::performance(pr_test, "prec", "rec")
  
  pr_val_df <- data.frame(
    recall = perf_val@x.values[[1]],
    precision = perf_val@y.values[[1]],
    split = "Validation"
  ) %>% dplyr::filter(!is.na(precision))
  
  pr_test_df <- data.frame(
    recall = perf_test@x.values[[1]],
    precision = perf_test@y.values[[1]],
    split = "Test"
  ) %>% dplyr::filter(!is.na(precision))
  
  pr_data <- dplyr::bind_rows(pr_val_df, pr_test_df)
  
  # Threshold points
  get_prec_rec_at_threshold <- function(pred, truth, threshold) {
    pred_class <- ifelse(pred >= threshold, 1, 0)
    tp <- sum(pred_class == 1 & truth == 1)
    fp <- sum(pred_class == 1 & truth == 0)
    fn <- sum(pred_class == 0 & truth == 1)
    prec <- tp / max(tp + fp, 1)
    rec <- tp / max(tp + fn, 1)
    c(precision = prec, recall = rec)
  }
  
  youden_val_pr <- get_prec_rec_at_threshold(pred_val, truth_val, youden_val$threshold)
  youden_test_pr <- get_prec_rec_at_threshold(pred_test, truth_test, youden_test$threshold)
  pr_opt_val_pr <- get_prec_rec_at_threshold(pred_val, truth_val, pr_opt_val$threshold)
  pr_opt_test_pr <- get_prec_rec_at_threshold(pred_test, truth_test, pr_opt_test$threshold)
  
  threshold_points <- data.frame(
    recall = c(youden_val_pr["recall"], youden_test_pr["recall"],
               pr_opt_val_pr["recall"], pr_opt_test_pr["recall"]),
    precision = c(youden_val_pr["precision"], youden_test_pr["precision"],
                  pr_opt_val_pr["precision"], pr_opt_test_pr["precision"]),
    threshold_type = c("Youden", "Youden", "PR-Opt", "PR-Opt"),
    split = c("Validation", "Test", "Validation", "Test"),
    threshold = c(youden_val$threshold, youden_test$threshold,
                  pr_opt_val$threshold, pr_opt_test$threshold)
  )
  
  # Baseline precision (proportion of positives)
  baseline <- mean(c(truth_val, truth_test))
  
  ggplot(pr_data, aes(x = recall, y = precision, color = split)) +
    geom_hline(yintercept = baseline, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(data = threshold_points, aes(shape = threshold_type), size = 4, stroke = 1.5) +
    scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
    scale_shape_manual(values = c("Youden" = 17, "PR-Opt" = 15)) +
    labs(
      title = sprintf("Precision-Recall Curve: %s", model_name),
      x = "Recall",
      y = "Precision",
      color = "Split",
      shape = "Threshold"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_nature() +
    theme(legend.position = "right")
}

#' Compute permutation importance for a model
compute_permutation_importance <- function(model, X, y, predict_fn, n_permutations = 10, seed = 202501) {
  set.seed(seed)
  baseline_pred <- predict_fn(model, X)
  baseline_auc <- compute_auc(baseline_pred, y)
  
  importance <- sapply(colnames(X), function(var) {
    perm_aucs <- replicate(n_permutations, {
      X_perm <- X
      X_perm[, var] <- sample(X_perm[, var])
      perm_pred <- predict_fn(model, X_perm)
      compute_auc(perm_pred, y)
    })
    baseline_auc - mean(perm_aucs, na.rm = TRUE)
  })
  
  data.frame(
    variable = names(importance),
    importance = as.numeric(importance),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(importance))
}

theme_nature <- function(base_size = 11, base_family = "Helvetica") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E0E0E0", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = "#4A4A4A"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    )
}

stratified_split <- function(ids, strata, train_frac = 0.7, val_frac = 0.15, seed = 202501) {
  stopifnot(length(ids) == length(strata))
  set.seed(seed)
  strata_factor <- factor(ifelse(is.na(strata), "missing", as.character(strata)))
  train_ids <- character(0)
  val_ids <- character(0)
  test_ids <- character(0)
  for (lvl in levels(strata_factor)) {
    lvl_ids <- ids[strata_factor == lvl]
    n <- length(lvl_ids)
    if (n == 0) next
    if (n <= 2) {
      train_ids <- c(train_ids, lvl_ids[1])
      if (n >= 2) {
        test_ids <- c(test_ids, lvl_ids[2])
      }
      next
    }
    shuffled <- sample(lvl_ids, size = n, replace = FALSE)
    n_train <- max(1, round(train_frac * n))
    n_val <- max(1, round(val_frac * n))
    while (n_train + n_val >= n) {
      if (n_val > 1) {
        n_val <- n_val - 1
      } else if (n_train > 1) {
        n_train <- n_train - 1
      } else {
        break
      }
    }
    n_test <- n - n_train - n_val
    if (n_test <= 0) {
      n_val <- max(1, n_val - 1)
      n_test <- n - n_train - n_val
    }
    if (n_test <= 0) {
      n_train <- max(1, n_train - 1)
      n_test <- n - n_train - n_val
    }
    idx_train <- seq_len(n_train)
    idx_val <- seq.int(n_train + 1, n_train + n_val)
    train_ids <- c(train_ids, shuffled[idx_train])
    if (n_val > 0) {
      val_ids <- c(val_ids, shuffled[idx_val])
    }
    remaining_idx <- setdiff(seq_len(n), c(idx_train, idx_val))
    if (length(remaining_idx)) {
      test_ids <- c(test_ids, shuffled[remaining_idx])
    }
  }
  list(
    train = unique(train_ids),
    val = unique(val_ids),
    test = unique(test_ids)
  )
}

coarsen_factor <- function(x, min_n = 5, other_label = "Other") {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr) | x_chr == ""] <- "missing"
  counts <- table(x_chr)
  rare_levels <- names(counts[counts < min_n])
  if (length(rare_levels)) {
    x_chr[x_chr %in% rare_levels] <- other_label
  }
  factor(x_chr)
}

build_strata_key <- function(df, vars, k, seed = 202501) {
  set.seed(seed)
  vars <- vars[vars %in% names(df)]
  if (!length(vars)) {
    return(list(key = factor(rep("all", nrow(df))), vars_used = character(0)))
  }
  df_work <- df
  for (v in vars) {
    df_work[[v]] <- coarsen_factor(df_work[[v]], min_n = k)
  }
  vars_used <- vars
  repeat {
    key <- interaction(df_work[vars_used], drop = TRUE, lex.order = TRUE)
    if (all(table(key) >= k) || length(vars_used) <= 1) {
      break
    }
    vars_used <- vars_used[-length(vars_used)]
  }
  list(key = key, vars_used = vars_used)
}

create_stratified_folds_ids <- function(ids, strata_key, k = 5, seed = 202501) {
  set.seed(seed)
  folds <- vector("list", k)
  strata_levels <- levels(factor(strata_key))
  for (lvl in strata_levels) {
    idx <- which(strata_key == lvl)
    if (!length(idx)) next
    fold_assign <- sample(rep_len(seq_len(k), length(idx)))
    for (f in seq_len(k)) {
      folds[[f]] <- c(folds[[f]], ids[idx][fold_assign == f])
    }
  }
  lapply(folds, unique)
}

create_stratified_folds_ids2 <- function(ids, strata, k = 5, seed = 202501) {
  stopifnot(length(ids) == length(strata))
  set.seed(seed)
  strata <- as.character(strata)
  if (any(is.na(strata))) stop("Strata contains NA after alignment. Fix ID matching first.")
  
  folds <- vector("list", k)
  for (lvl in sort(unique(strata))) {
    idx <- which(strata == lvl)
    fold_assign <- sample(rep_len(seq_len(k), length(idx)))
    for (f in seq_len(k)) {
      folds[[f]] <- c(folds[[f]], ids[idx][fold_assign == f])
    }
  }
  lapply(folds, unique)
}


create_nested_splits <- function(ids, df, strata_vars, k_outer = 5, k_inner = 5, seed = 202501) {
  outer_info <- build_strata_key(df, strata_vars, k = k_outer, seed = seed)
  outer_folds <- create_stratified_folds_ids(ids, outer_info$key, k = k_outer, seed = seed)
  nested <- vector("list", k_outer)
  for (o in seq_len(k_outer)) {
    outer_test <- outer_folds[[o]]
    outer_train <- setdiff(ids, outer_test)
    outer_df <- df[match(outer_train, df$patient_id), , drop = FALSE]
    inner_info <- build_strata_key(outer_df, strata_vars, k = k_inner, seed = seed + o)
    inner_folds <- create_stratified_folds_ids(outer_train, inner_info$key, k = k_inner, seed = seed + o)
    inner_splits <- lapply(inner_folds, function(val_ids) {
      list(train = setdiff(outer_train, val_ids), val = val_ids)
    })
    nested[[o]] <- list(
      outer = list(train = outer_train, test = outer_test, strata_vars = outer_info$vars_used),
      inner = inner_splits
    )
  }
  nested
}

fit_imputer <- function(df, variables) {
  df_sub <- as.data.frame(df)[, variables, drop = FALSE]
  numeric_mask <- vapply(df_sub, function(col) is.numeric(col) || is.integer(col), logical(1))
  num_cols <- variables[numeric_mask]
  cat_cols <- variables[!numeric_mask]
  num_stats <- if (length(num_cols)) {
    stats::setNames(vapply(num_cols, function(col) {
      median(df_sub[[col]], na.rm = TRUE)
    }, numeric(1)), num_cols)
  } else numeric(0)
  cat_modes <- if (length(cat_cols)) {
    stats::setNames(vapply(cat_cols, function(col) {
      vals <- df_sub[[col]]
      vals <- vals[!is.na(vals)]
      if (!length(vals)) {
        "missing"
      } else {
        names(sort(table(vals), decreasing = TRUE))[1]
      }
    }, character(1)), cat_cols)
  } else character(0)
  cat_levels <- lapply(cat_cols, function(col) {
    lvls <- levels(factor(df_sub[[col]]))
    unique(c(lvls, cat_modes[[col]]))
  })
  names(cat_levels) <- cat_cols
  list(
    variables = variables,
    num_cols = num_cols,
    cat_cols = cat_cols,
    num_stats = num_stats,
    cat_modes = cat_modes,
    cat_levels = cat_levels
  )
}

apply_imputer <- function(df, imputer) {
  out <- as.data.frame(df)[, imputer$variables, drop = FALSE]
  for (col in imputer$num_cols) {
    vals <- as.numeric(out[[col]])
    stat <- imputer$num_stats[[col]]
    if (is.na(stat)) stat <- 0
    vals[is.na(vals)] <- stat
    out[[col]] <- vals
  }
  for (col in imputer$cat_cols) {
    vals <- as.character(out[[col]])
    mode_val <- imputer$cat_modes[[col]]
    if (is.na(mode_val) || mode_val == "") {
      mode_val <- "missing"
    }
    vals[is.na(vals)] <- mode_val
    unknown_mask <- !vals %in% imputer$cat_levels[[col]]
    if (any(unknown_mask)) {
      vals[unknown_mask] <- mode_val
    }
    out[[col]] <- factor(vals, levels = imputer$cat_levels[[col]])
  }
  out
}

compute_auc <- function(pred, truth) {
  if (length(unique(truth)) < 2) {
    return(NA_real_)
  }
  roc_obj <- pROC::roc(response = truth, predictor = pred, quiet = TRUE, direction = "<")
  as.numeric(roc_obj$auc)
}

binary_c_index <- function(pred, truth) {
  valid <- !is.na(pred) & !is.na(truth)
  truth <- truth[valid]
  pred <- pred[valid]
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  ranks <- rank(pred, ties.method = "average")
  (sum(ranks[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

glmnet_learner <- function(family = c("gaussian", "binomial"), alpha = 1) {
  family <- match.arg(family)
  function(X, y) {
    X_mat <- as.matrix(X)
    y_vec <- as.numeric(y)
    fit <- glmnet::cv.glmnet(X_mat, y_vec, family = family, alpha = alpha, nfolds = 5)
    list(
      model = fit,
      predict = function(newX) {
        pred <- stats::predict(fit, newx = as.matrix(newX), s = "lambda.min", type = if (family == "binomial") "response" else "response")
        as.numeric(pred)
      }
    )
  }
}

ranger_learner <- function(task = c("regression", "classification"), num.trees = 1500, mtry = NULL) {
  task <- match.arg(task)
  function(X, y) {
    X_df <- as.data.frame(X)
    data <- cbind(y = y, X_df)
    if (task == "classification") {
      data$y <- factor(data$y, levels = c(0, 1))
    } else {
      data$y <- as.numeric(data$y)
    }
    fit <- ranger::ranger(
      y ~ .,
      data = data,
      num.trees = num.trees,
      mtry = ifelse(is.null(mtry), max(1, floor(sqrt(ncol(X_df)))), mtry),
      probability = task == "classification",
      respect.unordered.factors = "order",
      min.node.size = 10,
      seed = 202501
    )
    list(
      model = fit,
      predict = function(newX) {
        pred <- predict(fit, data = as.data.frame(newX))$predictions
        if (task == "classification") {
          if (is.matrix(pred)) {
            as.numeric(pred[, "1", drop = TRUE])
          } else {
            as.numeric(pred)
          }
        } else {
          as.numeric(pred)
        }
      }
    )
  }
}

dml_partial_linear <- function(Y, D, X, outcome_learner, treatment_learner, K = 5, seed = 202501) {
  n <- length(Y)
  stopifnot(length(D) == n, nrow(X) == n)
  set.seed(seed)
  folds <- sample(rep(1:K, length.out = n))
  y_tilde <- numeric(n)
  d_tilde <- numeric(n)
  for (k in seq_len(K)) {
    idx_train <- which(folds != k)
    idx_test <- which(folds == k)
    fit_g <- outcome_learner(X[idx_train, , drop = FALSE], Y[idx_train])
    g_hat <- fit_g$predict(X[idx_test, , drop = FALSE])
    fit_m <- treatment_learner(X[idx_train, , drop = FALSE], D[idx_train])
    m_hat <- fit_m$predict(X[idx_test, , drop = FALSE])
    y_tilde[idx_test] <- Y[idx_test] - g_hat
    d_tilde[idx_test] <- D[idx_test] - m_hat
  }
  dr_fit <- stats::lm(y_tilde ~ d_tilde)
  theta_hat <- stats::coef(dr_fit)["d_tilde"]
  se_hat <- summary(dr_fit)$coefficients["d_tilde", "Std. Error"]
  list(
    theta = theta_hat,
    se = se_hat,
    ci = c(theta_hat - 1.96 * se_hat, theta_hat + 1.96 * se_hat),
    folds = folds
  )
}

bart_dml_partial <- function(Y, D, X, split_label, K = 5, seed = 202501,
                             ntree = 150, ndpost = 750, nskip = 250) {
  n <- length(Y)
  if (!n || nrow(X) == 0) {
    return(dplyr::tibble())
  }
  stopifnot(length(D) == n, nrow(X) == n)
  set.seed(seed)
  folds <- sample(rep(1:K, length.out = n))
  y_tilde <- numeric(n)
  d_tilde <- numeric(n)
  for (k in seq_len(K)) {
    idx_train <- folds != k
    idx_test <- folds == k
    outcome_fit <- BART::pbart(
      x.train = X[idx_train, , drop = FALSE],
      y.train = Y[idx_train],
      x.test = X[idx_test, , drop = FALSE],
      ntree = ntree,
      ndpost = ndpost,
      nskip = nskip,
      sparse = TRUE,
      usequants = TRUE
    )
    treatment_fit <- BART::wbart(
      x.train = X[idx_train, , drop = FALSE],
      y.train = D[idx_train],
      x.test = X[idx_test, , drop = FALSE],
      ntree = ntree,
      ndpost = ndpost,
      nskip = nskip,
      sparse = TRUE,
      usequants = TRUE
    )
    g_hat <- as.numeric(outcome_fit$prob.test.mean)
    m_hat <- as.numeric(treatment_fit$yhat.test.mean)
    y_tilde[idx_test] <- Y[idx_test] - g_hat
    d_tilde[idx_test] <- D[idx_test] - m_hat
  }
  dr_fit <- stats::lm(y_tilde ~ d_tilde)
  theta_hat <- stats::coef(dr_fit)["d_tilde"]
  se_hat <- summary(dr_fit)$coefficients["d_tilde", "Std. Error"]
  dplyr::tibble(
    split = split_label,
    model = "DR-BART",
    theta = theta_hat,
    se = se_hat,
    ci_lower = theta_hat - 1.96 * se_hat,
    ci_upper = theta_hat + 1.96 * se_hat
  )
}

## Load the TCGA BRCA data files from data directory in DataFest root directory
cp  <- data.table::fread("data/Data_clinical_patient.txt")[5:.N] # drop 1-4th rows
cs   <- data.table::fread("data/Data_clinical_sample.txt")[5:.N] # drop 1-4th rows
mp     <- fread("data/Data_methylation_promoters_rrbs.txt") # maybe augment with PCA exp (90%)
mexp     <- data.table::fread("data/Data_mrna_illumina_microarray.txt")
mt      <- data.table::fread("data/data_mutations.txt")                         


################################################################################
# Clinical patient and sample data ----
################################################################################

# check any missing values in cp columns
# Find all "" and replace with NA
cs <- cs %>%
  dplyr::mutate_all(~dplyr::na_if(., ""))

## Rename the cs 
cs <- cs %>% dplyr::mutate(`HER2 Status` = dplyr::na_if(`HER2 Status`, "")) %>% 
  dplyr::rename(patient_id = `#Patient Identifier`,
              sample = `Sample Identifier`,
              cancer_type = `Cancer Type`,
              cancer_type_detailed = `Cancer Type Detailed`,
              er_status = `ER Status`,
              her2_status = `HER2 Status`,
              neoplasm_histologic_grade = `Neoplasm Histologic Grade`,
              oncotree_code = `Oncotree Code`,
              pr_status = `PR Status`,
              sample_type = `Sample Type`,
              tumor_size = `Tumor Size`,
              tumor_stg = `Tumor Stage`,
              tmb_nonsynonymous = `TMB (nonsynonymous)`) %>% 
  # add binarized tumor stage from levels levels(as.factor(cs$`Tumor Stage`)) "0" "1" "2" "3" "4"
  dplyr::mutate(tumor_stage = as.factor(case_when(
    str_detect(tumor_stg, "0|1")  ~ "early",
    str_detect(tumor_stg, "2|3|4") ~ "late",
    TRUE ~ NA_character_
  ))) %>% 
  # add observability mask for tumor stage (1 = observed, 0 = missing) to address missing outcome data
  dplyr::mutate(tumor_stage_mask = if_else(!is.na(tumor_stage), 1, 0)) %>% 
  # add inverse weights for tumor stage to address class imbalance
  dplyr::mutate(tumor_stage_invwt = case_when(
    tumor_stage == "early" ~ 1 / sum(tumor_stage == "early", na.rm = TRUE),
    tumor_stage == "late" ~ 1 / sum(tumor_stage == "late", na.rm = TRUE),
    TRUE ~ NA_real_
  )) %>% 
  # add ER status mask (1 = observed, 0 = missing) to address missing confounder data
  dplyr::mutate(er_status_mask = if_else(!is.na(er_status), 1, 0)) %>% 
  # convert covariates to factors
  dplyr::mutate(
    #cancer_type = as.factor(cancer_type), # all breast cancer hence not meaningful
    cancer_type_detailed = as.factor(cancer_type_detailed),
    er_status = as.factor(er_status),
    her2_status = as.factor(her2_status), 
    her2_status_mask = if_else(!is.na(her2_status), 1, 0),
    neoplasm_histologic_grade = as.factor(neoplasm_histologic_grade),
    neoplasm_histologic_grade_mask = if_else(!is.na(neoplasm_histologic_grade), 1, 0),
    oncotree_code = as.factor(oncotree_code), # no missing
    pr_status = as.factor(pr_status),
    pr_status_mask = if_else(!is.na(pr_status), 1, 0),
    # sample_type = as.factor(sample_type) # all primary hence not meaningful
    tumor_size = as.numeric(tumor_size),
    tumor_size_mask = if_else(!is.na(tumor_size), 1, 0),
    tumor_stg = as.numeric(tumor_stg)
  ) %>% 
  dplyr::select(-c(cancer_type, sample_type)) # drop non-informative covariates

sapply(cs , function(x) sum(is.na(x)) )

# rename the cp data file
# Find all "" and replace with NA
cp <- cp %>%
  dplyr::mutate_all(~dplyr::na_if(., ""))

cp <- cp %>%
  #dplyr::mutate(`Sex` = dplyr::na_if(`Sex`, "")) %>% all female hence not meaningful
  dplyr::rename(patient_id = `#Patient Identifier`,
              lymph_nodes_examined_positive = `Lymph nodes examined positive`,
              nottingham_prognostic_index = `Nottingham prognostic index`,
              cellularity = `Cellularity`,
              chemotherapy = `Chemotherapy`,
              cohort = `Cohort`,
              er_status_measured_by_ihc = `ER status measured by IHC`,
              her2_status_measured_by_snp6 = `HER2 status measured by SNP6`,
              hormone_therapy = `Hormone Therapy`,
              inferred_menopausal_state = `Inferred Menopausal State`,
              #sex = `Sex`, # not meaningful
              integrative_cluster = `Integrative Cluster`,
              age_at_diagnosis = `Age at Diagnosis`,
              overall_survival_months = `Overall Survival (Months)`,
              overall_survival_status = `Overall Survival Status`,
              pam50_claudin_low_subtype = `Pam50 + Claudin-low subtype`,
              three_gene_classifier_subtype = `3-Gene classifier subtype`,
              patients_vital_status = `Patient's Vital Status`,
              primary_tumor_laterality = `Primary Tumor Laterality`,
              radio_therapy = `Radio Therapy`,
              tumor_other_histologic_subtype = `Tumor Other Histologic Subtype`,
              type_of_breast_surgery = `Type of Breast Surgery`,
              relapse_free_status_months = `Relapse Free Status (Months)`,
              relapse_free_status = `Relapse Free Status`) %>% 
dplyr::mutate(
              # categorical features
              relapse_free_status = as.factor(relapse_free_status),
              overall_survival_status = as.factor(overall_survival_status),
              cellularity = as.factor(cellularity),
              chemotherapy = as.factor(chemotherapy),
              er_status_measured_by_ihc = as.factor(er_status_measured_by_ihc),
              her2_status_measured_by_snp6 = as.factor(her2_status_measured_by_snp6),
              hormone_therapy = as.factor(hormone_therapy),
              inferred_menopausal_state = as.factor(inferred_menopausal_state),
              integrative_cluster = as.factor(integrative_cluster),
              pam50_claudin_low_subtype = as.factor(pam50_claudin_low_subtype),
              three_gene_classifier_subtype = as.factor(three_gene_classifier_subtype),
              patients_vital_status = as.factor(patients_vital_status),
              primary_tumor_laterality = as.factor(primary_tumor_laterality),
              radio_therapy = as.factor(radio_therapy),
              tumor_other_histologic_subtype = as.factor(tumor_other_histologic_subtype),
              type_of_breast_surgery = as.factor(type_of_breast_surgery),
              # numeric features
              age_at_diagnosis = as.numeric(age_at_diagnosis), 
              lymph_nodes_examined_positive = as.numeric(lymph_nodes_examined_positive),
              nottingham_prognostic_index = as.numeric(nottingham_prognostic_index),
              cohort = as.numeric(cohort),
              overall_survival_months = as.numeric(overall_survival_months), 
              relapse_free_status_months = as.numeric(relapse_free_status_months) 
              ) #%>% 


# loop to check and create masks for all columns in cp with missing data
for (col in colnames(cp)) {
  if (any(is.na(cp[[col]]))) {
    mask_col <- paste0(col, "_mask")
    cp[[mask_col]] <- if_else(!is.na(cp[[col]]), 1, 0)
  }
}

sapply(cp, function(x) sum(is.na(x)) )

# merge clinical sample and patient data
dt <- merge(cs, cp, by = "patient_id", all.x = TRUE)
names(dt)
 
#################################################################################
# mRNA expression data ----
#################################################################################

# mexp: first column is gene symbol, second column is entrez id, rest are sample IDs 
# colnames are [1] "Hugo_Symbol"    "Entrez_Gene_Id" "MB-0362"        "MB-0346" etc
# not z-scores, likely log2 transformed intensities 
genes <- mexp$Hugo_Symbol
exp_mat <- as.matrix(mexp[ , -c(1,2), with = FALSE])

# simply drop genes with > 20% missing values
na_by_gene <- rowMeans(is.na(exp_mat))
exp_mat <- exp_mat[na_by_gene <= 0.2, , drop=FALSE]

# simply impute remaining NAs gene-wise by median
if (anyNA(exp_mat)) {
  meds <- rowMedians(exp_mat, na.rm = TRUE)
  idx <- which(is.na(exp_mat), arr.ind = TRUE)
  exp_mat[idx] <- meds[idx[,1]]
}

rownames(exp_mat) <- genes
colnames(exp_mat) <- colnames(mexp)[-c(1,2)]

# common samples between clinical and expression data and filter to only those samples
common_samples <- intersect(colnames(exp_mat), dt$patient_id) # 1980 samples
exp_mat <- exp_mat[ , common_samples, drop = FALSE]


dt2 <- dt %>%
  dplyr::filter(patient_id %in% common_samples) %>%
  dplyr::arrange(match(patient_id, colnames(exp_mat)))
# Keep tumor_stage (required for stratification) and don't drop other clinical columns yet
dt2 <- dt2[!is.na(dt2$tumor_stage), ] # filter to non-missing tumor stage
cat("dt2 after filtering to non-missing tumor stage: ", nrow(dt2), " rows\n")
cat("dt2 columns:", paste(colnames(dt2), collapse=", "), "\n")
exp_mat <- exp_mat[, dt2$patient_id, drop = FALSE] # keep only subjects with observed tumor stage to align ids
grp <- factor(dt2$tumor_stage, levels = c("early", "late")) # early = 1338; late = 128; n = 1466 samples
table(grp)

# check whether exp_mat are already normalized log2 (range is 0 to 16)
summary(as.numeric(exp_mat))

# omit all NA in exp_mat if any after median imputation
exp_mat <- exp_mat[complete.cases(exp_mat), , drop = FALSE]

# Pre-processing steps for mRNA expression data
## Nested CV split before preprocessing using stratification by tumor stage + covariates
sample_ids <- colnames(exp_mat)
n_total <- length(sample_ids)
K_outer <- 10
K_inner <- 9
# With the current fold construction (inner folds drawn from remaining outer folds),
# K_inner cannot exceed K_outer - 1.
if (K_inner > (K_outer - 1)) {
  stop(sprintf(
    "K_inner (%d) must be <= K_outer - 1 (%d) with the current fold construction.",
    K_inner, K_outer - 1
  ))
}
# strata_vars <- c("tumor_stage", "cancer_type_detailed", "er_status")
# strata_vars <- c("tumor_stage", "er_status")
strata_vars <- c("tumor_stage")

strata_info <- build_strata_key(dt2, strata_vars, k = K_outer, seed = 202501)
# cv_folds <- create_stratified_folds_ids(sample_ids, strata_info$key, k = K_outer, seed = 202501)
# if (length(cv_folds) != K_outer) {
#   stop("Failed to create nested CV folds with the requested K_outer.")
# }

# build folds using *aligned* strata
strata <- dt2$tumor_stage[match(sample_ids, dt2$patient_id)]
cv_folds <- create_stratified_folds_ids2(sample_ids, strata, k = K_outer, seed = 202501)

sapply(cv_folds, length)

total_fold_configs <- K_outer * K_inner
cat(sprintf(
  "Nested CV config: K_outer=%d, K_inner=%d (total fold configs: %d)\n",
  K_outer, K_inner, total_fold_configs
))

collapse_duplicate_genes <- function(mat) {
  rn <- rownames(mat)
  if (anyDuplicated(rn)) {
    warning("Detected duplicated gene symbols; collapsing by mean expression.")
    summed <- rowsum(mat, rn)
    counts <- as.numeric(table(rn)[rownames(summed)])
    mat <- sweep(summed, 1, counts, "/")
  }
  mat
}

get_gsva_top_pathways <- function(train_ids, val_ids, test_ids, fold_label, gsva_selection_split = "val", n_top = 20) {
  # This function extracts the top pathways from GSVA ranking without running the full analysis
  cat(fold_label, " | Extracting top pathways | n_train=", length(train_ids), " n_val=", length(val_ids), " n_test=", length(test_ids), "\n")
  
  # Split the expression matrix
  train_mat_raw <- exp_mat[, train_ids, drop = FALSE]
  val_mat_raw <- exp_mat[, val_ids, drop = FALSE]
  test_mat_raw <- exp_mat[, test_ids, drop = FALSE]
  
  expr_train <- collapse_duplicate_genes(train_mat_raw)
  expr_val <- collapse_duplicate_genes(val_mat_raw)
  expr_test <- collapse_duplicate_genes(test_mat_raw)
  
  expr_splits <- list(
    train = expr_train,
    val = expr_val,
    test = expr_test,
    split_ids = list(train = train_ids, val = val_ids, test = test_ids)
  )
  
  # Align clinical data splits with expression data splits
  clinical_splits <- lapply(expr_splits$split_ids, function(ids) align_clinical_split(dt2, ids))
  
  # Curate breast cancer relevant gene sets from MSigDB
  bc_terms <- msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "C2",
    subcollection = "CP"
  ) %>%
    dplyr::filter(grepl("breast|cancer|carcinoma|oncogenic|tumor", gs_name, ignore.case = TRUE))
  
  if (!nrow(bc_terms)) {
    bc_terms <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }
  
  bc_gene_sets <- split(bc_terms$gene_symbol, bc_terms$gs_name)
  
  if (!gsva_selection_split %in% names(expr_splits) || !ncol(expr_splits[[gsva_selection_split]])) {
    stop(sprintf("Invalid gsva_selection_split '%s'", gsva_selection_split))
  }
  
  gsva_select_expr <- expr_splits[[gsva_selection_split]]
  gsva_select_clin <- clinical_splits[[gsva_selection_split]]
  
  bc_gene_sets <- lapply(bc_gene_sets, function(genes) {
    intersect(genes, rownames(gsva_select_expr))
  })
  bc_gene_sets <- bc_gene_sets[lengths(bc_gene_sets) >= 10]
  
  if (!length(bc_gene_sets)) {
    stop("No gene sets remaining after filtering for minimum size.")
  }
  
  select_ids <- colnames(gsva_select_expr)
  select_stage <- gsva_select_clin$tumor_stage[match(select_ids, gsva_select_clin$patient_id)]
  stage_factor <- factor(select_stage)
  
  if (any(is.na(stage_factor)) || length(unique(stage_factor)) < 2) {
    stop(sprintf(
      "Fold %s: stage_factor has NA or < 2 levels. Cannot proceed with GSVA ranking.",
      fold_label
    ))
  }
  
  bc_gene_sets <- lapply(bc_gene_sets, function(genes) {
    intersect(genes, rownames(gsva_select_expr))
  })
  bc_gene_sets <- bc_gene_sets[lengths(bc_gene_sets) >= 10]
  
  if (!length(bc_gene_sets)) {
    stop("No gene sets remaining after second filtering.")
  }
  
  gsva_rank_param <- GSVA::gsvaParam(
    exprData = gsva_select_expr,
    geneSets = bc_gene_sets,
    minSize = 10,
    maxSize = 500,
    kcdf = "Gaussian",
    verbose = FALSE
  )
  
  gsva_rank_scores <- GSVA::gsva(gsva_rank_param)
  
  mean_early <- matrixStats::rowMeans2(gsva_rank_scores[, stage_factor == "early", drop = FALSE])
  mean_late <- matrixStats::rowMeans2(gsva_rank_scores[, stage_factor == "late", drop = FALSE])
  delta <- mean_late - mean_early
  
  pvals <- apply(gsva_rank_scores, 1, function(x) {
    stats::t.test(x[stage_factor == "late"], x[stage_factor == "early"])$p.value
  })
  
  gsva_rank_tbl <- data.frame(
    ID = names(bc_gene_sets),
    mean_early = mean_early,
    mean_late = mean_late,
    delta = delta,
    p_value = pvals,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(abs(delta)))
  
  top_pathways <- gsva_rank_tbl %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::pull(ID)
  
  return(list(
    top_pathways = top_pathways,
    gsva_rank_tbl = gsva_rank_tbl,
    fold_label = fold_label
  ))
}

run_fold <- function(train_ids, val_ids, test_ids, fold_label, results_dir, gsva_selection_split = "val", selected_pathway_id = NULL) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  cat(fold_label, " | n_train=", length(train_ids), " n_val=", length(val_ids), " n_test=", length(test_ids), "\n")

  split_stage_summary <- lapply(
    list(train = train_ids, val = val_ids, test = test_ids),
    function(ids) {
      stages <- dt2$tumor_stage[match(ids, dt2$patient_id)]
      table(factor(stages, levels = c("early", "late")), useNA = "ifany")
    }
  )
  print(sprintf("Tumor stage distribution by split (%s):", fold_label))
  print(split_stage_summary)
  capture.output(split_stage_summary, file = file.path(results_dir, "split_stage_summary.txt"))

  # Split the expression matrix
  train_mat_raw <- exp_mat[, train_ids, drop = FALSE]
  val_mat_raw <- exp_mat[, val_ids, drop = FALSE]
  test_mat_raw <- exp_mat[, test_ids, drop = FALSE]
  
  ## Preprocess
  # qn_reference <- compute_qn_reference(train_mat_raw)
  # train_qn <- apply_qn(train_mat_raw, qn_reference)
  # val_qn <- apply_qn(val_mat_raw, qn_reference)
  # test_qn <- apply_qn(test_mat_raw, qn_reference)
  # 
  # train_log2 <- safe_log2(train_qn)
  # val_log2 <- safe_log2(val_qn)
  # test_log2 <- safe_log2(test_qn)
  # 
  # train_log2 <- collapse_duplicate_genes(train_log2)
  # val_log2 <- collapse_duplicate_genes(val_log2)
  # test_log2 <- collapse_duplicate_genes(test_log2)
  # 
  # gene_means <- rowMeans(train_log2, na.rm = TRUE)
  # gene_sds <- apply(train_log2, 1, sd, na.rm = TRUE)
  # gene_sds[is.na(gene_sds) | gene_sds == 0] <- 1  # guard against zero variance
  
  #expr_train <- zscore_with_stats(train_log2, gene_means, gene_sds)
  #expr_val <- zscore_with_stats(val_log2, gene_means, gene_sds)
  #expr_test <- zscore_with_stats(test_log2, gene_means, gene_sds)
  
  
  # expr_train <- train_log2
  # expr_val <- val_log2
  # expr_test <- test_log2
  
  expr_train = collapse_duplicate_genes(train_mat_raw)
  expr_val = collapse_duplicate_genes(val_mat_raw)
  expr_test = collapse_duplicate_genes(test_mat_raw)
  
  expr_splits <- list(
    train = expr_train,
    val = expr_val,
    test = expr_test,
    split_ids = list(train = train_ids, val = val_ids, test = test_ids),
    stats = list() # list(qn_reference = qn_reference, means = gene_means, sds = gene_sds)
  )
  
  # Align clinical data splits with expression data splits
  clinical_splits <- lapply(expr_splits$split_ids, function(ids) align_clinical_split(dt2, ids))
  stopifnot(
    identical(colnames(expr_splits$train), clinical_splits$train$patient_id),
    identical(colnames(expr_splits$val), clinical_splits$val$patient_id),
    identical(colnames(expr_splits$test), clinical_splits$test$patient_id)
  )
  
  ## Diagnostic plots on train set
  train_plot_df <- data.frame(
    value = c(as.vector(train_mat_raw), as.vector(expr_train)),
    state = rep(c("Raw", "Processed"),
                times = c(length(train_mat_raw), length(expr_train)))
  ) %>%
    dplyr::filter(!is.na(value))
  
  # p_train_raw <- train_plot_df %>%
  #   dplyr::filter(state == "Raw") %>%
  #   ggplot2::ggplot(ggplot2::aes(x = value)) +
  #   ggplot2::geom_density(fill = "#F8766D", alpha = 0.6) +
  #   ggplot2::labs(title = "Train distribution before preprocessing",
  #                 x = "Expression", y = "Density") +
  #   ggplot2::theme_minimal()
  # 
  # p_train_processed <- train_plot_df %>%
  #   dplyr::filter(state != "Raw") %>%
  #   ggplot2::ggplot(ggplot2::aes(x = value)) +
  #   ggplot2::geom_density(fill = "#00BFC4", alpha = 0.6) +
  #   ggplot2::labs(title = "Train distribution after preprocessing",
  #                 x = "Z-score", y = "Density") +
  #   ggplot2::theme_minimal()
  # 
  # train_preprocessing_check <- cowplot::plot_grid(p_train_raw, p_train_processed, ncol = 2)
  
  #################################################################################
  ## DEG identification using limma-voom and edgeR on training set only
  ###################################################################################
  ## keep genes with at least 20% samples having expression > 1 (z-score)
  keep_genes <- apply(expr_train, 1, function(x) sum(x > 1) >= 0.2 * length(x))
  dt_train_all <- clinical_splits$train
  stage_complete <- !is.na(dt_train_all$tumor_stage)
  dt_train <- dt_train_all[stage_complete, , drop = FALSE]
  
  expr_train_stage <- expr_train[, dt_train$patient_id, drop = FALSE]
  expr_filt <- expr_train_stage[keep_genes, , drop = FALSE]
  
  grp_train <- factor(dt_train$tumor_stage, levels = c("early", "late"))
  stopifnot(length(grp_train) == ncol(expr_filt))
  
  design <- stats::model.matrix(~ grp_train) # keep groups for IDs in train set only
  fit <- limma::lmFit(expr_filt, design)
  fit <- limma::eBayes(fit)
  fit$genes <- data.frame(gene = rownames(expr_filt))
  
  # limma results: Late vs Early 
  res_limma <- limma::topTable(
    fit,
    coef = "grp_trainlate",
    number = Inf,
    adjust.method = "BH"
  ) %>%
    as.data.frame()
  
  # argument on cutoff selection. note the small estimated effects overall. The middle 50% of logFC is roughly 
  # [−0.024, 0.027], and the most extreme values are only about −0.95 to +1.01
  summary(res_limma)
  
  # thresholds
  alpha_fdr <- 0.10 # conservative for exploratory analysis
  lfc_cut   <- 0.2          # ~1.15-fold
  top_n_sig <- 10
  top_n_nsig <- 10
  
  # pull rownames in as gene IDs
  if (!"gene" %in% colnames(res_limma)) {
    res_limma <- res_limma %>% tibble::rownames_to_column("gene")
  }
  
  # annotate direction + DEG class on the same dataset
  res_annot <- res_limma %>%
    mutate(
      neglog10_adjP = -log10(adj.P.Val),
      direction = case_when(
        logFC >=  lfc_cut ~ "Up",
        logFC <= -lfc_cut ~ "Down",
        TRUE              ~ "None"
      ),
      deg_class = case_when(
        adj.P.Val < alpha_fdr & logFC >=  lfc_cut ~ "Upregulated",
        adj.P.Val < alpha_fdr & logFC <= -lfc_cut ~ "Downregulated",
        TRUE                                  ~ "Not significant"
      )
    )
  
  # select DEGs for downstream (GSVA / pathway enrichment / causal) 
  deg_limma_selected <- res_annot %>%
    filter(adj.P.Val < alpha_fdr, abs(logFC) >= lfc_cut) %>%
    arrange(adj.P.Val)
  
  if (!nrow(deg_limma_selected)) {
    warning("No DEGs passed the specified FDR/logFC thresholds; using top 50 genes by adjusted p-value as a fallback.")
    deg_limma_selected <- res_annot %>%
      arrange(adj.P.Val) %>%
      slice_head(n = 50)
  }
  
  deg_genes <- deg_limma_selected$gene
  length(deg_genes)
  
  # Separate up vs down gene sets (useful for enrichment or signed analyses)
  deg_genes_up   <- deg_limma_selected %>% filter(logFC >=  lfc_cut) %>% pull(gene)
  deg_genes_down <- deg_limma_selected %>% filter(logFC <= -lfc_cut) %>% pull(gene)
  
  # Table of up and down regulated DEGs
  deg_tbl_up <- deg_limma_selected %>%
    filter(deg_class == "Upregulated") %>%
    arrange(adj.P.Val, desc(logFC)) %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
    slice_head(n = 10)
  
  deg_tbl_down <- deg_limma_selected %>%
    filter(deg_class == "Downregulated") %>%
    arrange(adj.P.Val, logFC) %>%   # most negative first
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
    slice_head(n = 10)
  
  deg_tbl <- bind_rows(deg_tbl_up, deg_tbl_down)
  
  knitr::kable(
    deg_tbl,
    format = "latex",
    booktabs = TRUE,
    digits = 4,
    caption = "Top Upregulated and Downregulated DEGs identified by limma"
  )
  
  # Volcano plot with correct coloring for up, down, and non-sig genes with labels
  volcano_data <- res_annot
  # Label selection:
  # For Up/Down: top N by abs(logFC)
  # For Not significant: top N by a stable "extremeness" score
  label_data <- bind_rows(
    volcano_data %>%
      filter(deg_class %in% c("Upregulated", "Downregulated")) %>%
      arrange(desc(abs(logFC))) %>%
      group_by(deg_class) %>%
      slice_head(n = top_n_sig) %>%
      ungroup(),
    volcano_data %>%
      filter(deg_class == "Not significant") %>%
      mutate(extremeness = neglog10_adjP * abs(logFC)) %>%
      arrange(desc(extremeness)) %>%
      slice_head(n = top_n_nsig)
  )
  volcano_data$deg_class <- factor(
    volcano_data$deg_class,
    levels = c("Upregulated", "Downregulated", "Not significant")
  )
  
  volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neglog10_adjP, color = as.factor(deg_class))) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(alpha_fdr), linetype = "dashed", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,   
      box.padding = 0.3,
      point.padding = 0.2,
      segment.alpha = 0.6
    ) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not significant" = "black"
    )) +
    labs(
      title = "(a)", #Volcano plot of DEGs in late vs early breast cancer tumor stages
      x = "log2 Fold Change (Late vs Early)",
      y = "-log10(FDR)",
      color = "DEG class"
    ) +
    theme_nature()
  
  volcano_plot
  # save the 600 dpi image in results folder
  ggsave(
    filename = file.path(results_dir, "[a]limma_deg_volcano_plot.png"),
    plot = volcano_plot,
    width = 6,
    height = 5,
    dpi = 600
  )
  
  # Top DEG genes for downstream analyses (fallback to overall top-ranked genes if up/down tables are empty)
  top_gene_candidates <- unique(c(deg_tbl_up$gene, deg_tbl_down$gene))
  if (!length(top_gene_candidates)) {
    top_gene_candidates <- deg_limma_selected %>%
      arrange(adj.P.Val) %>%
      pull(gene)
  }
  top_genes <- top_gene_candidates[top_gene_candidates %in% rownames(expr_splits$train)]
  
  if (!length(top_genes)) {
    stop("No DEG genes overlap with expression matrix after filtering; cannot proceed to downstream analyses.")
  }
  
  ## Expression of top DEGs
  exp_top <- expr_splits$train[top_genes, dt_train$patient_id, drop = FALSE]
  exp_top_df <- as.data.frame(t(exp_top)) %>%
    dplyr::mutate(patient_id = rownames(.)) %>%
    dplyr::left_join(dt_train %>% dplyr::select(patient_id, tumor_stage), by = "patient_id") %>%
    tidyr::pivot_longer(cols = all_of(top_genes), names_to = "gene", values_to = "expression")
  
  ## Boxplots of top DEGs expression by tumor stage 
  boxplot_top_degs <- ggplot(exp_top_df, aes(x = tumor_stage, y = expression, fill = tumor_stage)) +
    geom_boxplot() +
    facet_wrap(~ gene, scales = "free_y") +
    labs(title = "(b)", # Expression of top DEGs by breast cancer tumor stage [not strictly necessary except to show distributions]
         x = "Tumor Stage",
         y = "Expression scores") +
    theme(legend.position = "none")+
    theme_nature()
  boxplot_top_degs
  # save the 600 dpi image in results folder
  ggsave(boxplot_top_degs,
         filename = file.path(results_dir, "[b]top_degs_boxplot.png"),
         width = 8,
         height = 6,
         dpi = 600
  )
  
  #################################################################################
  # GSVA-driven pathway ranking ----
  ################################################################################
  # In this section, we use GSVA's built-in ranking to score pathway activity directly
  # (no GSEA-derived gene ranks). We then select the pathway with the strongest
  # GSVA score contrast between tumor stages in the validation split (nested CV).
  
  # Curate breast cancer relevant gene sets from MSigDB (C2:CGP). Fallback to Hallmark if empty.
  # retrieve gene sets and their member genes from the human genome collections using the molecular signatures database (MSigDB) curated gene set canonical pathways (C2:CGP) related to breast cancer.
  bc_terms <- msigdbr::msigdbr(
    species = "Homo sapiens", 
    collection = "C2",
    subcollection = "CP" # canonical pathways: 
  ) %>%
    dplyr::filter(stringr::str_detect(gs_name, "BREAST|MAMMARY|BRCA|ERBB2|HER2|ESTROGEN|LUMINAL|BASAL", negate = FALSE)) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct()
  
  if (!nrow(bc_terms)) {
    message("No breast cancer-specific CGP sets found; defaulting to Hallmark collection.")
    bc_terms <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::distinct()
  }
  
  
  bc_gene_sets <- split(bc_terms$gene_symbol, bc_terms$gs_name)
  
  if (!gsva_selection_split %in% names(expr_splits) || !ncol(expr_splits[[gsva_selection_split]])) {
    gsva_selection_split <- "train"
  }
  gsva_select_expr <- expr_splits[[gsva_selection_split]]
  gsva_select_clin <- clinical_splits[[gsva_selection_split]]
  
  bc_gene_sets <- lapply(bc_gene_sets, function(genes) {
    intersect(genes, rownames(gsva_select_expr))
  })
  bc_gene_sets <- bc_gene_sets[lengths(bc_gene_sets) >= 10]
  if (!length(bc_gene_sets)) {
    message("No MSigDB gene sets overlap with the selection split; will retry after fallback.")
  }
  
  select_ids <- colnames(gsva_select_expr)
  select_stage <- gsva_select_clin$tumor_stage[match(select_ids, gsva_select_clin$patient_id)]
  stage_factor <- factor(select_stage)
  
  if (any(is.na(stage_factor)) || length(unique(stage_factor)) < 2) {
    if (gsva_selection_split != "train") {
      gsva_selection_split <- "train"
      gsva_select_expr <- expr_splits[[gsva_selection_split]]
      gsva_select_clin <- clinical_splits[[gsva_selection_split]]
      select_ids <- colnames(gsva_select_expr)
      select_stage <- gsva_select_clin$tumor_stage[match(select_ids, gsva_select_clin$patient_id)]
      stage_factor <- factor(select_stage)
    }
    if (any(is.na(stage_factor)) || length(unique(stage_factor)) < 2) {
      stop("Tumor stage labels are missing or have fewer than two groups; cannot rank pathways.")
    }
  }
  
  bc_gene_sets <- lapply(bc_gene_sets, function(genes) {
    intersect(genes, rownames(gsva_select_expr))
  })
  bc_gene_sets <- bc_gene_sets[lengths(bc_gene_sets) >= 10]
  
  if (!length(bc_gene_sets)) {
    stop("No MSigDB gene sets overlap with the expression matrix after filtering.")
  }
  
  gsva_rank_param <- GSVA::gsvaParam(
    exprData = gsva_select_expr,
    geneSets = bc_gene_sets,
    minSize = 10,
    maxSize = 500,
    kcdf = "Gaussian",
    verbose = FALSE
  )
  gsva_rank_scores <- GSVA::gsva(gsva_rank_param)
  
  mean_early <- matrixStats::rowMeans2(gsva_rank_scores[, stage_factor == "early", drop = FALSE])
  mean_late <- matrixStats::rowMeans2(gsva_rank_scores[, stage_factor == "late", drop = FALSE])
  delta <- mean_late - mean_early
  
  pvals <- apply(gsva_rank_scores, 1, function(x) {
    stats::t.test(x[stage_factor == "late"], x[stage_factor == "early"])$p.value
  })
  
  gsva_rank_tbl <- data.frame(
    ID = rownames(gsva_rank_scores),
    mean_early = mean_early,
    mean_late = mean_late,
    delta = delta,
    p_value = pvals,
    p_adjust = stats::p.adjust(pvals, method = "BH"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(p_adjust, desc(abs(delta)))
  
  gsva_top_tbl <- gsva_rank_tbl %>%
    dplyr::slice_head(n = 10)
  
  print("Top GSVA-ranked pathways:")
  print(gsva_top_tbl)
  
  gsva_plot_df <- gsva_rank_tbl %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::mutate(ID = stats::reorder(ID, delta))
  
  gsva_rank_plot <- ggplot(gsva_plot_df, aes(x = delta, y = ID, size = abs(delta), color = -log10(p_adjust))) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#97C4BC", high = "#1F78B4", name = "-log10 adj p") +
    labs(
      title = "Breast cancer relevant GSVA ranking",
      subtitle = "Top pathways by GSVA score contrast (late vs early)",
      x = "GSVA score delta (late - early)",
      y = NULL,
      size = "Absolute delta"
    ) +
    theme_nature()
  gsva_rank_plot
  
  ggsave(
    filename = file.path(results_dir, "[c]gsva_rank_plot.png"),
    plot = gsva_rank_plot,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  # Use pre-selected pathway if provided, otherwise select top pathway from GSVA ranking
  if (!is.null(selected_pathway_id)) {
    message(sprintf("Using pre-selected unified pathway: %s", selected_pathway_id))
    selected_pathway_label <- selected_pathway_id
  } else {
    selected_pathway <- gsva_rank_tbl %>% dplyr::slice_head(n = 1)
    
    if (!nrow(selected_pathway)) {
      stop("No pathways returned by GSVA ranking; cannot derive GSVA treatment.")
    }
    
    selected_pathway_id <- selected_pathway$ID[1]
    selected_pathway_label <- selected_pathway_id
    message(sprintf("Selected pathway for downstream GSVA-based treatment: %s", selected_pathway_label))
  }
  
  write.csv(gsva_rank_tbl, file.path(results_dir, "gsva_rank_table.csv"), row.names = FALSE)
  write.csv(gsva_top_tbl, file.path(results_dir, "gsva_top_table.csv"), row.names = FALSE)
  write.csv(selected_pathway, file.path(results_dir, "gsva_selected_pathway.csv"), row.names = FALSE)
  
  selected_gene_set <- bc_terms %>%
    dplyr::filter(gs_name == selected_pathway_id) %>%
    dplyr::pull(gene_symbol) %>%
    unique()
  
  selected_gene_set <- Reduce(
    intersect,
    list(selected_gene_set, rownames(expr_splits$train), rownames(expr_splits$val), rownames(expr_splits$test))
  )
  
  if (length(selected_gene_set) < 10) {
    stop("Selected pathway has fewer than 10 genes present in the expression matrix.")
  }
  
  pathway_gene_sets <- list(selected_pathway_id = selected_gene_set)
  
  # check which of the top DEGs are in the selected pathway and plot them
  # deg_in_pathway <- deg_limma_selected %>%
  #   dplyr::filter(gene %in% selected_gene_set) %>%
  #   dplyr::arrange(adj.P.Val)
  # print("DEGs in selected pathway:")
  # print(deg_in_pathway)
  # 
  # # plot the intersection of top DEGs and selected pathway genes using venn diagram
  # venn_data <- list(
  #   DEGs = deg_limma_selected$gene,
  #   Selected_Pathway = selected_gene_set
  # )
  # 
  # venn_plot <- ggvenn::ggvenn(
  #   venn_data,
  #   c("DEGs", "Selected_Pathway"),
  #   fill_color = c("#D55E00", "#0072B2"),
  #   stroke_size = 0.5,
  #   set_name_size = 4
  # ) +
  #   labs(
  #     title = "Overlap between DEGs and selected pathway genes",
  #     subtitle = sprintf("Pathway: %s", selected_pathway_label)
  #   ) +
  #   theme_nature()
  # venn_plot
  
  
    
  #################################################################################
  # Pathway activity scoring (GSVA) and covariate assembly ----
  ################################################################################
  
  gsva_inputs <- list(
    train = expr_splits$train,
    val = expr_splits$val,
    test = expr_splits$test
  )
  
  gsva_scores <- lapply(gsva_inputs, function(mat) {
    param <- GSVA::gsvaParam(
      exprData = mat,
      geneSets = pathway_gene_sets,
      minSize = 10,
      maxSize = 500,
      kcdf = "Gaussian",
      verbose = FALSE
    )
    GSVA::gsva(param)
  })
  
  pathway_scores <- lapply(names(gsva_scores), function(split_name) {
    score_mat <- gsva_scores[[split_name]]
    data.frame(
      patient_id = colnames(score_mat),
      pathway_score = as.numeric(score_mat[1, ]),
      split = split_name,
      stringsAsFactors = FALSE
    )
  })
  names(pathway_scores) <- names(gsva_scores)
  
  train_pathway_df <- pathway_scores$train %>%
    dplyr::left_join(clinical_splits$train %>% dplyr::select(patient_id, tumor_stage), by = "patient_id")
  
  gsva_density_plot <- ggplot(train_pathway_df, aes(x = pathway_score, fill = tumor_stage)) +
    geom_density(alpha = 0.65) +
    labs(
      title = sprintf("GSVA distribution for %s", selected_pathway_label),
      subtitle = "Treatment D (GSVA score)",
      x = "GSVA score",
      y = "Density",
      fill = "Tumor stage"
    ) +
    scale_fill_manual(values = c("early" = "#63A375", "late" = "#E4572E")) +
    theme_nature()
  gsva_density_plot
  
  ggsave(
    filename = file.path(results_dir, "[d]gsva_density_plot.png"),
    plot = gsva_density_plot,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  #################################################################################
  # GSVA Score Heatmap (samples × GSVA score with tumor stage annotation)
  #################################################################################
  # Combine all splits for heatmap
  all_pathway_scores <- dplyr::bind_rows(pathway_scores) %>%
    dplyr::left_join(
      dplyr::bind_rows(
        clinical_splits$train %>% dplyr::select(patient_id, tumor_stage) %>% dplyr::mutate(split = "train"),
        clinical_splits$val %>% dplyr::select(patient_id, tumor_stage) %>% dplyr::mutate(split = "val"),
        clinical_splits$test %>% dplyr::select(patient_id, tumor_stage) %>% dplyr::mutate(split = "test")
      ),
      by = "patient_id"
    ) %>%
    dplyr::arrange(tumor_stage, dplyr::desc(pathway_score)) %>%
    dplyr::mutate(sample_order = row_number())
  
  gsva_heatmap <- ggplot(all_pathway_scores, aes(x = sample_order, y = 1, fill = pathway_score)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", name = "GSVA\nScore") +
    facet_grid(. ~ tumor_stage, scales = "free_x", space = "free_x") +
    labs(
      title = sprintf("GSVA Score Heatmap: %s", selected_pathway_label),
      subtitle = "Samples ordered by GSVA score within tumor stage",
      x = "Sample",
      y = NULL
    ) +
    theme_nature() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(results_dir, "[d2]gsva_heatmap.png"),
    plot = gsva_heatmap,
    width = 12,
    height = 4,
    dpi = 600
  )
  
  # GSVA heatmap by split
  gsva_heatmap_by_split <- ggplot(
    all_pathway_scores %>% dplyr::mutate(split = factor(split.y, levels = c("train", "val", "test"))),
    aes(x = reorder(patient_id, pathway_score), y = split, fill = pathway_score)
  ) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", name = "GSVA\nScore") +
    facet_grid(. ~ tumor_stage, scales = "free_x", space = "free_x") +
    labs(
      title = sprintf("GSVA Score Heatmap by Split: %s", selected_pathway_label),
      x = "Sample (ordered by GSVA score)",
      y = "Data Split"
    ) +
    theme_nature() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(
    filename = file.path(results_dir, "[d3]gsva_heatmap_by_split.png"),
    plot = gsva_heatmap_by_split,
    width = 14,
    height = 5,
    dpi = 600
  )
  
  nonsynonymous_classes <- c(
    "Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site",
    "Translation_Start_Site", "In_Frame_Del", "In_Frame_Ins"
  )
  
  
  # prepare mutation and methylation summary statistics for downstream analyses
  mutation_summary <- mt %>%
    dplyr::mutate(
      patient_id = Tumor_Sample_Barcode,
      Variant_Classification = as.character(Variant_Classification),
      t_alt_count = suppressWarnings(as.numeric(t_alt_count)),
      t_ref_count = suppressWarnings(as.numeric(t_ref_count))
    ) %>%
    dplyr::filter(patient_id %in% dt2$patient_id) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(
      mutation_burden_total = dplyr::n(),
      mutation_nonsynonymous = sum(Variant_Classification %in% nonsynonymous_classes, na.rm = TRUE),
      mutation_missense_frac = mean(Variant_Classification == "Missense_Mutation", na.rm = TRUE),
      mean_alt_depth = mean(t_alt_count, na.rm = TRUE),
      mean_total_depth = mean(t_alt_count + t_ref_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mutation_missense_frac = ifelse(is.nan(mutation_missense_frac), NA_real_, mutation_missense_frac),
      mutation_nonsynonymous_frac = mutation_nonsynonymous / pmax(mutation_burden_total, 1)
    )
  
  mp_sample_cols <- intersect(colnames(mp)[-1], dt2$patient_id)
  mp_subset <- as.matrix(mp[, ..mp_sample_cols])
  rownames(mp_subset) <- mp$Hugo_Symbol
  
  methylation_summary <- data.frame(
    patient_id = mp_sample_cols,
    methylation_mean = colMeans(mp_subset, na.rm = TRUE),
    methylation_sd = matrixStats::colSds(mp_subset, na.rm = TRUE),
    methylation_hyper_frac = colMeans(mp_subset > 0.75, na.rm = TRUE),
    row.names = NULL
  )
  
  assemble_split_df <- function(clin_df, split_name) {
    clin_df %>%
      dplyr::left_join(pathway_scores[[split_name]] %>% dplyr::select(-split), by = "patient_id") %>%
      dplyr::left_join(mutation_summary, by = "patient_id") %>%
      dplyr::left_join(methylation_summary, by = "patient_id") %>%
      dplyr::mutate(
        tumor_stage_binary = dplyr::case_when(
          tumor_stage == "late" ~ 1,
          tumor_stage == "early" ~ 0,
          TRUE ~ NA_real_
        ),
        split = split_name
      ) %>%
      dplyr::filter(!is.na(tumor_stage_binary), !is.na(pathway_score))
  }
  
  train_df <- assemble_split_df(clinical_splits$train, "train")
  val_df <- assemble_split_df(clinical_splits$val, "val")
  test_df <- assemble_split_df(clinical_splits$test, "test")
  trainval_df <- dplyr::bind_rows(train_df, val_df)
  full_df <- dplyr::bind_rows(train_df, val_df, test_df)
  
  save_violin_plot <- function(df, split_label, results_dir, pathway_label) {
    if (!nrow(df) || !"tumor_stage" %in% names(df)) {
      return(NULL)
    }
    df$tumor_stage <- factor(df$tumor_stage, levels = c("early", "late"))
    p <- ggplot(df, aes(x = tumor_stage, y = pathway_score, fill = tumor_stage)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
      labs(
        title = sprintf("GSVA score by tumor stage (%s)", split_label),
        subtitle = sprintf("Pathway: %s", pathway_label),
        x = "Tumor stage",
        y = "GSVA score"
      ) +
      scale_fill_manual(values = c("early" = "#63A375", "late" = "#E4572E")) +
      theme_nature()
    if (length(unique(stats::na.omit(df$tumor_stage))) >= 2) {
      p <- p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif")
    }
    ggsave(
      filename = file.path(results_dir, sprintf("[v]%s_gsva_violin.png", split_label)),
      plot = p,
      width = 7,
      height = 5,
      dpi = 600
    )
    p
  }
  
  save_violin_plot(train_df, "train", results_dir, selected_pathway_label)
  save_violin_plot(val_df, "val", results_dir, selected_pathway_label)
  save_violin_plot(test_df, "test", results_dir, selected_pathway_label)
  
  adjustment_covariates <- c(
    "age_at_diagnosis",
    "tumor_size",
    "er_status",
    "pr_status",
    "her2_status",
    "neoplasm_histologic_grade",
    "lymph_nodes_examined_positive",
    "nottingham_prognostic_index",
    "tmb_nonsynonymous",
    "mutation_burden_total",
    "mutation_nonsynonymous_frac",
    "mutation_missense_frac",
    "methylation_mean",
    "methylation_sd",
    "methylation_hyper_frac"
  )
  
  assumption_summary <- train_df %>%
    dplyr::group_by(tumor_stage) %>%
    dplyr::summarise(
      mean_D = mean(pathway_score),
      sd_D = sd(pathway_score),
      min_D = min(pathway_score),
      max_D = max(pathway_score),
      n = dplyr::n(),
      .groups = "drop"
    )
  print("GSVA treatment summary by tumor stage (train split):")
  print(assumption_summary)
  
  summarize_stage_gsva <- function(df, split_label) {
    if (!nrow(df)) {
      return(dplyr::tibble())
    }
    df %>%
      dplyr::group_by(tumor_stage) %>%
      dplyr::summarise(
        mean_D = mean(pathway_score),
        sd_D = sd(pathway_score),
        min_D = min(pathway_score),
        max_D = max(pathway_score),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(split = split_label, .before = 1)
  }
  
  gsva_stage_summaries <- dplyr::bind_rows(
    summarize_stage_gsva(train_df, "train"),
    summarize_stage_gsva(val_df, "val"),
    summarize_stage_gsva(test_df, "test")
  )
  print("GSVA treatment summary by split (descriptive stage-wise means):")
  print(gsva_stage_summaries)
  # These split-specific GSVA contrasts are descriptive balance checks, not true causal treatment effects or significance tests--the DR estimators below address those.
  
  covariate_association <- lapply(adjustment_covariates, function(var) {
    values <- train_df[[var]]
    if (is.null(values)) {
      return(data.frame(variable = var, metric = NA_character_, value = NA_real_))
    }
    if (is.numeric(values) || is.integer(values)) {
      stat <- suppressWarnings(cor(train_df$pathway_score, values, use = "pairwise.complete.obs"))
      data.frame(variable = var, metric = "Pearson r", value = stat)
    } else {
      fit <- tryCatch(stats::lm(train_df$pathway_score ~ factor(values)), error = function(e) NULL)
      stat <- if (is.null(fit)) NA_real_ else summary(fit)$r.squared
      data.frame(variable = var, metric = "R^2", value = stat)
    }
  }) %>% dplyr::bind_rows()
  
  print("Association between GSVA treatment and covariates:")
  print(covariate_association)
  
  
  
  #################################################################################
  # Analysis methods ----
  ################################################################################
  
  
  set.seed(202501)
  
  cov_formula <- stats::as.formula(paste("~ 0 +", paste(adjustment_covariates, collapse = " + ")))
  run_split_dml <- function(df, split_label) {
    if (!nrow(df)) {
      return(dplyr::tibble())
    }
    cov_imp <- fit_imputer(df, adjustment_covariates)
    cov_df <- apply_imputer(df, cov_imp)
    cov_matrix <- model.matrix(cov_formula, data = cov_df)
  
    glmnet_fit <- dml_partial_linear(
      Y = df$tumor_stage_binary,
      D = df$pathway_score,
      X = cov_matrix,
      outcome_learner = glmnet_learner(family = "binomial", alpha = 0.5),
      treatment_learner = glmnet_learner(family = "gaussian", alpha = 0.5),
      K = 5,
      seed = 202501
    )
  
    rf_fit <- dml_partial_linear(
      Y = df$tumor_stage_binary,
      D = df$pathway_score,
      X = cov_df,
      outcome_learner = ranger_learner(task = "classification", num.trees = 1200),
      treatment_learner = ranger_learner(task = "regression", num.trees = 1200),
      K = 5,
      seed = 202501
    )
    bart_tbl <- bart_dml_partial(
      Y = df$tumor_stage_binary,
      D = df$pathway_score,
      X = cov_matrix,
      split_label = split_label,
      K = 5,
      seed = 202501
    )
  
    dplyr::bind_rows(
      dplyr::tibble(
        split = split_label,
        model = c("DR-GLMNET", "DR-RF"),
        theta = c(glmnet_fit$theta, rf_fit$theta),
        se = c(glmnet_fit$se, rf_fit$se),
        ci_lower = c(glmnet_fit$ci[1], rf_fit$ci[1]),
        ci_upper = c(glmnet_fit$ci[2], rf_fit$ci[2])
      ),
      bart_tbl
    )
  }
  
  dml_split_tbl <- dplyr::bind_rows(
    run_split_dml(train_df, "train"),
    run_split_dml(val_df, "val"),
    run_split_dml(test_df, "test")
  )
  
  print("Doubly robust estimates for GSVA treatment effect on tumor stage (all splits):")
  print(dml_split_tbl)
  # Best practice: When sharing split-wise DR effects, spell out the nuisance learners and fold structure so readers grasp how cross-fitting protected against overfitting.
  
  dml_full_tbl <- run_split_dml(full_df, "all")
  print("Doubly robust estimates when pooling train/val/test (exploratory full-data fit):")
  print(dml_full_tbl)
  # Best practice: Clearly label pooled fits as exploratory because no holdout remains—regulators expect confirmatory claims to rely on untouched test data.
  
  dr_effects_tbl <- dml_split_tbl %>% dplyr::filter(split == "train")
  
  create_stratified_folds <- function(y, k = 5, seed = 202501) {
    set.seed(seed)
    folds <- integer(length(y))
    classes <- unique(y)
    for (cls in classes) {
      idx <- which(y == cls)
      folds[idx] <- sample(rep_len(1:k, length(idx)))
    }
    folds
  }
  
  cov_imputer_trainval <- fit_imputer(trainval_df, adjustment_covariates)
  trainval_cov_df <- apply_imputer(trainval_df, cov_imputer_trainval)
  X_trainval_cov <- model.matrix(cov_formula, data = trainval_cov_df)
  X_trainval_full <- cbind(pathway_score = trainval_df$pathway_score, X_trainval_cov)
  
  foldid <- create_stratified_folds(trainval_df$tumor_stage_binary, k = 5, seed = 202501)
  glmnet_classifier <- glmnet::cv.glmnet(
    X_trainval_full,
    trainval_df$tumor_stage_binary,
    family = "binomial",
    alpha = 0.5,
    foldid = foldid
  )
  
  test_cov_df <- apply_imputer(test_df, cov_imputer_trainval)
  X_test_cov <- model.matrix(cov_formula, data = test_cov_df)
  X_test_full <- cbind(pathway_score = test_df$pathway_score, X_test_cov)
  glmnet_test_pred <- as.numeric(stats::predict(glmnet_classifier, newx = X_test_full, s = "lambda.min", type = "response"))
  
  glmnet_auc <- compute_auc(glmnet_test_pred, test_df$tumor_stage_binary)
  glmnet_cindex <- binary_c_index(glmnet_test_pred, test_df$tumor_stage_binary)
  
  rf_features <- c("pathway_score", adjustment_covariates)
  rf_imputer <- fit_imputer(trainval_df, rf_features)
  rf_train_df <- apply_imputer(trainval_df, rf_imputer)
  rf_train_df$tumor_stage_binary <- factor(trainval_df$tumor_stage_binary, levels = c(0, 1))
  
  rf_classifier <- ranger::ranger(
    tumor_stage_binary ~ .,
    data = rf_train_df,
    num.trees = 1500,
    mtry = max(1, floor(sqrt(ncol(rf_train_df) - 1))),
    probability = TRUE,
    respect.unordered.factors = "order",
    seed = 202501
  )
  
  rf_test_df <- apply_imputer(test_df, rf_imputer)
  rf_test_pred_mat <- predict(rf_classifier, data = rf_test_df)$predictions
  rf_test_pred <- if (is.matrix(rf_test_pred_mat)) {
    as.numeric(rf_test_pred_mat[, "1", drop = TRUE])
  } else {
    as.numeric(rf_test_pred_mat)
  }
  
  rf_auc <- compute_auc(rf_test_pred, test_df$tumor_stage_binary)
  rf_cindex <- binary_c_index(rf_test_pred, test_df$tumor_stage_binary)
  
  # Bayesian Additive Regression Trees provide a flexible, uncertainty-aware classifier we can
  # compare against the discriminative DR models for predictive sensitivity analysis.
  bart_seed <- 202501
  set.seed(bart_seed)
  bart_fit <- BART::pbart(
    x.train = X_trainval_full,
    y.train = trainval_df$tumor_stage_binary,
    x.test = X_test_full,
    ntree = 200,
    ndpost = 1500,
    nskip = 500,
    usequants = TRUE,
    sparse = TRUE
  )
  bart_test_pred <- as.numeric(bart_fit$prob.test.mean)
  bart_auc <- compute_auc(bart_test_pred, test_df$tumor_stage_binary)
  bart_cindex <- binary_c_index(bart_test_pred, test_df$tumor_stage_binary)
  
  model_metrics <- dplyr::tibble(
    model = c("DR-GLMNET", "DR-RF", "DR-BART"),
    test_auc = c(glmnet_auc, rf_auc, bart_auc),
    test_cindex = c(glmnet_cindex, rf_cindex, bart_cindex)
  )
  
  print("Test-set predictive performance (AUC and concordance index):")
  print(model_metrics)
  
  #################################################################################
  # Comprehensive test metrics for all models
  #################################################################################
  glmnet_metrics <- compute_classification_metrics(glmnet_test_pred, test_df$tumor_stage_binary)
  rf_metrics <- compute_classification_metrics(rf_test_pred, test_df$tumor_stage_binary)
  bart_metrics <- compute_classification_metrics(bart_test_pred, test_df$tumor_stage_binary)
  
  comprehensive_metrics <- dplyr::bind_rows(
    glmnet_metrics %>% dplyr::mutate(model = "DR-GLMNET"),
    rf_metrics %>% dplyr::mutate(model = "DR-RF"),
    bart_metrics %>% dplyr::mutate(model = "DR-BART")
  ) %>%
    dplyr::select(model, dplyr::everything())
  
  write.csv(comprehensive_metrics, file.path(results_dir, "comprehensive_test_metrics.csv"), row.names = FALSE)
  
  #################################################################################
  # Optimal Threshold Analysis (Youden's J and PR-optimal)
  #################################################################################
  
  # Compute validation predictions (needed for threshold optimization)
  cal_val_glmnet_early <- apply_imputer(val_df, cov_imputer_trainval)
  X_val_cov_early <- model.matrix(cov_formula, data = cal_val_glmnet_early)
  X_val_full_early <- cbind(pathway_score = val_df$pathway_score, X_val_cov_early)
  glmnet_val_pred_early <- as.numeric(stats::predict(glmnet_classifier, newx = X_val_full_early, s = "lambda.min", type = "response"))
  
  rf_val_df_early <- apply_imputer(val_df, rf_imputer)
  rf_val_pred_mat_early <- predict(rf_classifier, data = rf_val_df_early)$predictions
  rf_val_pred_early <- if (is.matrix(rf_val_pred_mat_early)) as.numeric(rf_val_pred_mat_early[, "1"]) else as.numeric(rf_val_pred_mat_early)
  
  bart_val_pred_early <- as.numeric(stats::predict(bart_fit, newdata = X_val_full_early)$prob.test.mean)
  
  # Compute optimal thresholds for each model using validation set
  # GLMNET thresholds
  glmnet_youden_val <- compute_youden_threshold(glmnet_val_pred_early, val_df$tumor_stage_binary)
  glmnet_pr_opt_val <- compute_pr_optimal_threshold(glmnet_val_pred_early, val_df$tumor_stage_binary)
  glmnet_youden_test <- compute_youden_threshold(glmnet_test_pred, test_df$tumor_stage_binary)
  glmnet_pr_opt_test <- compute_pr_optimal_threshold(glmnet_test_pred, test_df$tumor_stage_binary)
  
  # RF thresholds
  rf_youden_val <- compute_youden_threshold(rf_val_pred_early, val_df$tumor_stage_binary)
  rf_pr_opt_val <- compute_pr_optimal_threshold(rf_val_pred_early, val_df$tumor_stage_binary)
  rf_youden_test <- compute_youden_threshold(rf_test_pred, test_df$tumor_stage_binary)
  rf_pr_opt_test <- compute_pr_optimal_threshold(rf_test_pred, test_df$tumor_stage_binary)
  
  # BART thresholds
  bart_youden_val <- compute_youden_threshold(bart_val_pred_early, val_df$tumor_stage_binary)
  bart_pr_opt_val <- compute_pr_optimal_threshold(bart_val_pred_early, val_df$tumor_stage_binary)
  bart_youden_test <- compute_youden_threshold(bart_test_pred, test_df$tumor_stage_binary)
  bart_pr_opt_test <- compute_pr_optimal_threshold(bart_test_pred, test_df$tumor_stage_binary)
  
  # Compute metrics at each threshold for validation set
  val_metrics_all_thresholds <- dplyr::bind_rows(
    # Standard 0.5 threshold
    compute_metrics_at_threshold(glmnet_val_pred_early, val_df$tumor_stage_binary, 0.5, "Standard (0.5)") %>% dplyr::mutate(model = "DR-GLMNET", split = "Validation"),
    compute_metrics_at_threshold(rf_val_pred_early, val_df$tumor_stage_binary, 0.5, "Standard (0.5)") %>% dplyr::mutate(model = "DR-RF", split = "Validation"),
    compute_metrics_at_threshold(bart_val_pred_early, val_df$tumor_stage_binary, 0.5, "Standard (0.5)") %>% dplyr::mutate(model = "DR-BART", split = "Validation"),
    # Youden's J threshold
    compute_metrics_at_threshold(glmnet_val_pred_early, val_df$tumor_stage_binary, glmnet_youden_val$threshold, "Youden") %>% dplyr::mutate(model = "DR-GLMNET", split = "Validation"),
    compute_metrics_at_threshold(rf_val_pred_early, val_df$tumor_stage_binary, rf_youden_val$threshold, "Youden") %>% dplyr::mutate(model = "DR-RF", split = "Validation"),
    compute_metrics_at_threshold(bart_val_pred_early, val_df$tumor_stage_binary, bart_youden_val$threshold, "Youden") %>% dplyr::mutate(model = "DR-BART", split = "Validation"),
    # PR-optimal (F1-optimal) threshold
    compute_metrics_at_threshold(glmnet_val_pred_early, val_df$tumor_stage_binary, glmnet_pr_opt_val$threshold, "PR-Optimal") %>% dplyr::mutate(model = "DR-GLMNET", split = "Validation"),
    compute_metrics_at_threshold(rf_val_pred_early, val_df$tumor_stage_binary, rf_pr_opt_val$threshold, "PR-Optimal") %>% dplyr::mutate(model = "DR-RF", split = "Validation"),
    compute_metrics_at_threshold(bart_val_pred_early, val_df$tumor_stage_binary, bart_pr_opt_val$threshold, "PR-Optimal") %>% dplyr::mutate(model = "DR-BART", split = "Validation")
  )
  
  # Compute metrics at each threshold for test set
  test_metrics_all_thresholds <- dplyr::bind_rows(
    # Standard 0.5 threshold
    compute_metrics_at_threshold(glmnet_test_pred, test_df$tumor_stage_binary, 0.5, "Standard (0.5)") %>% dplyr::mutate(model = "DR-GLMNET", split = "Test"),
    compute_metrics_at_threshold(rf_test_pred, test_df$tumor_stage_binary, 0.5, "Standard (0.5)") %>% dplyr::mutate(model = "DR-RF", split = "Test"),
    compute_metrics_at_threshold(bart_test_pred, test_df$tumor_stage_binary, 0.5, "Standard (0.5)") %>% dplyr::mutate(model = "DR-BART", split = "Test"),
    # Youden's J threshold (using val-optimized threshold)
    compute_metrics_at_threshold(glmnet_test_pred, test_df$tumor_stage_binary, glmnet_youden_val$threshold, "Youden (Val-tuned)") %>% dplyr::mutate(model = "DR-GLMNET", split = "Test"),
    compute_metrics_at_threshold(rf_test_pred, test_df$tumor_stage_binary, rf_youden_val$threshold, "Youden (Val-tuned)") %>% dplyr::mutate(model = "DR-RF", split = "Test"),
    compute_metrics_at_threshold(bart_test_pred, test_df$tumor_stage_binary, bart_youden_val$threshold, "Youden (Val-tuned)") %>% dplyr::mutate(model = "DR-BART", split = "Test"),
    # PR-optimal threshold (using val-optimized threshold)
    compute_metrics_at_threshold(glmnet_test_pred, test_df$tumor_stage_binary, glmnet_pr_opt_val$threshold, "PR-Optimal (Val-tuned)") %>% dplyr::mutate(model = "DR-GLMNET", split = "Test"),
    compute_metrics_at_threshold(rf_test_pred, test_df$tumor_stage_binary, rf_pr_opt_val$threshold, "PR-Optimal (Val-tuned)") %>% dplyr::mutate(model = "DR-RF", split = "Test"),
    compute_metrics_at_threshold(bart_test_pred, test_df$tumor_stage_binary, bart_pr_opt_val$threshold, "PR-Optimal (Val-tuned)") %>% dplyr::mutate(model = "DR-BART", split = "Test"),
    # Test-set optimal thresholds (oracle, for reference only)
    compute_metrics_at_threshold(glmnet_test_pred, test_df$tumor_stage_binary, glmnet_youden_test$threshold, "Youden (Test-Oracle)") %>% dplyr::mutate(model = "DR-GLMNET", split = "Test"),
    compute_metrics_at_threshold(rf_test_pred, test_df$tumor_stage_binary, rf_youden_test$threshold, "Youden (Test-Oracle)") %>% dplyr::mutate(model = "DR-RF", split = "Test"),
    compute_metrics_at_threshold(bart_test_pred, test_df$tumor_stage_binary, bart_youden_test$threshold, "Youden (Test-Oracle)") %>% dplyr::mutate(model = "DR-BART", split = "Test")
  )
  
  all_threshold_metrics <- dplyr::bind_rows(val_metrics_all_thresholds, test_metrics_all_thresholds)
  write.csv(all_threshold_metrics, file.path(results_dir, "metrics_by_threshold.csv"), row.names = FALSE)
  
  # Save threshold summary
  threshold_summary <- data.frame(
    model = rep(c("DR-GLMNET", "DR-RF", "DR-BART"), each = 4),
    threshold_type = rep(c("Youden-Val", "Youden-Test", "PR-Opt-Val", "PR-Opt-Test"), 3),
    threshold = c(
      glmnet_youden_val$threshold, glmnet_youden_test$threshold, glmnet_pr_opt_val$threshold, glmnet_pr_opt_test$threshold,
      rf_youden_val$threshold, rf_youden_test$threshold, rf_pr_opt_val$threshold, rf_pr_opt_test$threshold,
      bart_youden_val$threshold, bart_youden_test$threshold, bart_pr_opt_val$threshold, bart_pr_opt_test$threshold
    ),
    stringsAsFactors = FALSE
  )
  write.csv(threshold_summary, file.path(results_dir, "optimal_thresholds.csv"), row.names = FALSE)
  
  # Store predictions for later visualization
  test_predictions <- data.frame(
    patient_id = test_df$patient_id,
    true_label = test_df$tumor_stage_binary,
    glmnet_pred = glmnet_test_pred,
    rf_pred = rf_test_pred,
    bart_pred = bart_test_pred,
    fold_label = fold_label,
    stringsAsFactors = FALSE
  )
  
  # Add validation predictions
  val_predictions <- data.frame(
    patient_id = val_df$patient_id,
    true_label = val_df$tumor_stage_binary,
    glmnet_pred = glmnet_val_pred_early,
    rf_pred = rf_val_pred_early,
    bart_pred = bart_val_pred_early,
    fold_label = fold_label,
    stringsAsFactors = FALSE
  )
  
  write.csv(test_predictions, file.path(results_dir, "test_predictions.csv"), row.names = FALSE)
  write.csv(val_predictions, file.path(results_dir, "val_predictions.csv"), row.names = FALSE)
  
  #################################################################################
  # GSVA Score Ridge Plots (Train/Val/Test)
  #################################################################################
  ridge_data <- dplyr::bind_rows(
    train_df %>% dplyr::select(patient_id, pathway_score, tumor_stage) %>% dplyr::mutate(split = "Train"),
    val_df %>% dplyr::select(patient_id, pathway_score, tumor_stage) %>% dplyr::mutate(split = "Validation"),
    test_df %>% dplyr::select(patient_id, pathway_score, tumor_stage) %>% dplyr::mutate(split = "Test")
  ) %>%
    dplyr::mutate(
      split = factor(split, levels = c("Train", "Validation", "Test")),
      tumor_stage = factor(tumor_stage, levels = c("early", "late"))
    )
  
  ridge_plot <- ggplot(ridge_data, aes(x = pathway_score, y = split, fill = tumor_stage)) +
    ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2) +
    scale_fill_manual(values = c("early" = "#63A375", "late" = "#E4572E")) +
    labs(
      title = sprintf("GSVA Score Distribution: %s", selected_pathway_label),
      subtitle = "Ridge plot by data split and tumor stage",
      x = "GSVA Pathway Score",
      y = "Data Split",
      fill = "Tumor Stage"
    ) +
    theme_nature()
  
  ggsave(
    filename = file.path(results_dir, "[f]gsva_ridge_plot.png"),
    plot = ridge_plot,
    width = 10,
    height = 6,
    dpi = 600
  )
  
  #################################################################################
  # Calibration Curves for Val and Test
  #################################################################################
  cal_val_glmnet <- apply_imputer(val_df, cov_imputer_trainval)
  X_val_cov <- model.matrix(cov_formula, data = cal_val_glmnet)
  X_val_full <- cbind(pathway_score = val_df$pathway_score, X_val_cov)
  glmnet_val_pred <- as.numeric(stats::predict(glmnet_classifier, newx = X_val_full, s = "lambda.min", type = "response"))
  
  rf_val_df <- apply_imputer(val_df, rf_imputer)
  rf_val_pred_mat <- predict(rf_classifier, data = rf_val_df)$predictions
  rf_val_pred <- if (is.matrix(rf_val_pred_mat)) as.numeric(rf_val_pred_mat[, "1"]) else as.numeric(rf_val_pred_mat)
  
  bart_val_pred <- as.numeric(stats::predict(bart_fit, newdata = X_val_full)$prob.test.mean)
  
  # Create calibration data
  cal_test_data <- dplyr::bind_rows(
    create_calibration_data(glmnet_test_pred, test_df$tumor_stage_binary) %>% dplyr::mutate(model = "DR-GLMNET", split = "Test"),
    create_calibration_data(rf_test_pred, test_df$tumor_stage_binary) %>% dplyr::mutate(model = "DR-RF", split = "Test"),
    create_calibration_data(bart_test_pred, test_df$tumor_stage_binary) %>% dplyr::mutate(model = "DR-BART", split = "Test")
  )
  
  cal_val_data <- dplyr::bind_rows(
    create_calibration_data(glmnet_val_pred, val_df$tumor_stage_binary) %>% dplyr::mutate(model = "DR-GLMNET", split = "Validation"),
    create_calibration_data(rf_val_pred, val_df$tumor_stage_binary) %>% dplyr::mutate(model = "DR-RF", split = "Validation"),
    create_calibration_data(bart_val_pred, val_df$tumor_stage_binary) %>% dplyr::mutate(model = "DR-BART", split = "Validation")
  )
  
  cal_all_data <- dplyr::bind_rows(cal_test_data, cal_val_data)
  
  calibration_plot <- ggplot(cal_all_data, aes(x = mean_pred, y = mean_obs, color = model, shape = split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3, alpha = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, mean_obs - 1.96 * se), ymax = pmin(1, mean_obs + 1.96 * se)), width = 0.02, alpha = 0.6) +
    geom_line(aes(group = interaction(model, split)), alpha = 0.5) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "Calibration Curves",
      subtitle = "Predicted probability vs observed proportion",
      x = "Mean Predicted Probability",
      y = "Observed Proportion",
      color = "Model",
      shape = "Split"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_nature()
  
  ggsave(
    filename = file.path(results_dir, "[g]calibration_curves.png"),
    plot = calibration_plot,
    width = 8,
    height = 7,
    dpi = 600
  )
  
  #################################################################################
  # ROC Curves with Confidence Bands
  #################################################################################
  roc_glmnet_test <- pROC::roc(test_df$tumor_stage_binary, glmnet_test_pred, quiet = TRUE)
  roc_rf_test <- pROC::roc(test_df$tumor_stage_binary, rf_test_pred, quiet = TRUE)
  roc_bart_test <- pROC::roc(test_df$tumor_stage_binary, bart_test_pred, quiet = TRUE)
  
  roc_glmnet_val <- pROC::roc(val_df$tumor_stage_binary, glmnet_val_pred, quiet = TRUE)
  roc_rf_val <- pROC::roc(val_df$tumor_stage_binary, rf_val_pred, quiet = TRUE)
  roc_bart_val <- pROC::roc(val_df$tumor_stage_binary, bart_val_pred, quiet = TRUE)
  
  # Function to extract ROC curve data
  roc_to_df <- function(roc_obj, model_name, split_name) {
    data.frame(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      model = model_name,
      split = split_name,
      auc = as.numeric(roc_obj$auc),
      stringsAsFactors = FALSE
    )
  }
  
  roc_data <- dplyr::bind_rows(
    roc_to_df(roc_glmnet_test, "DR-GLMNET", "Test"),
    roc_to_df(roc_rf_test, "DR-RF", "Test"),
    roc_to_df(roc_bart_test, "DR-BART", "Test"),
    roc_to_df(roc_glmnet_val, "DR-GLMNET", "Validation"),
    roc_to_df(roc_rf_val, "DR-RF", "Validation"),
    roc_to_df(roc_bart_val, "DR-BART", "Validation")
  )
  
  # Calculate CI for AUC using bootstrap
  ci_glmnet <- tryCatch(pROC::ci.auc(roc_glmnet_test, method = "bootstrap", boot.n = 500, quiet = TRUE), error = function(e) c(NA, NA, NA))
  ci_rf <- tryCatch(pROC::ci.auc(roc_rf_test, method = "bootstrap", boot.n = 500, quiet = TRUE), error = function(e) c(NA, NA, NA))
  ci_bart <- tryCatch(pROC::ci.auc(roc_bart_test, method = "bootstrap", boot.n = 500, quiet = TRUE), error = function(e) c(NA, NA, NA))
  
  auc_labels <- data.frame(
    model = c("DR-GLMNET", "DR-RF", "DR-BART"),
    auc = c(as.numeric(roc_glmnet_test$auc), as.numeric(roc_rf_test$auc), as.numeric(roc_bart_test$auc)),
    ci_lower = c(ci_glmnet[1], ci_rf[1], ci_bart[1]),
    ci_upper = c(ci_glmnet[3], ci_rf[3], ci_bart[3]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(label = sprintf("%s: %.3f [%.3f-%.3f]", model, auc, ci_lower, ci_upper))
  
  roc_plot <- ggplot(roc_data %>% dplyr::filter(split == "Test"), aes(x = 1 - specificity, y = sensitivity, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(
      values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A"),
      labels = auc_labels$label
    ) +
    labs(
      title = "ROC Curves (Test Set)",
      subtitle = "With bootstrap 95% CI for AUC",
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Model (AUC [95% CI])"
    ) +
    coord_equal() +
    theme_nature() +
    theme(legend.position = c(0.7, 0.25), legend.background = element_rect(fill = "white", color = "gray80"))
  
  ggsave(
    filename = file.path(results_dir, "[h]roc_curves_test.png"),
    plot = roc_plot,
    width = 8,
    height = 7,
    dpi = 600
  )
  
  # Validation ROC
  auc_labels_val <- data.frame(
    model = c("DR-GLMNET", "DR-RF", "DR-BART"),
    auc = c(as.numeric(roc_glmnet_val$auc), as.numeric(roc_rf_val$auc), as.numeric(roc_bart_val$auc)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(label = sprintf("%s: %.3f", model, auc))
  
  roc_plot_val <- ggplot(roc_data %>% dplyr::filter(split == "Validation"), aes(x = 1 - specificity, y = sensitivity, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A"), labels = auc_labels_val$label) +
    labs(title = "ROC Curves (Validation Set)", x = "1 - Specificity", y = "Sensitivity", color = "Model (AUC)") +
    coord_equal() +
    theme_nature() +
    theme(legend.position = c(0.7, 0.25), legend.background = element_rect(fill = "white", color = "gray80"))
  
  ggsave(filename = file.path(results_dir, "[h]roc_curves_val.png"), plot = roc_plot_val, width = 8, height = 7, dpi = 600)
  
  #################################################################################
  # Combined ROC Curves: All 3 Models on Same Plot (Val + Test)
  #################################################################################
  combined_roc_plot <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = model, linetype = split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.0, alpha = 0.9) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    scale_linetype_manual(values = c("Validation" = "dashed", "Test" = "solid")) +
    labs(
      title = sprintf("ROC Curves: All Models (%s)", fold_label),
      subtitle = sprintf("Test AUC - GLMNET: %.3f, RF: %.3f, BART: %.3f",
                        as.numeric(roc_glmnet_test$auc), as.numeric(roc_rf_test$auc), as.numeric(roc_bart_test$auc)),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Model",
      linetype = "Split"
    ) +
    coord_equal() +
    theme_nature()
  
  ggsave(filename = file.path(results_dir, "[h2]roc_curves_combined_3models.png"), plot = combined_roc_plot, width = 10, height = 8, dpi = 600)
  
  #################################################################################
  # PR Curves for All Models
  #################################################################################
  pr_to_df <- function(pred, truth, model_name, split_name) {
    pr_obj <- ROCR::prediction(pred, truth)
    perf <- ROCR::performance(pr_obj, "prec", "rec")
    prauc <- tryCatch({
      as.numeric(ROCR::performance(pr_obj, "aucpr")@y.values[[1]])
    }, error = function(e) NA_real_)
    
    data.frame(
      recall = perf@x.values[[1]],
      precision = perf@y.values[[1]],
      model = model_name,
      split = split_name,
      pr_auc = prauc,
      stringsAsFactors = FALSE
    ) %>% dplyr::filter(!is.na(precision))
  }
  
  pr_data <- dplyr::bind_rows(
    pr_to_df(glmnet_test_pred, test_df$tumor_stage_binary, "DR-GLMNET", "Test"),
    pr_to_df(rf_test_pred, test_df$tumor_stage_binary, "DR-RF", "Test"),
    pr_to_df(bart_test_pred, test_df$tumor_stage_binary, "DR-BART", "Test"),
    pr_to_df(glmnet_val_pred, val_df$tumor_stage_binary, "DR-GLMNET", "Validation"),
    pr_to_df(rf_val_pred, val_df$tumor_stage_binary, "DR-RF", "Validation"),
    pr_to_df(bart_val_pred, val_df$tumor_stage_binary, "DR-BART", "Validation")
  )
  
  # Baseline for PR curve
  baseline_precision <- mean(test_df$tumor_stage_binary)
  
  # Get PR-AUC for each model
  pr_aucs <- pr_data %>% dplyr::group_by(model, split) %>% dplyr::summarise(pr_auc = first(pr_auc), .groups = "drop")
  
  combined_pr_plot <- ggplot(pr_data, aes(x = recall, y = precision, color = model, linetype = split)) +
    geom_hline(yintercept = baseline_precision, linetype = "dotted", color = "gray50") +
    geom_line(linewidth = 1.0, alpha = 0.9) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    scale_linetype_manual(values = c("Validation" = "dashed", "Test" = "solid")) +
    labs(
      title = sprintf("Precision-Recall Curves: All Models (%s)", fold_label),
      subtitle = sprintf("Test PR-AUC - GLMNET: %.3f, RF: %.3f, BART: %.3f",
                        pr_aucs %>% dplyr::filter(model == "DR-GLMNET", split == "Test") %>% dplyr::pull(pr_auc),
                        pr_aucs %>% dplyr::filter(model == "DR-RF", split == "Test") %>% dplyr::pull(pr_auc),
                        pr_aucs %>% dplyr::filter(model == "DR-BART", split == "Test") %>% dplyr::pull(pr_auc)),
      x = "Recall",
      y = "Precision",
      color = "Model",
      linetype = "Split"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_nature()
  
  ggsave(filename = file.path(results_dir, "[h3]pr_curves_combined_3models.png"), plot = combined_pr_plot, width = 10, height = 8, dpi = 600)
  
  #################################################################################
  # Per-Model ROC Curves with Val/Test + Threshold Points
  #################################################################################
  
  # Helper to create single model ROC plot with thresholds
  create_single_model_roc_with_thresholds <- function(roc_val, roc_test, youden_val, youden_test, pr_opt_val, pr_opt_test, model_name) {
    roc_val_df <- data.frame(fpr = 1 - roc_val$specificities, tpr = roc_val$sensitivities, split = "Validation")
    roc_test_df <- data.frame(fpr = 1 - roc_test$specificities, tpr = roc_test$sensitivities, split = "Test")
    roc_combined <- dplyr::bind_rows(roc_val_df, roc_test_df)
    
    # Get coordinates for threshold points
    get_roc_coords <- function(roc_obj, threshold) {
      coords <- tryCatch({
        pROC::coords(roc_obj, threshold, input = "threshold", ret = c("sensitivity", "specificity"))
      }, error = function(e) c(sensitivity = NA, specificity = NA))
      c(fpr = 1 - as.numeric(coords["specificity"]), tpr = as.numeric(coords["sensitivity"]))
    }
    
    threshold_points <- data.frame(
      fpr = c(get_roc_coords(roc_val, youden_val$threshold)["fpr"],
              get_roc_coords(roc_test, youden_val$threshold)["fpr"],  # Apply val threshold to test
              get_roc_coords(roc_val, pr_opt_val$threshold)["fpr"],
              get_roc_coords(roc_test, pr_opt_val$threshold)["fpr"]),
      tpr = c(get_roc_coords(roc_val, youden_val$threshold)["tpr"],
              get_roc_coords(roc_test, youden_val$threshold)["tpr"],
              get_roc_coords(roc_val, pr_opt_val$threshold)["tpr"],
              get_roc_coords(roc_test, pr_opt_val$threshold)["tpr"]),
      threshold_type = c("Youden", "Youden", "PR-Opt", "PR-Opt"),
      split = c("Validation", "Test", "Validation", "Test"),
      threshold = c(youden_val$threshold, youden_val$threshold, pr_opt_val$threshold, pr_opt_val$threshold)
    )
    
    ggplot(roc_combined, aes(x = fpr, y = tpr, color = split)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(linewidth = 1.2, alpha = 0.9) +
      geom_point(data = threshold_points, aes(shape = threshold_type), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
      scale_shape_manual(values = c("Youden" = 17, "PR-Opt" = 15)) +
      labs(
        title = sprintf("ROC Curve with Threshold Points: %s", model_name),
        subtitle = sprintf("Val AUC: %.3f | Test AUC: %.3f | Youden thresh: %.3f | PR-Opt thresh: %.3f",
                          as.numeric(roc_val$auc), as.numeric(roc_test$auc), youden_val$threshold, pr_opt_val$threshold),
        x = "1 - Specificity (False Positive Rate)",
        y = "Sensitivity (True Positive Rate)",
        color = "Split",
        shape = "Threshold"
      ) +
      coord_equal() +
      theme_nature()
  }
  
  # Create per-model ROC plots with thresholds
  glmnet_roc_thresh <- create_single_model_roc_with_thresholds(
    roc_glmnet_val, roc_glmnet_test, glmnet_youden_val, glmnet_youden_test, glmnet_pr_opt_val, glmnet_pr_opt_test, "DR-GLMNET")
  rf_roc_thresh <- create_single_model_roc_with_thresholds(
    roc_rf_val, roc_rf_test, rf_youden_val, rf_youden_test, rf_pr_opt_val, rf_pr_opt_test, "DR-RF")
  bart_roc_thresh <- create_single_model_roc_with_thresholds(
    roc_bart_val, roc_bart_test, bart_youden_val, bart_youden_test, bart_pr_opt_val, bart_pr_opt_test, "DR-BART")
  
  ggsave(filename = file.path(results_dir, "[h4]roc_glmnet_with_thresholds.png"), plot = glmnet_roc_thresh, width = 9, height = 7, dpi = 600)
  ggsave(filename = file.path(results_dir, "[h4]roc_rf_with_thresholds.png"), plot = rf_roc_thresh, width = 9, height = 7, dpi = 600)
  ggsave(filename = file.path(results_dir, "[h4]roc_bart_with_thresholds.png"), plot = bart_roc_thresh, width = 9, height = 7, dpi = 600)
  
  #################################################################################
  # Per-Model PR Curves with Val/Test + Threshold Points
  #################################################################################
  
  create_single_model_pr_with_thresholds <- function(pred_val, truth_val, pred_test, truth_test, 
                                                      youden_val, pr_opt_val, model_name) {
    pr_val_obj <- ROCR::prediction(pred_val, truth_val)
    pr_test_obj <- ROCR::prediction(pred_test, truth_test)
    
    perf_val <- ROCR::performance(pr_val_obj, "prec", "rec")
    perf_test <- ROCR::performance(pr_test_obj, "prec", "rec")
    
    prauc_val <- tryCatch(as.numeric(ROCR::performance(pr_val_obj, "aucpr")@y.values[[1]]), error = function(e) NA)
    prauc_test <- tryCatch(as.numeric(ROCR::performance(pr_test_obj, "aucpr")@y.values[[1]]), error = function(e) NA)
    
    pr_val_df <- data.frame(recall = perf_val@x.values[[1]], precision = perf_val@y.values[[1]], split = "Validation") %>% dplyr::filter(!is.na(precision))
    pr_test_df <- data.frame(recall = perf_test@x.values[[1]], precision = perf_test@y.values[[1]], split = "Test") %>% dplyr::filter(!is.na(precision))
    pr_combined <- dplyr::bind_rows(pr_val_df, pr_test_df)
    
    # Get precision/recall at thresholds
    get_pr_at_thresh <- function(pred, truth, thresh) {
      pred_class <- ifelse(pred >= thresh, 1, 0)
      tp <- sum(pred_class == 1 & truth == 1)
      fp <- sum(pred_class == 1 & truth == 0)
      fn <- sum(pred_class == 0 & truth == 1)
      c(precision = tp / max(tp + fp, 1), recall = tp / max(tp + fn, 1))
    }
    
    threshold_points <- data.frame(
      recall = c(get_pr_at_thresh(pred_val, truth_val, youden_val$threshold)["recall"],
                 get_pr_at_thresh(pred_test, truth_test, youden_val$threshold)["recall"],
                 get_pr_at_thresh(pred_val, truth_val, pr_opt_val$threshold)["recall"],
                 get_pr_at_thresh(pred_test, truth_test, pr_opt_val$threshold)["recall"]),
      precision = c(get_pr_at_thresh(pred_val, truth_val, youden_val$threshold)["precision"],
                    get_pr_at_thresh(pred_test, truth_test, youden_val$threshold)["precision"],
                    get_pr_at_thresh(pred_val, truth_val, pr_opt_val$threshold)["precision"],
                    get_pr_at_thresh(pred_test, truth_test, pr_opt_val$threshold)["precision"]),
      threshold_type = c("Youden", "Youden", "PR-Opt", "PR-Opt"),
      split = c("Validation", "Test", "Validation", "Test")
    )
    
    baseline <- mean(c(truth_val, truth_test))
    
    ggplot(pr_combined, aes(x = recall, y = precision, color = split)) +
      geom_hline(yintercept = baseline, linetype = "dotted", color = "gray50") +
      geom_line(linewidth = 1.2, alpha = 0.9) +
      geom_point(data = threshold_points, aes(shape = threshold_type), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
      scale_shape_manual(values = c("Youden" = 17, "PR-Opt" = 15)) +
      labs(
        title = sprintf("PR Curve with Threshold Points: %s", model_name),
        subtitle = sprintf("Val PR-AUC: %.3f | Test PR-AUC: %.3f", prauc_val, prauc_test),
        x = "Recall",
        y = "Precision",
        color = "Split",
        shape = "Threshold"
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_nature()
  }
  
  glmnet_pr_thresh <- create_single_model_pr_with_thresholds(
    glmnet_val_pred, val_df$tumor_stage_binary, glmnet_test_pred, test_df$tumor_stage_binary, 
    glmnet_youden_val, glmnet_pr_opt_val, "DR-GLMNET")
  rf_pr_thresh <- create_single_model_pr_with_thresholds(
    rf_val_pred, val_df$tumor_stage_binary, rf_test_pred, test_df$tumor_stage_binary,
    rf_youden_val, rf_pr_opt_val, "DR-RF")
  bart_pr_thresh <- create_single_model_pr_with_thresholds(
    bart_val_pred, val_df$tumor_stage_binary, bart_test_pred, test_df$tumor_stage_binary,
    bart_youden_val, bart_pr_opt_val, "DR-BART")
  
  ggsave(filename = file.path(results_dir, "[h5]pr_glmnet_with_thresholds.png"), plot = glmnet_pr_thresh, width = 9, height = 7, dpi = 600)
  ggsave(filename = file.path(results_dir, "[h5]pr_rf_with_thresholds.png"), plot = rf_pr_thresh, width = 9, height = 7, dpi = 600)
  ggsave(filename = file.path(results_dir, "[h5]pr_bart_with_thresholds.png"), plot = bart_pr_thresh, width = 9, height = 7, dpi = 600)

  #################################################################################
  # Predicted vs True Visualization (Test Set)
  #################################################################################
  pred_vs_true_data <- test_predictions %>%
    tidyr::pivot_longer(cols = c(glmnet_pred, rf_pred, bart_pred), names_to = "model", values_to = "predicted") %>%
    dplyr::mutate(
      model = dplyr::case_when(
        model == "glmnet_pred" ~ "DR-GLMNET",
        model == "rf_pred" ~ "DR-RF",
        model == "bart_pred" ~ "DR-BART"
      ),
      true_label = factor(true_label, levels = c(0, 1), labels = c("Early", "Late"))
    )
  
  pred_vs_true_plot <- ggplot(pred_vs_true_data, aes(x = true_label, y = predicted, fill = true_label)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~ model, ncol = 3) +
    scale_fill_manual(values = c("Early" = "#63A375", "Late" = "#E4572E")) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
    labs(
      title = "Predicted Probabilities vs True Labels (Test Set)",
      x = "True Tumor Stage",
      y = "Predicted Probability of Late Stage",
      fill = "True Stage"
    ) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(
    filename = file.path(results_dir, "[i]predicted_vs_true_test.png"),
    plot = pred_vs_true_plot,
    width = 12,
    height = 5,
    dpi = 600
  )
  
  #################################################################################
  # Variable Importance (GLMNET coefficients, RF importance, BART varcount)
  #################################################################################
  # GLMNET coefficients
  glmnet_coefs <- as.matrix(coef(glmnet_classifier, s = "lambda.min"))
  glmnet_importance <- data.frame(
    variable = rownames(glmnet_coefs),
    importance = abs(as.numeric(glmnet_coefs)),
    model = "DR-GLMNET",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(variable != "(Intercept)", importance > 0) %>%
    dplyr::arrange(dplyr::desc(importance)) %>%
    dplyr::slice_head(n = 15)
  
  # RF importance
  rf_importance <- data.frame(
    variable = names(rf_classifier$variable.importance),
    importance = as.numeric(rf_classifier$variable.importance),
    model = "DR-RF",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(importance)) %>%
    dplyr::slice_head(n = 15)
  
  # BART variable usage (varcount)
  bart_varcount <- colMeans(bart_fit$varcount)
  bart_importance <- data.frame(
    variable = colnames(X_trainval_full),
    importance = bart_varcount,
    model = "DR-BART",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(importance)) %>%
    dplyr::slice_head(n = 15)
  
  # Combine and normalize importance
  all_importance <- dplyr::bind_rows(glmnet_importance, rf_importance, bart_importance) %>%
    dplyr::group_by(model) %>%
    dplyr::mutate(importance_scaled = importance / max(importance, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  # Get top variables across all models
  top_vars <- all_importance %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(total_imp = sum(importance_scaled, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_imp)) %>%
    dplyr::slice_head(n = 12) %>%
    dplyr::pull(variable)
  
  importance_plot_data <- all_importance %>%
    dplyr::filter(variable %in% top_vars) %>%
    dplyr::mutate(variable = factor(variable, levels = rev(top_vars)))
  
  importance_plot <- ggplot(importance_plot_data, aes(x = importance_scaled, y = variable, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "Variable Importance Comparison",
      subtitle = "Normalized importance (GLMNET: |coef|, RF: impurity, BART: varcount)",
      x = "Scaled Importance",
      y = NULL,
      fill = "Model"
    ) +
    theme_nature() +
    theme(legend.position = "top")
  
  ggsave(
    filename = file.path(results_dir, "[j]variable_importance_comparison.png"),
    plot = importance_plot,
    width = 10,
    height = 8,
    dpi = 600
  )
  
  #################################################################################
  # Permutation Importance for RF (optional, can be slow)
  #################################################################################
  # Compute permutation importance for RF on test set
  rf_perm_importance <- tryCatch({
    compute_permutation_importance(
      model = rf_classifier,
      X = rf_test_df,
      y = test_df$tumor_stage_binary,
      predict_fn = function(m, x) {
        pred_mat <- predict(m, data = x)$predictions
        if (is.matrix(pred_mat)) as.numeric(pred_mat[, "1"]) else as.numeric(pred_mat)
      },
      n_permutations = 5,
      seed = 202501
    )
  }, error = function(e) NULL)
  
  if (!is.null(rf_perm_importance)) {
    perm_importance_plot <- ggplot(
      rf_perm_importance %>% dplyr::slice_head(n = 15) %>% dplyr::mutate(variable = forcats::fct_reorder(variable, importance)),
      aes(x = importance, y = variable, fill = importance)
    ) +
      geom_col() +
      scale_fill_viridis_c(option = "plasma") +
      labs(
        title = "Permutation Importance (DR-RF, Test Set)",
        subtitle = "Drop in AUC when feature is permuted",
        x = "Importance (AUC drop)",
        y = NULL
      ) +
      theme_nature() +
      theme(legend.position = "none")
    
    ggsave(
      filename = file.path(results_dir, "[k]permutation_importance_rf.png"),
      plot = perm_importance_plot,
      width = 9,
      height = 7,
      dpi = 600
    )
    
    write.csv(rf_perm_importance, file.path(results_dir, "permutation_importance_rf.csv"), row.names = FALSE)
  }
  
  #################################################################################
  # DEG-Pathway Overlap Venn Diagram
  #################################################################################
  deg_genes_set <- unique(deg_limma_selected$gene)
  pathway_genes_set <- unique(selected_gene_set)
  
  venn_data <- list(
    DEGs = deg_genes_set,
    Pathway = pathway_genes_set
  )
  
  venn_plot <- tryCatch({
    ggvenn::ggvenn(
      venn_data,
      fill_color = c("#D55E00", "#0072B2"),
      stroke_size = 0.5,
      set_name_size = 5,
      text_size = 4
    ) +
      labs(
        title = "DEG-Pathway Gene Overlap",
        subtitle = sprintf("Pathway: %s", selected_pathway_label)
      ) +
      theme_nature()
  }, error = function(e) {
    message("Venn diagram failed: ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(venn_plot)) {
    ggsave(
      filename = file.path(results_dir, "[l]deg_pathway_venn.png"),
      plot = venn_plot,
      width = 8,
      height = 7,
      dpi = 600
    )
  }
  
  # Save overlap genes
  overlap_genes <- intersect(deg_genes_set, pathway_genes_set)
  if (length(overlap_genes)) {
    write.csv(data.frame(gene = overlap_genes), file.path(results_dir, "deg_pathway_overlap_genes.csv"), row.names = FALSE)
  }

  reporting_candidates <- c("DR-GLMNET", "DR-RF", "DR-BART")
  best_model <- model_metrics %>%
    dplyr::filter(model %in% reporting_candidates) %>%
    dplyr::arrange(dplyr::desc(test_auc), dplyr::desc(test_cindex)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(model)
  
  best_effect <- dr_effects_tbl %>% dplyr::filter(model == best_model)
  
  message(sprintf("Selected model for reporting (based on test AUC): %s", best_model))
  print(best_effect)
  
  #################################################################################
  # Confusion Matrix for Best Model (at different thresholds)
  #################################################################################
  
  # Get predictions for best model
  best_test_pred <- switch(best_model,
    "DR-GLMNET" = glmnet_test_pred,
    "DR-RF" = rf_test_pred,
    "DR-BART" = bart_test_pred
  )
  best_val_pred <- switch(best_model,
    "DR-GLMNET" = glmnet_val_pred,
    "DR-RF" = rf_val_pred,
    "DR-BART" = bart_val_pred
  )
  best_youden <- switch(best_model,
    "DR-GLMNET" = glmnet_youden_val,
    "DR-RF" = rf_youden_val,
    "DR-BART" = bart_youden_val
  )
  best_pr_opt <- switch(best_model,
    "DR-GLMNET" = glmnet_pr_opt_val,
    "DR-RF" = rf_pr_opt_val,
    "DR-BART" = bart_pr_opt_val
  )
  
  # Confusion matrix at standard threshold (0.5)
  cm_standard <- create_confusion_matrix_plot(best_test_pred, test_df$tumor_stage_binary, 0.5, best_model, "Test", "Standard")
  ggsave(filename = file.path(results_dir, "[m1]confusion_matrix_best_model_standard.png"), plot = cm_standard, width = 6, height = 5, dpi = 600)
  
  # Confusion matrix at Youden's J threshold
  cm_youden <- create_confusion_matrix_plot(best_test_pred, test_df$tumor_stage_binary, best_youden$threshold, best_model, "Test", sprintf("Youden (%.3f)", best_youden$threshold))
  ggsave(filename = file.path(results_dir, "[m2]confusion_matrix_best_model_youden.png"), plot = cm_youden, width = 6, height = 5, dpi = 600)
  
  # Confusion matrix at PR-optimal threshold
  cm_pr_opt <- create_confusion_matrix_plot(best_test_pred, test_df$tumor_stage_binary, best_pr_opt$threshold, best_model, "Test", sprintf("PR-Opt (%.3f)", best_pr_opt$threshold))
  ggsave(filename = file.path(results_dir, "[m3]confusion_matrix_best_model_pr_opt.png"), plot = cm_pr_opt, width = 6, height = 5, dpi = 600)
  
  #################################################################################
  # Reliability Diagram (Calibration Plot) for Best Model
  #################################################################################
  
  best_cal_test <- create_calibration_data(best_test_pred, test_df$tumor_stage_binary, n_bins = 10)
  best_cal_val <- create_calibration_data(best_val_pred, val_df$tumor_stage_binary, n_bins = 10)
  
  best_cal_data <- dplyr::bind_rows(
    best_cal_test %>% dplyr::mutate(split = "Test"),
    best_cal_val %>% dplyr::mutate(split = "Validation")
  )
  
  reliability_diagram <- ggplot(best_cal_data, aes(x = mean_pred, y = mean_obs, color = split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
    geom_point(aes(size = n), alpha = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, mean_obs - 1.96 * se), ymax = pmin(1, mean_obs + 1.96 * se)), width = 0.02, alpha = 0.6) +
    geom_line(alpha = 0.7, linewidth = 1) +
    scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
    scale_size_continuous(range = c(2, 8), name = "N obs") +
    labs(
      title = sprintf("Reliability Diagram: %s", best_model),
      subtitle = "Calibration of predicted probabilities (95% CI)",
      x = "Mean Predicted Probability",
      y = "Observed Proportion (Actual Late Stage)",
      color = "Split"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_nature()
  
  ggsave(filename = file.path(results_dir, "[n]reliability_diagram_best_model.png"), plot = reliability_diagram, width = 8, height = 7, dpi = 600)
  
  # Add histogram of predictions for calibration context
  pred_histogram_data <- data.frame(
    pred = c(best_val_pred, best_test_pred),
    split = c(rep("Validation", length(best_val_pred)), rep("Test", length(best_test_pred))),
    true_label = c(val_df$tumor_stage_binary, test_df$tumor_stage_binary)
  )
  
  pred_histogram <- ggplot(pred_histogram_data, aes(x = pred, fill = factor(true_label))) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    facet_wrap(~ split, ncol = 2) +
    scale_fill_manual(values = c("0" = "#63A375", "1" = "#E4572E"), labels = c("Early", "Late"), name = "True Stage") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    geom_vline(xintercept = best_youden$threshold, linetype = "dotted", color = "blue", linewidth = 0.8) +
    geom_vline(xintercept = best_pr_opt$threshold, linetype = "dotdash", color = "purple", linewidth = 0.8) +
    labs(
      title = sprintf("Prediction Distribution: %s", best_model),
      subtitle = sprintf("Thresholds: 0.5 (black), Youden %.3f (blue), PR-Opt %.3f (purple)",
                        best_youden$threshold, best_pr_opt$threshold),
      x = "Predicted Probability",
      y = "Count"
    ) +
    theme_nature()
  
  ggsave(filename = file.path(results_dir, "[n2]prediction_distribution_best_model.png"), plot = pred_histogram, width = 10, height = 5, dpi = 600)
  
  
  #################################################################################
  # Sensitivity analysis: remove genomic covariates ----
  ################################################################################
  
  genomic_covariates <- c(
    "mutation_burden_total",
    "mutation_nonsynonymous_frac",
    "mutation_missense_frac",
    "methylation_mean",
    "methylation_sd",
    "methylation_hyper_frac"
  )
  clinical_covariates <- setdiff(adjustment_covariates, genomic_covariates)
  clinical_formula <- stats::as.formula(paste("~ 0 +", paste(clinical_covariates, collapse = " + ")))
  
  clinical_imp <- fit_imputer(train_df, clinical_covariates)
  clinical_cov_df <- apply_imputer(train_df, clinical_imp)
  clinical_matrix <- model.matrix(clinical_formula, data = clinical_cov_df)
  
  clinical_glmnet <- dml_partial_linear(
    Y = train_df$tumor_stage_binary,
    D = train_df$pathway_score,
    X = clinical_matrix,
    outcome_learner = glmnet_learner(family = "binomial", alpha = 0.5),
    treatment_learner = glmnet_learner(family = "gaussian", alpha = 0.5),
    K = 5,
    seed = 202501
  )
  
  clinical_rf <- dml_partial_linear(
    Y = train_df$tumor_stage_binary,
    D = train_df$pathway_score,
    X = clinical_cov_df,
    outcome_learner = ranger_learner(task = "classification", num.trees = 1200),
    treatment_learner = ranger_learner(task = "regression", num.trees = 1200),
    K = 5,
    seed = 202501
  )
  
  clinical_sensitivity_tbl <- dplyr::tibble(
    model = c("DR-GLMNET", "DR-RF"),
    theta = c(clinical_glmnet$theta, clinical_rf$theta),
    ci_lower = c(clinical_glmnet$ci[1], clinical_rf$ci[1]),
    ci_upper = c(clinical_glmnet$ci[2], clinical_rf$ci[2])
  )
  print("Sensitivity (clinical-only covariates) DR estimates:")
  print(clinical_sensitivity_tbl)
  
  write.csv(dml_split_tbl, file.path(results_dir, "dml_split_table.csv"), row.names = FALSE)
  write.csv(model_metrics, file.path(results_dir, "test_model_metrics.csv"), row.names = FALSE)
  write.csv(best_effect, file.path(results_dir, "best_model_effect.csv"), row.names = FALSE)
  write.csv(clinical_sensitivity_tbl, file.path(results_dir, "clinical_sensitivity_table.csv"), row.names = FALSE)
  
  saveRDS(
    list(
      fold_label = fold_label,
      selected_pathway = selected_pathway_id,
      dml_split_tbl = dml_split_tbl,
      model_metrics = model_metrics,
      best_effect = best_effect,
      clinical_sensitivity_tbl = clinical_sensitivity_tbl
    ),
    file.path(results_dir, "fold_results.rds")
  )
  
    list(
      fold_label = fold_label,
      selected_pathway = selected_pathway_id,
      dml_split_tbl = dml_split_tbl,
      model_metrics = model_metrics,
      comprehensive_metrics = comprehensive_metrics,
      all_threshold_metrics = all_threshold_metrics,
      test_predictions = test_predictions,
      val_predictions = val_predictions,
      best_model = best_model,
      best_effect = best_effect,
      clinical_sensitivity_tbl = clinical_sensitivity_tbl,
      threshold_summary = threshold_summary,
      roc_auc_test = c("DR-GLMNET" = as.numeric(roc_glmnet_test$auc), 
                       "DR-RF" = as.numeric(roc_rf_test$auc), 
                       "DR-BART" = as.numeric(roc_bart_test$auc)),
      roc_auc_val = c("DR-GLMNET" = as.numeric(roc_glmnet_val$auc), 
                      "DR-RF" = as.numeric(roc_rf_val$auc), 
                      "DR-BART" = as.numeric(roc_bart_val$auc)),
      pr_auc_test = setNames(
        sapply(list(glmnet_test_pred, rf_test_pred, bart_test_pred), function(p) {
          tryCatch(as.numeric(ROCR::performance(ROCR::prediction(p, test_df$tumor_stage_binary), "aucpr")@y.values[[1]]), error = function(e) NA)
        }),
        c("DR-GLMNET", "DR-RF", "DR-BART")
      ),
      pr_auc_val = setNames(
        sapply(list(glmnet_val_pred, rf_val_pred, bart_val_pred), function(p) {
          tryCatch(as.numeric(ROCR::performance(ROCR::prediction(p, val_df$tumor_stage_binary), "aucpr")@y.values[[1]]), error = function(e) NA)
        }),
        c("DR-GLMNET", "DR-RF", "DR-BART")
      )
    )
}

dir.create("results", "final_results", showWarnings = FALSE)

#################################################################################
# Phase 1: Collect top pathways from all folds for unified pathway selection
#################################################################################

cat("\n=== Phase 1: Collecting top pathways from all folds ===\n\n")
cat(sprintf(
  "Outer folds: %d | Inner folds per outer: %d | Total configs: %d\n",
  K_outer, K_inner, total_fold_configs
))

all_top_pathways <- list()
pathway_collection_errors <- list()

get_inner_folds_for_outer <- function(outer_fold, k_outer, k_inner) {
  remaining <- setdiff(seq_len(k_outer), outer_fold)
  if (k_inner > length(remaining)) {
    stop(sprintf(
      "K_inner (%d) exceeds available folds (%d) for outer fold %d.",
      k_inner, length(remaining), outer_fold
    ))
  }
  remaining[seq_len(k_inner)]
}

for (outer_fold in seq_len(K_outer)) {
  inner_folds <- get_inner_folds_for_outer(outer_fold, K_outer, K_inner)
  
  for (inner_fold in inner_folds) {
    train_folds <- setdiff(seq_len(K_outer), c(outer_fold, inner_fold))
    train_ids <- unlist(cv_folds[train_folds], use.names = FALSE)
    val_ids <- cv_folds[[inner_fold]]
    test_ids <- cv_folds[[outer_fold]]
    
    if (!length(train_ids) || !length(val_ids) || !length(test_ids)) {
      next
    }
    
    fold_label <- sprintf(
      "Train%s_Val%d_Test%d",
      paste(train_folds, collapse = ""),
      inner_fold,
      outer_fold
    )
    
    pathway_res <- tryCatch(
      get_gsva_top_pathways(train_ids, val_ids, test_ids, fold_label, gsva_selection_split = "val", n_top = 20),
      error = function(e) {
        pathway_collection_errors[[fold_label]] <<- conditionMessage(e)
        NULL
      }
    )
    
    if (!is.null(pathway_res)) {
      all_top_pathways[[fold_label]] <- pathway_res
    }
  }
}

if (length(pathway_collection_errors)) {
  cat("\nWarning: Some folds failed during pathway collection:\n")
  print(pathway_collection_errors)
}

# Perform rank-weighted voting to select unified pathway
cat("\n=== Performing rank-weighted voting for unified pathway ===\n\n")

if (length(all_top_pathways) == 0) {
  stop("No pathways were collected from any fold. Cannot proceed.")
}

# Collect all top pathways with fold-specific ranks
all_pathways_by_fold <- lapply(names(all_top_pathways), function(fold_label) {
  top_paths <- all_top_pathways[[fold_label]]$top_pathways
  data.frame(
    fold = fold_label,
    pathway = top_paths,
    rank = seq_along(top_paths),
    max_rank = length(top_paths),
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

# Aggregate votes and rank scores (higher score = higher rank across folds)
pathway_votes_df <- all_pathways_by_fold %>%
  dplyr::mutate(rank_score = max_rank - rank + 1) %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(
    votes = dplyr::n(),
    rank_score = sum(rank_score),
    mean_rank = mean(rank),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(rank_score), dplyr::desc(votes), mean_rank, pathway)

cat("\nPathway voting results (top 20):\n")
print(head(pathway_votes_df, 20))

# Select the pathway with the best rank-weighted score
unified_pathway <- pathway_votes_df$pathway[1]
cat(sprintf("\n*** UNIFIED PATHWAY SELECTED: %s (rank_score: %d; votes: %d/%d folds) ***\n\n",
            unified_pathway, pathway_votes_df$rank_score[1], pathway_votes_df$votes[1], length(all_top_pathways)))

# Save pathway selection results
pathway_selection_dir <- file.path("results", "final_results", "summary_nested_cv")
dir.create(pathway_selection_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  pathway_votes_df,
  file.path(pathway_selection_dir, "pathway_voting_results.csv"),
  row.names = FALSE
)

# Save top pathways from each fold (define before using in plots)

write.csv(
  all_pathways_by_fold,
  file.path(pathway_selection_dir, "all-20-selected-pathways-based-on-internal-validation-folds.csv"),
  row.names = FALSE
)

#################################################################################
# Pathway Vote Bar Chart
#################################################################################
pathway_vote_plot <- ggplot(
  pathway_votes_df %>% dplyr::slice_head(n = 25) %>% dplyr::mutate(pathway = forcats::fct_reorder(pathway, votes)),
  aes(x = votes, y = pathway, fill = votes)
) +
  geom_col() +
  geom_vline(xintercept = pathway_votes_df$votes[1], linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    title = "Pathway Voting Results Across Folds",
    subtitle = sprintf("Selected pathway: %s (votes: %d)", unified_pathway, pathway_votes_df$votes[1]),
    x = "Number of Votes (Top-20 appearances across folds)",
    y = NULL,
    fill = "Votes"
  ) +
  theme_nature() +
  theme(legend.position = "none")

ggsave(
  filename = file.path(pathway_selection_dir, "pathway_vote_bar_chart.png"),
  plot = pathway_vote_plot,
  width = 12,
  height = 10,
  dpi = 600
)

#################################################################################
# Pathway Selection Heatmap (Fold × Pathway)
#################################################################################
# Create presence/rank matrix
pathway_presence_matrix <- all_pathways_by_fold %>%
  dplyr::mutate(present = 1) %>%
  tidyr::pivot_wider(
    id_cols = fold,
    names_from = pathway,
    values_from = rank,
    values_fill = NA
  ) %>%
  tibble::column_to_rownames("fold")

# Filter to top pathways only
top_pathways_to_show <- pathway_votes_df %>% dplyr::slice_head(n = 20) %>% dplyr::pull(pathway)
pathway_presence_matrix <- pathway_presence_matrix[, intersect(top_pathways_to_show, colnames(pathway_presence_matrix)), drop = FALSE]

pathway_heatmap_df <- all_pathways_by_fold %>%
  dplyr::filter(pathway %in% top_pathways_to_show) %>%
  dplyr::mutate(
    pathway = factor(pathway, levels = rev(top_pathways_to_show)),
    inv_rank = 21 - rank  # higher = better rank
  )

pathway_heatmap_plot <- ggplot(pathway_heatmap_df, aes(x = fold, y = pathway, fill = inv_rank)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", na.value = "grey90", name = "Rank\n(20=top)") +
  labs(
    title = "Pathway Ranking Heatmap Across Folds",
    subtitle = "Cell value indicates rank (darker = higher rank)",
    x = "Fold",
    y = "Pathway"
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(pathway_selection_dir, "pathway_selection_heatmap.png"),
  plot = pathway_heatmap_plot,
  width = 16,
  height = 10,
  dpi = 600
)

#################################################################################
# Phase 2: Run full analysis with unified pathway
#################################################################################

cat("\n=== Phase 2: Running full analysis with unified pathway ===\n\n")
cat(sprintf(
  "Outer folds: %d | Inner folds per outer: %d | Total configs: %d\n",
  K_outer, K_inner, total_fold_configs
))

all_fold_results <- list()
all_dml_results <- list()
all_model_metrics <- list()
all_comprehensive_metrics <- list()
all_threshold_metrics <- list()
all_test_predictions <- list()
fold_errors <- list()

for (outer_fold in seq_len(K_outer)) {
  inner_folds <- get_inner_folds_for_outer(outer_fold, K_outer, K_inner)
  inner_effects <- list()
  inner_dirs <- character(0)

  for (inner_fold in inner_folds) {
    train_folds <- setdiff(seq_len(K_outer), c(outer_fold, inner_fold))
    train_ids <- unlist(cv_folds[train_folds], use.names = FALSE)
    val_ids <- cv_folds[[inner_fold]]
    test_ids <- cv_folds[[outer_fold]]

    if (!length(train_ids) || !length(val_ids) || !length(test_ids)) {
      next
    }

    fold_label <- sprintf(
      "Train%s_Val%d_Test%d",
      paste(train_folds, collapse = ""),
      inner_fold,
      outer_fold
    )
    results_dir <- file.path("results", "final_results", fold_label)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    res <- tryCatch(
      run_fold(train_ids, val_ids, test_ids, fold_label, results_dir, 
               gsva_selection_split = "val", selected_pathway_id = unified_pathway),
      error = function(e) {
        fold_errors[[fold_label]] <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(res)) next

    all_fold_results[[fold_label]] <- res
    all_dml_results[[fold_label]] <- res$dml_split_tbl %>%
      dplyr::mutate(outer_fold = outer_fold, inner_fold = inner_fold, fold_label = fold_label)
    all_model_metrics[[fold_label]] <- res$model_metrics %>%
      dplyr::mutate(outer_fold = outer_fold, inner_fold = inner_fold, fold_label = fold_label)
    
    # Collect comprehensive metrics and predictions
    if (!is.null(res$comprehensive_metrics)) {
      all_comprehensive_metrics[[fold_label]] <- res$comprehensive_metrics %>%
        dplyr::mutate(outer_fold = outer_fold, inner_fold = inner_fold, fold_label = fold_label)
    }
    if (!is.null(res$all_threshold_metrics)) {
      all_threshold_metrics[[fold_label]] <- res$all_threshold_metrics %>%
        dplyr::mutate(outer_fold = outer_fold, inner_fold = inner_fold, fold_label = fold_label)
    }
    if (!is.null(res$test_predictions)) {
      all_test_predictions[[fold_label]] <- res$test_predictions
    }

    inner_effects[[fold_label]] <- res$dml_split_tbl %>%
      dplyr::filter(split == "train") %>%
      dplyr::mutate(inner_fold = inner_fold, fold_label = fold_label)
    inner_dirs <- c(inner_dirs, results_dir)
  }

  inner_effects_df <- dplyr::bind_rows(inner_effects)
  if (nrow(inner_effects_df)) {
    inner_effects_plot <- ggplot(
      inner_effects_df,
      aes(x = model, y = theta, fill = factor(inner_fold))
    ) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_errorbar(
        aes(ymin = theta - 1.96 * se, ymax = theta + 1.96 * se),
        position = position_dodge(width = 0.7),
        width = 0.2
      ) +
      labs(
        title = sprintf("Inner-fold causal effects (outer fold %d)", outer_fold),
        x = "Model",
        y = "Estimated effect (theta)",
        fill = "Inner fold"
      ) +
      theme_nature()

    for (dir_path in inner_dirs) {
      ggsave(
        filename = file.path(dir_path, "[e]inner_fold_causal_effects.png"),
        plot = inner_effects_plot,
        width = 9,
        height = 6,
        dpi = 600
      )
    }
  }
}

if (length(fold_errors)) {
  write.csv(
    data.frame(fold_label = names(fold_errors), error = unlist(fold_errors)),
    file.path("results", "final_results", "summary_nested_cv", "fold_errors.csv"),
    row.names = FALSE
  )
}

summary_dir <- file.path("results", "final_results", "summary_nested_cv")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

all_dml_df <- dplyr::bind_rows(all_dml_results)
all_model_metrics_df <- dplyr::bind_rows(all_model_metrics)

if (nrow(all_dml_df)) {
  write.csv(all_dml_df, file.path(summary_dir, "all_dml_results.csv"), row.names = FALSE)
}

if (nrow(all_model_metrics_df)) {
  write.csv(all_model_metrics_df, file.path(summary_dir, "all_model_metrics.csv"), row.names = FALSE)
}

if (nrow(all_model_metrics_df)) {
  model_metric_summary <- all_model_metrics_df %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      mean_auc = mean(test_auc, na.rm = TRUE),
      sd_auc = sd(test_auc, na.rm = TRUE),
      mean_cindex = mean(test_cindex, na.rm = TRUE),
      sd_cindex = sd(test_cindex, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(model_metric_summary, file.path(summary_dir, "model_metrics_summary.csv"), row.names = FALSE)
}

#################################################################################
# Comprehensive Summary Visualizations
#################################################################################

cat("\n=== Generating comprehensive summary plots ===\n\n")

# Collect all comprehensive metrics
all_comprehensive_metrics_df <- dplyr::bind_rows(all_comprehensive_metrics)
all_threshold_metrics_df <- dplyr::bind_rows(all_threshold_metrics)
all_test_predictions_df <- dplyr::bind_rows(all_test_predictions)

if (nrow(all_comprehensive_metrics_df)) {
  write.csv(all_comprehensive_metrics_df, file.path(summary_dir, "all_comprehensive_test_metrics.csv"), row.names = FALSE)
}

if (nrow(all_threshold_metrics_df)) {
  write.csv(all_threshold_metrics_df, file.path(summary_dir, "all_threshold_metrics.csv"), row.names = FALSE)
}

if (nrow(all_test_predictions_df)) {
  write.csv(all_test_predictions_df, file.path(summary_dir, "all_test_predictions.csv"), row.names = FALSE)
}

#################################################################################
# Effect Size Forest Plot (DR theta across all test folds)
#################################################################################
if (nrow(all_dml_df)) {
  # Filter to test split effects for forest plot
  test_effects <- all_dml_df %>%
    dplyr::filter(split == "test") %>%
    dplyr::mutate(
      fold_short = gsub("Train|_Val|_Test", "", fold_label),
      fold_idx = as.numeric(factor(fold_label))
    )
  
  effect_forest_plot <- ggplot(test_effects, aes(x = theta, y = reorder(fold_label, fold_idx), color = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    geom_errorbarh(
      aes(xmin = ci_lower, xmax = ci_upper),
      position = position_dodge(width = 0.6),
      height = 0.3,
      alpha = 0.7
    ) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "Causal Effect Estimates Across All Folds (Forest Plot)",
      subtitle = "Doubly robust theta with 95% CI (test set)",
      x = "Estimated Treatment Effect (theta)",
      y = "Fold",
      color = "Model"
    ) +
    theme_nature() +
    theme(axis.text.y = element_text(size = 6))
  
  ggsave(
    filename = file.path(summary_dir, "effect_size_forest_plot.png"),
    plot = effect_forest_plot,
    width = 14,
    height = 12,
    dpi = 600
  )
  
  #################################################################################
  # Fold-wise Waterfall Plot for Causal Effects
  #################################################################################
  waterfall_data <- test_effects %>%
    dplyr::group_by(model) %>%
    dplyr::arrange(theta) %>%
    dplyr::mutate(rank = row_number()) %>%
    dplyr::ungroup()
  
  waterfall_plot <- ggplot(waterfall_data, aes(x = rank, y = theta, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "Waterfall Plot: Ordered Causal Effect Estimates",
      subtitle = "Each bar represents a fold, ordered by effect magnitude",
      x = "Fold Rank (ordered by theta)",
      y = "Estimated Treatment Effect (theta)",
      fill = "Model"
    ) +
    theme_nature()
  
  ggsave(
    filename = file.path(summary_dir, "causal_effect_waterfall.png"),
    plot = waterfall_plot,
    width = 14,
    height = 7,
    dpi = 600
  )
  
  #################################################################################
  # Causal Effect Violin Plot (Distribution across all test folds)
  #################################################################################
  causal_violin_plot <- ggplot(test_effects, aes(x = model, y = theta, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = sprintf("Distribution of Causal Effect Estimates Across %d Test Folds", total_fold_configs),
      subtitle = "Doubly robust theta estimates (test set)",
      x = "Model",
      y = "Estimated Treatment Effect (theta)",
      fill = "Model"
    ) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(
    filename = file.path(summary_dir, "causal_effect_violin.png"),
    plot = causal_violin_plot,
    width = 10,
    height = 7,
    dpi = 600
  )
}

#################################################################################
# C-index Forest Plot (Model performance across folds)
#################################################################################
if (nrow(all_model_metrics_df)) {
  cindex_forest_plot <- ggplot(
    all_model_metrics_df %>% dplyr::mutate(fold_idx = as.numeric(factor(fold_label))),
    aes(x = test_cindex, y = reorder(fold_label, fold_idx), color = model)
  ) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "C-index Performance Across All Test Folds",
      subtitle = "Concordance index (higher = better discrimination)",
      x = "C-index",
      y = "Fold",
      color = "Model"
    ) +
    theme_nature() +
    theme(axis.text.y = element_text(size = 6))
  
  ggsave(
    filename = file.path(summary_dir, "cindex_forest_plot.png"),
    plot = cindex_forest_plot,
    width = 12,
    height = 12,
    dpi = 600
  )
}

#################################################################################
# Final Comparison: 3 Models × Fold Grouped Bars (Test Fold Performance)
#################################################################################
if (nrow(all_model_metrics_df)) {
  final_comparison_plot <- ggplot(
    all_model_metrics_df %>% dplyr::mutate(fold_idx = as.numeric(factor(fold_label))),
    aes(x = factor(fold_idx), y = test_auc, fill = model)
  ) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8, alpha = 0.85) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = sprintf("Test Set AUC Performance: 3 Models Across %d Folds", total_fold_configs),
      subtitle = "Held-out test fold performance comparison",
      x = "Fold Index",
      y = "Test AUC",
      fill = "Model"
    ) +
    theme_nature() +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(summary_dir, sprintf("final_comparison_%dfolds_grouped_bars.png", total_fold_configs)),
    plot = final_comparison_plot,
    width = 16,
    height = 8,
    dpi = 600
  )
}

#################################################################################
# Comprehensive Metrics Violin Plot (All metrics distribution)
#################################################################################
if (nrow(all_comprehensive_metrics_df)) {
  metrics_long <- all_comprehensive_metrics_df %>%
    tidyr::pivot_longer(
      cols = c(accuracy, balanced_accuracy, precision, recall, sensitivity, specificity, f1, roc_auc, pr_auc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("accuracy", "balanced_accuracy", "precision", "recall", 
                                          "sensitivity", "specificity", "f1", "roc_auc", "pr_auc"))
    )
  
  metrics_violin_plot <- ggplot(metrics_long, aes(x = metric, y = value, fill = model)) +
    geom_violin(alpha = 0.6, position = position_dodge(width = 0.8), trim = FALSE) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), alpha = 0.8, outlier.size = 0.5) +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "Distribution of All Classification Metrics Across 20 Test Folds",
      subtitle = "Comprehensive performance evaluation",
      x = "Metric",
      y = "Value",
      fill = "Model"
    ) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(summary_dir, "all_metrics_violin_plot.png"),
    plot = metrics_violin_plot,
    width = 14,
    height = 8,
    dpi = 600
  )
  
  # Save summary statistics for all metrics
  metrics_summary <- metrics_long %>%
    dplyr::group_by(model, metric) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      n = sum(!is.na(value)),
      .groups = "drop"
    )
  
  write.csv(metrics_summary, file.path(summary_dir, "comprehensive_metrics_summary.csv"), row.names = FALSE)
}

#################################################################################
# GSVA Score Heatmap (if we have predictions with GSVA scores)
#################################################################################
if (nrow(all_test_predictions_df) && "pathway_score" %in% names(all_test_predictions_df)) {
  # Create GSVA heatmap data
  gsva_heatmap_data <- all_test_predictions_df %>%
    dplyr::select(patient_id, fold_label, true_label) %>%
    dplyr::distinct()
  
  # This would need pathway_score in predictions - add if available
}

#################################################################################
# Predicted vs True Summary (Combined across folds)
#################################################################################
if (nrow(all_test_predictions_df)) {
  pred_summary <- all_test_predictions_df %>%
    tidyr::pivot_longer(
      cols = c(glmnet_pred, rf_pred, bart_pred),
      names_to = "model",
      values_to = "predicted"
    ) %>%
    dplyr::mutate(
      model = dplyr::case_when(
        model == "glmnet_pred" ~ "DR-GLMNET",
        model == "rf_pred" ~ "DR-RF",
        model == "bart_pred" ~ "DR-BART"
      ),
      true_label_factor = factor(true_label, levels = c(0, 1), labels = c("Early", "Late"))
    )
  
  pred_vs_true_summary_plot <- ggplot(pred_summary, aes(x = true_label_factor, y = predicted, fill = true_label_factor)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    facet_wrap(~ model, ncol = 3) +
    scale_fill_manual(values = c("Early" = "#63A375", "Late" = "#E4572E")) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
    labs(
      title = "Predicted Probabilities vs True Labels (All Test Folds Combined)",
      subtitle = sprintf("N = %d predictions across %d folds", nrow(pred_summary) / 3, total_fold_configs),
      x = "True Tumor Stage",
      y = "Predicted Probability of Late Stage",
      fill = "True Stage"
    ) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(
    filename = file.path(summary_dir, "predicted_vs_true_all_folds.png"),
    plot = pred_vs_true_summary_plot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

#################################################################################
# Overall Averaged ROC-AUC and PR-AUC Plots with Uncertainty Intervals
#################################################################################

cat("\n=== Generating overall averaged AUC plots across all folds ===\n\n")

# Collect AUC values from all folds
all_auc_data <- lapply(names(all_fold_results), function(fold_name) {
  res <- all_fold_results[[fold_name]]
  if (is.null(res$roc_auc_test) || is.null(res$pr_auc_test)) return(NULL)
  
  data.frame(
    fold = fold_name,
    model = names(res$roc_auc_test),
    roc_auc_test = as.numeric(res$roc_auc_test),
    roc_auc_val = as.numeric(res$roc_auc_val),
    pr_auc_test = as.numeric(res$pr_auc_test),
    pr_auc_val = as.numeric(res$pr_auc_val),
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

if (nrow(all_auc_data)) {
  # Save AUC data
  write.csv(all_auc_data, file.path(summary_dir, "all_auc_values_by_fold.csv"), row.names = FALSE)
  
  # Compute summary statistics for AUC values
  auc_summary <- all_auc_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      n_folds = dplyr::n(),
      # ROC-AUC Test
      roc_auc_test_mean = mean(roc_auc_test, na.rm = TRUE),
      roc_auc_test_sd = sd(roc_auc_test, na.rm = TRUE),
      roc_auc_test_se = sd(roc_auc_test, na.rm = TRUE) / sqrt(sum(!is.na(roc_auc_test))),
      roc_auc_test_ci_lower = roc_auc_test_mean - 1.96 * roc_auc_test_se,
      roc_auc_test_ci_upper = roc_auc_test_mean + 1.96 * roc_auc_test_se,
      # PR-AUC Test
      pr_auc_test_mean = mean(pr_auc_test, na.rm = TRUE),
      pr_auc_test_sd = sd(pr_auc_test, na.rm = TRUE),
      pr_auc_test_se = sd(pr_auc_test, na.rm = TRUE) / sqrt(sum(!is.na(pr_auc_test))),
      pr_auc_test_ci_lower = pr_auc_test_mean - 1.96 * pr_auc_test_se,
      pr_auc_test_ci_upper = pr_auc_test_mean + 1.96 * pr_auc_test_se,
      .groups = "drop"
    )
  
  write.csv(auc_summary, file.path(summary_dir, "auc_summary_statistics.csv"), row.names = FALSE)
  
  # Print summary
  cat("\nROC-AUC Summary (Test Set) across", total_fold_configs, "folds:\n")
  print(auc_summary %>% dplyr::select(model, n_folds, roc_auc_test_mean, roc_auc_test_sd, roc_auc_test_ci_lower, roc_auc_test_ci_upper))
  
  cat("\nPR-AUC Summary (Test Set) across", total_fold_configs, "folds:\n")
  print(auc_summary %>% dplyr::select(model, n_folds, pr_auc_test_mean, pr_auc_test_sd, pr_auc_test_ci_lower, pr_auc_test_ci_upper))
  
  #################################################################################
  # Averaged ROC-AUC Bar Plot with 95% CI
  #################################################################################
  roc_auc_bar_plot <- ggplot(auc_summary, aes(x = model, y = roc_auc_test_mean, fill = model)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_errorbar(aes(ymin = pmax(0, roc_auc_test_ci_lower), ymax = pmin(1, roc_auc_test_ci_upper)), 
                  width = 0.2, linewidth = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = sprintf("Average ROC-AUC Across %d Test Folds", total_fold_configs),
      subtitle = "Mean ± 95% CI (standard error-based)",
      x = "Model",
      y = "ROC-AUC",
      fill = "Model"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(summary_dir, "averaged_roc_auc_bar_plot.png"), plot = roc_auc_bar_plot, width = 8, height = 6, dpi = 600)
  
  #################################################################################
  # Averaged PR-AUC Bar Plot with 95% CI
  #################################################################################
  pr_auc_bar_plot <- ggplot(auc_summary, aes(x = model, y = pr_auc_test_mean, fill = model)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_errorbar(aes(ymin = pmax(0, pr_auc_test_ci_lower), ymax = pmin(1, pr_auc_test_ci_upper)), 
                  width = 0.2, linewidth = 0.8) +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = sprintf("Average PR-AUC Across %d Test Folds", total_fold_configs),
      subtitle = "Mean ± 95% CI (standard error-based)",
      x = "Model",
      y = "PR-AUC",
      fill = "Model"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(summary_dir, "averaged_pr_auc_bar_plot.png"), plot = pr_auc_bar_plot, width = 8, height = 6, dpi = 600)
  
  #################################################################################
  # Combined AUC Forest Plot (ROC and PR-AUC side by side)
  #################################################################################
  auc_forest_data <- auc_summary %>%
    tidyr::pivot_longer(
      cols = c(roc_auc_test_mean, pr_auc_test_mean),
      names_to = "metric_type",
      values_to = "mean_auc"
    ) %>%
    dplyr::mutate(
      ci_lower = dplyr::case_when(
        metric_type == "roc_auc_test_mean" ~ roc_auc_test_ci_lower,
        metric_type == "pr_auc_test_mean" ~ pr_auc_test_ci_lower
      ),
      ci_upper = dplyr::case_when(
        metric_type == "roc_auc_test_mean" ~ roc_auc_test_ci_upper,
        metric_type == "pr_auc_test_mean" ~ pr_auc_test_ci_upper
      ),
      metric_type = dplyr::case_when(
        metric_type == "roc_auc_test_mean" ~ "ROC-AUC",
        metric_type == "pr_auc_test_mean" ~ "PR-AUC"
      )
    )
  
  auc_forest_plot <- ggplot(auc_forest_data, aes(x = mean_auc, y = model, color = model)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 1) +
    facet_wrap(~ metric_type, ncol = 2, scales = "free_x") +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = sprintf("Model Performance Summary: Mean AUC ± 95%% CI (%d Folds)", total_fold_configs),
      x = "AUC Value",
      y = "Model",
      color = "Model"
    ) +
    coord_cartesian(xlim = c(0.3, 1)) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(summary_dir, "auc_forest_plot_summary.png"), plot = auc_forest_plot, width = 12, height = 5, dpi = 600)
  
  #################################################################################
  # AUC Distribution Box Plot (showing variability across folds)
  #################################################################################
  auc_long <- all_auc_data %>%
    tidyr::pivot_longer(
      cols = c(roc_auc_test, pr_auc_test),
      names_to = "metric",
      values_to = "auc"
    ) %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        metric == "roc_auc_test" ~ "ROC-AUC",
        metric == "pr_auc_test" ~ "PR-AUC"
      )
    )
  
  auc_distribution_plot <- ggplot(auc_long, aes(x = model, y = auc, fill = model)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~ metric, ncol = 2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = sprintf("AUC Distribution Across %d Test Folds", total_fold_configs),
      subtitle = "Each point represents one fold's test set performance",
      x = "Model",
      y = "AUC Value",
      fill = "Model"
    ) +
    coord_cartesian(ylim = c(0.3, 1)) +
    theme_nature() +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(summary_dir, "auc_distribution_boxplot.png"), plot = auc_distribution_plot, width = 10, height = 6, dpi = 600)
  
  #################################################################################
  # Pairwise Model Comparison (Fold-by-Fold)
  #################################################################################
  auc_wide <- all_auc_data %>%
    dplyr::select(fold, model, roc_auc_test) %>%
    tidyr::pivot_wider(names_from = model, values_from = roc_auc_test)
  
  if (all(c("DR-GLMNET", "DR-RF", "DR-BART") %in% names(auc_wide))) {
    pairwise_plot <- ggplot(auc_wide, aes(x = `DR-GLMNET`, y = `DR-RF`)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(size = 2.5, alpha = 0.7, color = "#377EB8") +
      labs(
        title = "Fold-by-Fold ROC-AUC Comparison: GLMNET vs RF",
        subtitle = "Points above diagonal indicate RF outperforms GLMNET",
        x = "DR-GLMNET ROC-AUC",
        y = "DR-RF ROC-AUC"
      ) +
      coord_equal(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
      theme_nature()
    
    ggsave(filename = file.path(summary_dir, "pairwise_auc_glmnet_vs_rf.png"), plot = pairwise_plot, width = 7, height = 7, dpi = 600)
  }
}

cat("\n=== Summary plots saved to:", summary_dir, "===\n\n")


#################################################################################
# UNIFIED PATHWAY SELECTION APPROACH (Updated Jan 2026)
#################################################################################
# This analysis now uses a unified pathway selection strategy across all folds:
# 
# Phase 1: Pathway Collection & Voting
# - For each fold configuration, extract the top 20 pathways based on GSVA ranking
#   using the internal validation set
# - Perform unweighted voting across all fold configurations
# - Select the pathway with the most votes as the unified pathway
# - Results saved to: results/summary_nested_cv/pathway_voting_results.csv
#
# Phase 2: Full Analysis with Unified Pathway  
# - Run complete analysis pipeline using the unified pathway across all folds
# - This ensures consistent pathway-based treatment definition across all splits
# - Enables more robust and comparable causal effect estimates
#
# Key changes:
# - Added get_gsva_top_pathways() function to extract pathways without full analysis
# - Modified run_fold() to accept optional selected_pathway_id parameter
# - Split main execution into two phases for pathway selection and full analysis
#################################################################################

#################################################################################
# Q: Describe the model fitting approach described in this pipeline. What are the key components and steps involved in the analysis?
# A: In this pipeline, the model fitting approach involves several key components and steps to analyze the relationship between gene expression data and breast cancer tumor stages. The main steps include: 
# 1. Data Preprocessing: The gene expression data is preprocessed using z-score normalization to standardize the expression levels across samples. This step ensures that the data is on a comparable scale for downstream analyses.
# 2. Differential Expression Analysis: The limma-voom method is employed to identify differentially expressed genes (DEGs) between early and late tumor stages. Genes are filtered based on expression levels, and statistical tests are performed to determine significant DEGs.
# 3. GSVA-driven pathway ranking: GSVA scores are computed for MSigDB gene sets using GSVA's built-in ranking,
#    and pathways are ranked by the GSVA score contrast between tumor stages.
# 4. Pathway Activity Scoring: The GSVA method is used to compute pathway activity scores for each sample based
#    on the selected pathway from the GSVA ranking. This step quantifies the activity of specific pathways
#    associated with tumor stages.
# 5. Covariate Assembly: Clinical and genomic covariates are assembled alongside the pathway activity scores. This includes variables such as age, tumor size, hormone receptor status, mutation burden, and methylation levels.
# 6. Doubly Robust Estimation: The DML (Doubly Machine Learning) approach is applied to estimate the causal effect of pathway activity on tumor stage while adjusting for covariates. Different machine learning models (e.g., GLMNET, Random Forest, BART) are used as nuisance learners to model the outcome and treatment. 
# 7. Sensitivity Analysis: A sensitivity analysis is conducted by removing genomic covariates and re-estimating the treatment effect using only clinical covariates. This step assesses the robustness of the findings to different sets of covariates. 
# Overall, the pipeline integrates statistical modeling, machine learning, and biological interpretation to investigate the relationship between gene expression, pathway activity, and breast cancer tumor stages.

# Explain the data splitting strategy used in this analysis. How are the training, validation, and test sets defined and utilized?
# A: The data splitting strategy in this analysis involves dividing the dataset into three distinct subsets: training, validation, and test sets.
# 1. Training Set: The training set is used to fit the models and estimate the parameters. It includes a portion of the data that is representative of the overall dataset. The training set is utilized for tasks such as differential expression analysis, fitting machine learning models for outcome and treatment, and estimating the causal effects using the DML approach.
# 2. Validation Set: The validation set is used to tune model hyperparameters and select the best-performing models. It serves as an intermediate evaluation step to prevent overfitting to the training data. The validation set is employed to assess the performance of different models (e.g., GLM NET, Random Forest) and select the optimal model based on metrics such as AUC and concordance index. 
# 3. Test Set: The test set is reserved for the final evaluation of the selected models. It is not used during the model fitting or selection process, ensuring an unbiased assessment of model performance. The test set is utilized to compute predictive performance metrics (AUC and concordance index) for the final models and to report the doubly robust estimates of the treatment effect on tumor stage. 
# The splitting is typically done randomly while ensuring that the distribution of key variables (e.g., tumor stage) is maintained across the subsets. Each subset is used for specific purposes in the analysis pipeline, allowing for robust model training, validation, and testing to ensure reliable and generalizable results. All models fitted adopt cross-fitting to mitigate overfitting and enhance causal effect estimation. Stratified folds are created based on the binary tumor stage outcome to ensure balanced representation in each fold during cross-validation. This approach helps to ensure that the models are trained and evaluated on diverse subsets of the data, leading to more robust and reliable findings.
# 4. Reporting and Interpretation: The results from each split are reported separately, allowing for a comprehensive understanding of the model's performance and the estimated treatment effects across different data subsets. This approach provides insights into the generalizability of the findings and helps to identify any potential discrepancies between the splits.
# 5. Exploratory Full-Data Fit: An additional exploratory analysis is conducted by pooling all data (train, validation, and test) to fit the models and estimate treatment effects. This pooled analysis is clearly labeled as exploratory since it does not maintain an untouched test set for confirmatory claims. The results from this full-data fit provide additional insights but should be interpreted with caution due to the lack of a separate test set.
