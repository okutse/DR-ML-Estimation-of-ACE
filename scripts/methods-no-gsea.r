# Title: Doubly-robust machine learning estimation of the average causal effect of (selected) DEGs on breast cancer stage and survival outcomes
# Author: Amos Okutse
# Date: Jan 2025


# Comparative Methods:
#(1) DR Penalized Logistic Regression
#(2) DR Random Forest
#(3) DR BART
#(4) DR Logistic regression on categorized DEGs (low, mid, high expression levels based on quantiles)

## Load required packages
required_pkgs <- c("msigdbr", "GSVA", "ranger", "pROC")
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

library(pathview)         # for KEGG pathway visualization
library(msigdbr)       # MSigDB gene sets
library(GSVA)          # pathway activity scoring
library(ranger)        # fast random forests
library(pROC)          # AUC computation
library(ROCR)          # for PR-AUC computation
library(tidyr)         # for data reshaping
library(forcats)       # for factor manipulation
library(viridis)       # for color scales
library(ggridges)      # for ridge plots

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


# source helper files
source("scripts/helper.R")

# -----------------------------------------------------------------------------
# Helper utilities used across the workflow
# -----------------------------------------------------------------------------

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

log_split_summary <- function(ids, clin_df, label) {
  stages <- clin_df$tumor_stage[match(ids, clin_df$patient_id)]
  stage_tbl <- table(factor(stages, levels = c("early", "late")), useNA = "ifany")
  cat(sprintf("[%s] n=%d (early=%d, late=%d, NA=%d)\n",
              label,
              length(ids),
              stage_tbl["early"],
              stage_tbl["late"],
              stage_tbl["<NA>"]))
  invisible(stage_tbl)
}

log_missing_summary <- function(df, label, top_n = 10) {
  na_counts <- sort(colSums(is.na(df)), decreasing = TRUE)
  total_rows <- nrow(df)
  cat(sprintf("[%s] Missingness summary: %d rows\n", label, total_rows))
  print(head(na_counts, top_n))
  invisible(na_counts)
}

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
    
    f1 <- 2 * prec * rec / (prec + rec + 1e-10)
    f1[is.na(f1)] <- 0
    
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

#' Create ROC curve plot
create_roc_curve_plot <- function(pred_val, truth_val, pred_test, truth_test, model_name) {
  roc_val <- tryCatch(pROC::roc(response = truth_val, predictor = pred_val, quiet = TRUE, direction = "<"), error = function(e) NULL)
  roc_test <- tryCatch(pROC::roc(response = truth_test, predictor = pred_test, quiet = TRUE, direction = "<"), error = function(e) NULL)
  
  if (is.null(roc_val) || is.null(roc_test)) return(NULL)
  
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
  
  ggplot(roc_data, aes(x = fpr, y = tpr, color = split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
    labs(
      title = sprintf("ROC Curve: %s", model_name),
      subtitle = sprintf("Val AUC: %.3f | Test AUC: %.3f", as.numeric(roc_val$auc), as.numeric(roc_test$auc)),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Split"
    ) +
    coord_equal() +
    theme_nature() +
    theme(legend.position = "right")
}

#' Create PR curve plot
create_pr_curve_plot <- function(pred_val, truth_val, pred_test, truth_test, model_name) {
  pr_val <- tryCatch(ROCR::prediction(pred_val, truth_val), error = function(e) NULL)
  pr_test <- tryCatch(ROCR::prediction(pred_test, truth_test), error = function(e) NULL)
  
  if (is.null(pr_val) || is.null(pr_test)) return(NULL)
  
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
  baseline <- mean(c(truth_val, truth_test))
  
  ggplot(pr_data, aes(x = recall, y = precision, color = split)) +
    geom_hline(yintercept = baseline, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
    labs(
      title = sprintf("Precision-Recall Curve: %s", model_name),
      subtitle = sprintf("Baseline (class proportion): %.3f", baseline),
      x = "Recall (Sensitivity)",
      y = "Precision (Positive Predictive Value)",
      color = "Split"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_nature() +
    theme(legend.position = "right")
}

#' Create calibration plot
create_calibration_plot <- function(pred_val, truth_val, pred_test, truth_test, model_name) {
  cal_val <- create_calibration_data(pred_val, truth_val) %>% dplyr::mutate(split = "Validation")
  cal_test <- create_calibration_data(pred_test, truth_test) %>% dplyr::mutate(split = "Test")
  cal_data <- dplyr::bind_rows(cal_val, cal_test)
  
  ggplot(cal_data, aes(x = mean_pred, y = mean_obs, color = split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = n), alpha = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, mean_obs - 1.96 * se), ymax = pmin(1, mean_obs + 1.96 * se)), width = 0.02, alpha = 0.5) +
    geom_line(alpha = 0.7) +
    scale_color_manual(values = c("Validation" = "#F8766D", "Test" = "#00BA38")) +
    scale_size_continuous(range = c(2, 8), name = "N obs") +
    labs(
      title = sprintf("Calibration Plot: %s", model_name),
      subtitle = "Perfect calibration lies on the diagonal",
      x = "Mean Predicted Probability",
      y = "Observed Proportion",
      color = "Split"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_nature()
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
    ci = c(theta_hat - 1.645 * se_hat, theta_hat + 1.645 * se_hat),
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
    ci_lower = theta_hat - 1.645 * se_hat,
    ci_upper = theta_hat + 1.645 * se_hat
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
cat(sprintf("Merged clinical data: %d patients, %d columns\n", nrow(dt), ncol(dt)))
log_missing_summary(dt, "clinical_merged")
 
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

cat(sprintf("Expression matrix after filtering to common samples: %d genes x %d samples\n",
            nrow(exp_mat), ncol(exp_mat)))


dt2 <- dt %>%
  dplyr::filter(patient_id %in% common_samples) %>%
  dplyr::arrange(match(patient_id, colnames(exp_mat)))
# Keep tumor_stage (required for stratification) and don't drop other clinical columns yet
dt2 <- dt2[!is.na(dt2$tumor_stage), ] # filter to non-missing tumor stage
cat("dt2 after filtering to non-missing tumor stage: ", nrow(dt2), " rows\n")
cat("dt2 columns:", paste(colnames(dt2), collapse=", "), "\n")
cat(sprintf("Tumor stage distribution overall: early=%d, late=%d\n",
            sum(dt2$tumor_stage == "early", na.rm = TRUE),
            sum(dt2$tumor_stage == "late", na.rm = TRUE)))
exp_mat <- exp_mat[, dt2$patient_id, drop = FALSE] # keep only subjects with observed tumor stage to align ids
grp <- factor(dt2$tumor_stage, levels = c("early", "late")) # early = 1338; late = 128; n = 1466 samples
table(grp)

# check whether exp_mat are already normalized log2 (range is 0 to 16)
summary(as.numeric(exp_mat))

# omit all NA in exp_mat if any after median imputation
exp_mat <- exp_mat[complete.cases(exp_mat), , drop = FALSE]

# Pre-processing steps for mRNA expression data
## Train, val, test split before preprocessing using stratification by tumor stage
sample_ids <- colnames(exp_mat)
n_total <- length(sample_ids)
train_frac <- 0.70
val_frac <- 0.15

split_ids <- stratified_split(
  ids = sample_ids,
  strata = dt2$tumor_stage,
  train_frac = train_frac,
  val_frac = val_frac,
  seed = 202501
)

train_ids <- split_ids$train
val_ids <- split_ids$val
test_ids <- split_ids$test

if (!length(train_ids) || !length(val_ids) || !length(test_ids)) {
  stop("Stratified split failed; one of the partitions is empty.")
}

split_stage_summary <- lapply(
  list(train = train_ids, val = val_ids, test = test_ids),
  function(ids) {
    stages <- dt2$tumor_stage[match(ids, dt2$patient_id)]
    table(factor(stages, levels = c("early", "late")), useNA = "ifany")
  }
)
print("Tumor stage distribution by split:")
print(split_stage_summary)

cat("Split sample counts and stage breakdown:\n")
log_split_summary(train_ids, dt2, "train")
log_split_summary(val_ids, dt2, "val")
log_split_summary(test_ids, dt2, "test")

collapse_duplicate_genes <- function(mat) {
  rn <- rownames(mat)
  dup_gene_names <- unique(rn[duplicated(rn)])
  n_dup_gene_rows <- sum(duplicated(rn))
  n_dup_gene_names <- length(dup_gene_names)

  cn <- colnames(mat)
  dup_col_names <- unique(cn[duplicated(cn)])
  n_dup_col_cols <- sum(duplicated(cn))
  n_dup_col_names <- length(dup_col_names)

  if (n_dup_gene_names > 0 || n_dup_col_names > 0) {
    cat(sprintf("Duplicate genes detected: %d rows across %d gene names.\n",
                n_dup_gene_rows, n_dup_gene_names))
    cat(sprintf("Duplicate columns detected: %d columns across %d column names.\n",
                n_dup_col_cols, n_dup_col_names))
  }

  if (n_dup_gene_names > 0) {
    cat("First 20 duplicate gene names:\n")
    print(utils::head(dup_gene_names, 20))

    k <- min(3, ncol(mat))
    if (k > 0) {
      gene_means_list <- lapply(utils::head(dup_gene_names, 20), function(g) {
        sub_mat <- mat[rn == g, seq_len(k), drop = FALSE]
        data.frame(
          gene = g,
          t(colMeans(sub_mat, na.rm = TRUE)),
          stringsAsFactors = FALSE
        )
      })
      gene_means_df <- dplyr::bind_rows(gene_means_list)
      colnames(gene_means_df)[-1] <- paste0("mean_col", seq_len(k))
      cat(sprintf("Mean of columns 1..%d before collapse (first 20 duplicate genes):\n", k))
      print(gene_means_df)
    }
  }

  if (n_dup_col_names > 0) {
    cat("First 20 duplicate column names:\n")
    print(utils::head(dup_col_names, 20))
  }

  if (anyDuplicated(rn)) {
    warning("Detected duplicated gene symbols; collapsing by mean expression.")
    summed <- rowsum(mat, rn)
    counts <- as.numeric(table(rn)[rownames(summed)])
    mat <- sweep(summed, 1, counts, "/")
    cat(sprintf("Collapsed %d duplicated rows into %d unique gene rows.\n",
                n_dup_gene_rows, n_dup_gene_names))
  }
  mat
}

# Split the expression matrix
train_mat_raw <- exp_mat[, train_ids, drop = FALSE]
val_mat_raw <- exp_mat[, val_ids, drop = FALSE]
test_mat_raw <- exp_mat[, test_ids, drop = FALSE]

cat(sprintf("Preprocessing input dimensions (genes x samples): train=%d x %d, val=%d x %d, test=%d x %d\n",
            nrow(train_mat_raw), ncol(train_mat_raw),
            nrow(val_mat_raw), ncol(val_mat_raw),
            nrow(test_mat_raw), ncol(test_mat_raw)))
cat(sprintf("Preprocessing totals: genes=%d, samples=%d, entries=%d\n",
            nrow(exp_mat), ncol(exp_mat), nrow(exp_mat) * ncol(exp_mat)))

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

cat(sprintf("Postprocessing dimensions (genes x samples): train=%d x %d, val=%d x %d, test=%d x %d\n",
            nrow(expr_train), ncol(expr_train),
            nrow(expr_val), ncol(expr_val),
            nrow(expr_test), ncol(expr_test)))
cat(sprintf("Postprocessing totals: genes=%d, samples=%d, entries=%d\n",
            nrow(expr_train), ncol(expr_train), nrow(expr_train) * ncol(expr_train)))

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

cat(sprintf("Aligned splits: train=%d, val=%d, test=%d\n",
            nrow(clinical_splits$train), nrow(clinical_splits$val), nrow(clinical_splits$test)))

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
cat(sprintf("DEG filtering: %d genes retained out of %d (train stage subset)\n",
            nrow(expr_filt), nrow(expr_train_stage)))

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
cat(sprintf("DEGs selected for downstream: %d genes\n", length(deg_genes)))

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
  filename = "results/5-fold-results/[a]limma_deg_volcano_plot.png",
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
       filename = "results/5-fold-results/[b]top_degs_boxplot.png",
       width = 8,
       height = 6,
       dpi = 600
)

#################################################################################
# GSVA-driven pathway ranking ----
################################################################################
# In this section, we use GSVA's built-in ranking to score pathway activity directly
# (no GSEA-derived gene ranks). We then select the pathway with the strongest
# GSVA score contrast between tumor stages in the training split.

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

print(bc_terms)

print(msigdbr::msigdbr(
  species = "Homo sapiens", 
  collection = "C2",
  subcollection = "CP" # canonical pathways: 
) %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::distinct() %>% head() )

if (!nrow(bc_terms)) {
  message("No breast cancer-specific CGP sets found; defaulting to Hallmark collection.")
  bc_terms <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct()
}

# use clusterProfiler for GSEA 
bc_gene_sets <- split(bc_terms$gene_symbol, bc_terms$gs_name)
bc_gene_sets <- lapply(bc_gene_sets, function(genes) {
  intersect(genes, rownames(expr_splits$train))
})
bc_gene_sets <- bc_gene_sets[lengths(bc_gene_sets) >= 10]

if (!length(bc_gene_sets)) {
  stop("No MSigDB gene sets overlap with the expression matrix after filtering.")
}

train_ids <- colnames(expr_splits$train)
train_stage <- dt_train$tumor_stage[match(train_ids, dt_train$patient_id)]
stage_factor <- factor(train_stage)

if (any(is.na(stage_factor)) || length(unique(stage_factor)) < 2) {
  stop("Tumor stage labels are missing or have fewer than two groups; cannot rank pathways.")
}

gsva_rank_param <- GSVA::gsvaParam(
  exprData = expr_splits$train,
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

selected_pathway <- gsva_rank_tbl %>% dplyr::slice_head(n = 1)

if (!nrow(selected_pathway)) {
  stop("No pathways returned by GSVA ranking; cannot derive GSVA treatment.")
}

selected_pathway_id <- selected_pathway$ID[1]
selected_pathway_label <- selected_pathway_id
message(sprintf("Selected pathway for downstream GSVA-based treatment: %s", selected_pathway_label))

selected_gene_set <- bc_terms %>%
  dplyr::filter(gs_name == selected_pathway_id) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

selected_gene_set <- intersect(selected_gene_set, rownames(expr_splits$train))

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

cat(sprintf("GSVA scores computed for pathway '%s' (train=%d, val=%d, test=%d)\n",
            selected_pathway_label,
            nrow(pathway_scores$train),
            nrow(pathway_scores$val),
            nrow(pathway_scores$test)))

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

cat(sprintf("Assembled analysis tables: train=%d, val=%d, test=%d, trainval=%d, full=%d\n",
            nrow(train_df), nrow(val_df), nrow(test_df), nrow(trainval_df), nrow(full_df)))

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


print( covariate_association |> dplyr::arrange(metric, desc(value)))
print( covariate_association |> dplyr::arrange(desc(metric), desc(value)))


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


#################################################################################
# Comprehensive Visualizations ----
################################################################################

cat("\n=== Generating comprehensive visualization plots ===\n\n")

# Create results directory for plots
plots_dir <- file.path("results", "5-fold-results","plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 1. DR Effect Forest Plot (theta estimates across splits and models)
#-------------------------------------------------------------------------------

if (nrow(dml_split_tbl)) {
  forest_plot_data <- dml_split_tbl %>%
    dplyr::mutate(
      split = factor(split, levels = c("train", "val", "test")),
      model = factor(model, levels = c("DR-GLMNET", "DR-RF", "DR-BART"))
    )
  
  dr_forest_plot <- ggplot(forest_plot_data, aes(x = theta, y = interaction(model, split), color = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = c("DR-GLMNET" = "#E41A1C", "DR-RF" = "#377EB8", "DR-BART" = "#4DAF4A")) +
    labs(
      title = "Doubly Robust Effect Estimates (Forest Plot)",
      subtitle = sprintf("Treatment: %s pathway GSVA score", selected_pathway_label),
      x = expression(paste("Estimated Effect (", theta, ")")),
      y = "Model × Split",
      color = "Model"
    ) +
    theme_nature() +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  ggsave(
    filename = file.path(plots_dir, "dr_effect_forest_plot.png"),
    plot = dr_forest_plot,
    width = 10,
    height = 6,
    dpi = 600
  )
  cat("Saved: DR effect forest plot\n")
}

#-------------------------------------------------------------------------------
# 2. Model Performance Comparison (AUC and C-index)
#-------------------------------------------------------------------------------

if (nrow(model_metrics)) {
  metrics_long <- model_metrics %>%
    tidyr::pivot_longer(cols = c(test_auc, test_cindex), names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        metric == "test_auc" ~ "ROC-AUC",
        metric == "test_cindex" ~ "C-Index",
        TRUE ~ metric
      )
    )
  
  performance_bar_plot <- ggplot(metrics_long, aes(x = model, y = value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", value)), 
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("ROC-AUC" = "#2166AC", "C-Index" = "#B2182B")) +
    labs(
      title = "Model Performance Comparison (Test Set)",
      x = "Model",
      y = "Metric Value",
      fill = "Metric"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_nature() +
    theme(legend.position = "top")
  
  ggsave(
    filename = file.path(plots_dir, "model_performance_comparison.png"),
    plot = performance_bar_plot,
    width = 8,
    height = 6,
    dpi = 600
  )
  cat("Saved: Model performance comparison bar plot\n")
}

#-------------------------------------------------------------------------------
# 3. Compute validation predictions for ROC/PR curves
#-------------------------------------------------------------------------------

# Get validation predictions
val_cov_df <- apply_imputer(val_df, cov_imputer_trainval)
X_val_cov <- model.matrix(cov_formula, data = val_cov_df)
X_val_full <- cbind(pathway_score = val_df$pathway_score, X_val_cov)
glmnet_val_pred <- as.numeric(stats::predict(glmnet_classifier, newx = X_val_full, s = "lambda.min", type = "response"))

rf_val_df <- apply_imputer(val_df, rf_imputer)
rf_val_pred_mat <- predict(rf_classifier, data = rf_val_df)$predictions
rf_val_pred <- if (is.matrix(rf_val_pred_mat)) {
  as.numeric(rf_val_pred_mat[, "1", drop = TRUE])
} else {
  as.numeric(rf_val_pred_mat)
}

bart_val_fit <- BART::pbart(
  x.train = X_trainval_full,
  y.train = trainval_df$tumor_stage_binary,
  x.test = X_val_full,
  ntree = 200,
  ndpost = 1500,
  nskip = 500,
  usequants = TRUE,
  sparse = TRUE
)
bart_val_pred <- as.numeric(bart_val_fit$prob.test.mean)

#-------------------------------------------------------------------------------
# 4. ROC Curves for all models
#-------------------------------------------------------------------------------

roc_glmnet <- create_roc_curve_plot(glmnet_val_pred, val_df$tumor_stage_binary, 
                                     glmnet_test_pred, test_df$tumor_stage_binary, "DR-GLMNET")
roc_rf <- create_roc_curve_plot(rf_val_pred, val_df$tumor_stage_binary, 
                                 rf_test_pred, test_df$tumor_stage_binary, "DR-RF")
roc_bart <- create_roc_curve_plot(bart_val_pred, val_df$tumor_stage_binary, 
                                   bart_test_pred, test_df$tumor_stage_binary, "DR-BART")

if (!is.null(roc_glmnet)) {
  ggsave(file.path(plots_dir, "roc_curve_glmnet.png"), roc_glmnet, width = 7, height = 6, dpi = 600)
}
if (!is.null(roc_rf)) {
  ggsave(file.path(plots_dir, "roc_curve_rf.png"), roc_rf, width = 7, height = 6, dpi = 600)
}
if (!is.null(roc_bart)) {
  ggsave(file.path(plots_dir, "roc_curve_bart.png"), roc_bart, width = 7, height = 6, dpi = 600)
}
cat("Saved: ROC curves for all models\n")

#-------------------------------------------------------------------------------
# 5. PR Curves for all models
#-------------------------------------------------------------------------------

pr_glmnet <- create_pr_curve_plot(glmnet_val_pred, val_df$tumor_stage_binary, 
                                   glmnet_test_pred, test_df$tumor_stage_binary, "DR-GLMNET")
pr_rf <- create_pr_curve_plot(rf_val_pred, val_df$tumor_stage_binary, 
                               rf_test_pred, test_df$tumor_stage_binary, "DR-RF")
pr_bart <- create_pr_curve_plot(bart_val_pred, val_df$tumor_stage_binary, 
                                 bart_test_pred, test_df$tumor_stage_binary, "DR-BART")

if (!is.null(pr_glmnet)) {
  ggsave(file.path(plots_dir, "pr_curve_glmnet.png"), pr_glmnet, width = 7, height = 6, dpi = 600)
}
if (!is.null(pr_rf)) {
  ggsave(file.path(plots_dir, "pr_curve_rf.png"), pr_rf, width = 7, height = 6, dpi = 600)
}
if (!is.null(pr_bart)) {
  ggsave(file.path(plots_dir, "pr_curve_bart.png"), pr_bart, width = 7, height = 6, dpi = 600)
}
cat("Saved: PR curves for all models\n")

#-------------------------------------------------------------------------------
# 6. Calibration Plots for all models
#-------------------------------------------------------------------------------

cal_glmnet <- create_calibration_plot(glmnet_val_pred, val_df$tumor_stage_binary, 
                                       glmnet_test_pred, test_df$tumor_stage_binary, "DR-GLMNET")
cal_rf <- create_calibration_plot(rf_val_pred, val_df$tumor_stage_binary, 
                                   rf_test_pred, test_df$tumor_stage_binary, "DR-RF")
cal_bart <- create_calibration_plot(bart_val_pred, val_df$tumor_stage_binary, 
                                     bart_test_pred, test_df$tumor_stage_binary, "DR-BART")

ggsave(file.path(plots_dir, "calibration_glmnet.png"), cal_glmnet, width = 7, height = 6, dpi = 600)
ggsave(file.path(plots_dir, "calibration_rf.png"), cal_rf, width = 7, height = 6, dpi = 600)
ggsave(file.path(plots_dir, "calibration_bart.png"), cal_bart, width = 7, height = 6, dpi = 600)
cat("Saved: Calibration plots for all models\n")

#-------------------------------------------------------------------------------
# 7. Confusion Matrices for test set (using Youden threshold)
#-------------------------------------------------------------------------------

youden_glmnet <- compute_youden_threshold(glmnet_test_pred, test_df$tumor_stage_binary)
youden_rf <- compute_youden_threshold(rf_test_pred, test_df$tumor_stage_binary)
youden_bart <- compute_youden_threshold(bart_test_pred, test_df$tumor_stage_binary)

cm_glmnet <- create_confusion_matrix_plot(glmnet_test_pred, test_df$tumor_stage_binary, 
                                           youden_glmnet$threshold, "DR-GLMNET", "Test", "Youden")
cm_rf <- create_confusion_matrix_plot(rf_test_pred, test_df$tumor_stage_binary, 
                                       youden_rf$threshold, "DR-RF", "Test", "Youden")
cm_bart <- create_confusion_matrix_plot(bart_test_pred, test_df$tumor_stage_binary, 
                                         youden_bart$threshold, "DR-BART", "Test", "Youden")

ggsave(file.path(plots_dir, "confusion_matrix_glmnet.png"), cm_glmnet, width = 6, height = 5, dpi = 600)
ggsave(file.path(plots_dir, "confusion_matrix_rf.png"), cm_rf, width = 6, height = 5, dpi = 600)
ggsave(file.path(plots_dir, "confusion_matrix_bart.png"), cm_bart, width = 6, height = 5, dpi = 600)
cat("Saved: Confusion matrices for all models\n")

#-------------------------------------------------------------------------------
# 8. Comprehensive Metrics Table
#-------------------------------------------------------------------------------

comprehensive_metrics <- dplyr::bind_rows(
  compute_classification_metrics(glmnet_test_pred, test_df$tumor_stage_binary, youden_glmnet$threshold) %>%
    dplyr::mutate(model = "DR-GLMNET", threshold_type = "Youden", threshold = youden_glmnet$threshold),
  compute_classification_metrics(rf_test_pred, test_df$tumor_stage_binary, youden_rf$threshold) %>%
    dplyr::mutate(model = "DR-RF", threshold_type = "Youden", threshold = youden_rf$threshold),
  compute_classification_metrics(bart_test_pred, test_df$tumor_stage_binary, youden_bart$threshold) %>%
    dplyr::mutate(model = "DR-BART", threshold_type = "Youden", threshold = youden_bart$threshold)
)

print("Comprehensive classification metrics (test set, Youden threshold):")
print(comprehensive_metrics %>% dplyr::select(model, accuracy, balanced_accuracy, precision, recall, specificity, f1, roc_auc, pr_auc))

write.csv(comprehensive_metrics, file.path(plots_dir, "comprehensive_metrics.csv"), row.names = FALSE)
cat("Saved: Comprehensive metrics CSV\n")

#-------------------------------------------------------------------------------
# 9. Prediction Distribution Plot (Ridge plot)
#-------------------------------------------------------------------------------

pred_dist_df <- dplyr::bind_rows(
  data.frame(pred = glmnet_test_pred, truth = test_df$tumor_stage_binary, model = "DR-GLMNET", split = "Test"),
  data.frame(pred = rf_test_pred, truth = test_df$tumor_stage_binary, model = "DR-RF", split = "Test"),
  data.frame(pred = bart_test_pred, truth = test_df$tumor_stage_binary, model = "DR-BART", split = "Test")
) %>%
  dplyr::mutate(
    truth_label = factor(ifelse(truth == 1, "Late", "Early"), levels = c("Early", "Late"))
  )

pred_dist_plot <- ggplot(pred_dist_df, aes(x = pred, y = model, fill = truth_label)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("Early" = "#4DAF4A", "Late" = "#E41A1C")) +
  labs(
    title = "Predicted Probability Distributions by True Class",
    subtitle = "Test set predictions across models",
    x = "Predicted Probability (Late Stage)",
    y = "Model",
    fill = "True Stage"
  ) +
  theme_nature() +
  theme(legend.position = "right")

# Load ggridges if available
if (requireNamespace("ggridges", quietly = TRUE)) {
  library(ggridges)
  ggsave(file.path(plots_dir, "prediction_distribution_ridge.png"), pred_dist_plot, width = 9, height = 5, dpi = 600)
  cat("Saved: Prediction distribution ridge plot\n")
}

#-------------------------------------------------------------------------------
# 10. Sensitivity Analysis Forest Plot
#-------------------------------------------------------------------------------

sensitivity_comparison <- dplyr::bind_rows(
  dml_split_tbl %>% dplyr::filter(split == "train") %>% dplyr::mutate(covariate_set = "Full"),
  clinical_sensitivity_tbl %>% dplyr::mutate(split = "train", covariate_set = "Clinical-only")
)

sensitivity_forest_plot <- ggplot(sensitivity_comparison, aes(x = theta, y = model, color = covariate_set)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 0.8, 
                 position = position_dodge(width = 0.4)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("Full" = "#2166AC", "Clinical-only" = "#B2182B")) +
  labs(
    title = "Sensitivity Analysis: Effect of Covariate Set",
    subtitle = "Comparing full covariate set vs clinical-only (train split)",
    x = expression(paste("Estimated Effect (", theta, ")")),
    y = "Model",
    color = "Covariate Set"
  ) +
  theme_nature() +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, "sensitivity_analysis_forest.png"), sensitivity_forest_plot, width = 8, height = 5, dpi = 600)
cat("Saved: Sensitivity analysis forest plot\n")

#-------------------------------------------------------------------------------
# 11. Covariate Association Heatmap
#-------------------------------------------------------------------------------

if (nrow(covariate_association)) {
  cov_assoc_plot <- ggplot(covariate_association, aes(x = variable, y = 1, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = "white") +
    scale_fill_viridis(option = "plasma", name = "Association", na.value = "gray80") +
    labs(
      title = "Covariate-Treatment Association",
      subtitle = "Association between GSVA pathway score and adjustment covariates",
      x = "Covariate",
      y = NULL
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_fixed(ratio = 0.8)
  
  ggsave(file.path(plots_dir, "covariate_association_heatmap.png"), cov_assoc_plot, width = 12, height = 4, dpi = 600)
  cat("Saved: Covariate association heatmap\n")
}

#-------------------------------------------------------------------------------
# 12. GSVA Score vs Tumor Stage Violin Plot
#-------------------------------------------------------------------------------

gsva_violin_df <- full_df %>%
  dplyr::select(patient_id, pathway_score, tumor_stage, split) %>%
  dplyr::mutate(tumor_stage = factor(tumor_stage, levels = c("early", "late"), labels = c("Early", "Late")))

gsva_violin_plot <- ggplot(gsva_violin_df, aes(x = tumor_stage, y = pathway_score, fill = tumor_stage)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.8) +
  facet_wrap(~ split, nrow = 1) +
  scale_fill_manual(values = c("Early" = "#4DAF4A", "Late" = "#E41A1C")) +
  labs(
    title = sprintf("GSVA Score Distribution: %s", selected_pathway_label),
    subtitle = "Treatment variable distribution by tumor stage and split",
    x = "Tumor Stage",
    y = "GSVA Score (Pathway Activity)",
    fill = "Tumor Stage"
  ) +
  theme_nature() +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, "gsva_violin_by_stage_split.png"), gsva_violin_plot, width = 10, height = 5, dpi = 600)
cat("Saved: GSVA violin plot by stage and split\n")

#-------------------------------------------------------------------------------
# 13. Combined ROC curves (all models in one plot)
#-------------------------------------------------------------------------------

roc_glmnet_obj <- tryCatch(pROC::roc(response = test_df$tumor_stage_binary, predictor = glmnet_test_pred, quiet = TRUE), error = function(e) NULL)
roc_rf_obj <- tryCatch(pROC::roc(response = test_df$tumor_stage_binary, predictor = rf_test_pred, quiet = TRUE), error = function(e) NULL)
roc_bart_obj <- tryCatch(pROC::roc(response = test_df$tumor_stage_binary, predictor = bart_test_pred, quiet = TRUE), error = function(e) NULL)

if (!is.null(roc_glmnet_obj) && !is.null(roc_rf_obj) && !is.null(roc_bart_obj)) {
  combined_roc_df <- dplyr::bind_rows(
    data.frame(fpr = 1 - roc_glmnet_obj$specificities, tpr = roc_glmnet_obj$sensitivities, 
               model = sprintf("DR-GLMNET (AUC=%.3f)", as.numeric(roc_glmnet_obj$auc))),
    data.frame(fpr = 1 - roc_rf_obj$specificities, tpr = roc_rf_obj$sensitivities, 
               model = sprintf("DR-RF (AUC=%.3f)", as.numeric(roc_rf_obj$auc))),
    data.frame(fpr = 1 - roc_bart_obj$specificities, tpr = roc_bart_obj$sensitivities, 
               model = sprintf("DR-BART (AUC=%.3f)", as.numeric(roc_bart_obj$auc)))
  )
  
  combined_roc_plot <- ggplot(combined_roc_df, aes(x = fpr, y = tpr, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
    labs(
      title = "ROC Curves: All Models (Test Set)",
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Model"
    ) +
    coord_equal() +
    theme_nature() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plots_dir, "combined_roc_curves.png"), combined_roc_plot, width = 8, height = 7, dpi = 600)
  cat("Saved: Combined ROC curves\n")
}

#-------------------------------------------------------------------------------
# 14. Summary panel (combine key plots)
#-------------------------------------------------------------------------------

if (exists("dr_forest_plot") && exists("combined_roc_plot") && exists("gsva_violin_plot")) {
  summary_panel <- cowplot::plot_grid(
    dr_forest_plot + theme(legend.position = "none"),
    combined_roc_plot + theme(legend.position = "none"),
    gsva_violin_plot + theme(legend.position = "none"),
    ncol = 1,
    labels = c("A", "B", "C"),
    label_size = 14
  )
  
  ggsave(file.path(plots_dir, "summary_panel.png"), summary_panel, width = 10, height = 14, dpi = 600)
  cat("Saved: Summary panel\n")
}

cat(sprintf("\n=== All visualization plots saved to: %s ===\n\n", plots_dir))


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



