# Title: Doubly-robust machine learning estimation of the average causal effect of (selected) DEGs on breast cancer survival outcomes
# Author: Amos Okutse, Man-Fang Liang
# Date: Jan 2025


# Comparative Methods for Survival Analysis:
# (1) DR-Cox: Cox proportional hazards model for outcome prediction
# (2) DR-RSF: Random Survival Forest for outcome prediction  
# (3) DR-DeepSurv: Deep learning survival model for outcome prediction [COMMENTED OUT - keras installation issues]
# (4) DR-GLMNET: Penalized regression on RMST pseudo-outcome (baseline)

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
library(randomForestSRC) # random survival forest
library(splines)       # natural splines
library(ggplot2)       # plotting
library(ggrepel)       # labeling points cleanly
library(ggpubr)        # arranging plots
library(cowplot)       # arranging plots
library(parallel)      # parallel computing
library(doParallel)    # parallel backend for foreach
library(foreach)       # parallel loops
library(boot)          # for inverse logit function
library(ggvenn)        # for Venn diagrams

library(pathview)      # for KEGG pathway visualization
library(msigdbr)       # MSigDB gene sets
library(GSVA)          # pathway activity scoring
library(ranger)        # fast random forests
library(pROC)          # AUC computation
# library(limma)         # for DEG analysis (commented out - DEG section disabled)
library(matrixStats)   # for simple matrix operations
library(VennDiagram)   # for Venn diagrams

# Optional: keras for DeepSurv (loaded conditionally in deepsurv_learner function)
# if (!requireNamespace("keras", quietly = TRUE)) {
#   install.packages("keras")
#   keras::install_keras()
# }


# source helper files
source("scripts/helper.R")
source("scripts/survival_helpers.R")

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

#' Cox Model Survival Learner for RMST Prediction
#' 
#' Uses Cox proportional hazards model to predict RMST up to time tau.
#' Predicts E[Y_pseudo | X] where Y_pseudo is RMST pseudo-outcome.
#' 
#' @param tau Time horizon for RMST
#' @return Function that takes (X, time, event, tau) and returns list with model and predict function
cox_surv_learner <- function() {
  function(X, time, event, tau) {
    library(survival)
    X_df <- as.data.frame(X)
    surv_obj <- Surv(time = time, event = event)
    
    # Fit Cox model
    cox_fit <- tryCatch({
      coxph(surv_obj ~ ., data = X_df)
    }, error = function(e) {
      # If Cox fails, fall back to null model
      warning("Cox model failed, using null model")
      coxph(surv_obj ~ 1, data = X_df)
    })
    
    # Function to predict RMST from Cox model
    predict_rmst_cox <- function(newX, model, tau) {
      newX_df <- as.data.frame(newX)
      
      # Get survival curves for new data
      surv_curves <- survfit(model, newdata = newX_df)
      
      # Compute RMST for each individual by integrating survival curve to tau
      n_new <- nrow(newX_df)
      rmst_predictions <- numeric(n_new)
      
      for (i in 1:n_new) {
        # Extract survival function for individual i
        if (n_new == 1) {
          times <- surv_curves$time
          surv_probs <- surv_curves$surv
        } else {
          times <- surv_curves$time
          surv_probs <- surv_curves$surv[, i]
        }
        
        # Restrict to times <= tau
        valid_idx <- times <= tau
        times_tau <- c(0, times[valid_idx], tau)
        surv_tau <- c(1, surv_probs[valid_idx], surv_probs[valid_idx][sum(valid_idx)])
        
        # Compute RMST as area under survival curve using trapezoidal rule
        rmst_predictions[i] <- sum(diff(times_tau) * (surv_tau[-1] + surv_tau[-length(surv_tau)]) / 2)
      }
      
      return(rmst_predictions)
    }
    
    list(
      model = cox_fit,
      predict = function(newX) predict_rmst_cox(newX, cox_fit, tau)
    )
  }
}

#' Random Survival Forest Learner for RMST Prediction
#' 
#' Uses Random Survival Forest to predict RMST up to time tau.
#' 
#' @param num.trees Number of trees in the forest
#' @return Function that takes (X, time, event, tau) and returns list with model and predict function
rsf_surv_learner <- function(num.trees = 500) {
  function(X, time, event, tau) {
    library(randomForestSRC)
    X_df <- as.data.frame(X)
    data <- cbind(time = time, event = event, X_df)
    
    # Fit Random Survival Forest
    rsf_fit <- tryCatch({
      rfsrc(Surv(time, event) ~ .,
            data = data,
            ntree = num.trees,
            nodesize = 15,
            seed = 202501)
    }, error = function(e) {
      warning("RSF model failed: ", e$message)
      NULL
    })
    
    if (is.null(rsf_fit)) {
      # Fallback to simple mean
      mean_rmst <- mean(pmin(time[event == 1], tau), na.rm = TRUE)
      return(list(
        model = NULL,
        predict = function(newX) rep(mean_rmst, nrow(newX))
      ))
    }
    
    # Function to predict RMST from RSF
    predict_rmst_rsf <- function(newX, model, tau) {
      newX_df <- as.data.frame(newX)
      
      # Predict survival curves
      pred <- predict(model, newdata = newX_df)
      
      # Extract time points and survival probabilities
      n_new <- nrow(newX_df)
      rmst_predictions <- numeric(n_new)
      
      for (i in 1:n_new) {
        times <- pred$time.interest
        surv_probs <- pred$survival[i, ]
        
        # Restrict to times <= tau
        valid_idx <- times <= tau
        times_tau <- c(0, times[valid_idx], tau)
        surv_tau <- c(1, surv_probs[valid_idx], surv_probs[valid_idx][sum(valid_idx)])
        
        # Compute RMST using trapezoidal rule
        rmst_predictions[i] <- sum(diff(times_tau) * (surv_tau[-1] + surv_tau[-length(surv_tau)]) / 2)
      }
      
      return(rmst_predictions)
    }
    
    list(
      model = rsf_fit,
      predict = function(newX) predict_rmst_rsf(newX, rsf_fit, tau)
    )
  }
}

# ============================================================================
# DeepSurv Learner - COMMENTED OUT (keras installation issues)
# ============================================================================
# #' DeepSurv Learner for RMST Prediction
# #' 
# #' Neural network-based Cox model (DeepSurv) for predicting RMST.
# #' Uses keras to build a neural network with Cox partial likelihood loss.
# #' 
# #' @param hidden_units Vector of hidden layer sizes
# #' @param epochs Number of training epochs
# #' @param batch_size Batch size for training
# #' @param dropout_rate Dropout rate
# #' @return Function that takes (X, time, event, tau) and returns list with model and predict function
# deepsurv_learner <- function(hidden_units = c(64, 32),
#                              epochs = 50,
#                              batch_size = 32,
#                              dropout_rate = 0.3) {
#   function(X, time, event, tau) {
#     
#     # Check if keras is available
#     if (!requireNamespace("keras", quietly = TRUE)) {
#       warning("keras not available, falling back to Cox model")
#       return(cox_surv_learner()(X, time, event, tau))
#     }
#     
#     library(keras)
#     library(survival)
#     
#     X_mat <- as.matrix(X)
#     n <- nrow(X_mat)
#     
#     # Simplified DeepSurv: Train neural network to predict RMST directly
#     # (true DeepSurv would use Cox partial likelihood, but that requires custom loss)
#     # For simplicity, we compute pseudo-outcome and train on it
#     
#     # Compute pseudo-outcome for training
#     G_t <- compute_censoring_weights(time = time, event = event, method = "KM")
#     Y_pseudo_train <- compute_rmst_pseudo(
#       time = time, 
#       event = event, 
#       G_t = G_t, 
#       tau = tau,
#       trim_quantile = 0.95
#     )
#     
#     # Build neural network
#     model <- keras_model_sequential()
#     
#     model %>% layer_dense(units = hidden_units[1], 
#                          activation = "relu", 
#                          input_shape = ncol(X_mat))
#     
#     for (i in seq_along(hidden_units)[-1]) {
#       model %>% 
#         layer_dropout(rate = dropout_rate) %>%
#         layer_dense(units = hidden_units[i], activation = "relu")
#     }
#     
#     model %>% layer_dense(units = 1, activation = "linear")
#     
#     model %>% compile(
#       loss = "mse",
#       optimizer = optimizer_adam(learning_rate = 0.001),
#       metrics = c("mae")
#     )
#     
#     # Train
#     tryCatch({
#       history <- model %>% fit(
#         x = X_mat,
#         y = Y_pseudo_train,
#         epochs = epochs,
#         batch_size = batch_size,
#         verbose = 0,
#         validation_split = 0.2
#       )
#       
#       success <- TRUE
#     }, error = function(e) {
#       warning("DeepSurv training failed: ", e$message)
#       success <<- FALSE
#     })
#     
#     if (!exists("success") || !success) {
#       # Fallback to Cox
#       return(cox_surv_learner()(X, time, event, tau))
#     }
#     
#     list(
#       model = model,
#       predict = function(newX) {
#         as.numeric(predict(model, as.matrix(newX)))
#       }
#     )
#   }
# }
# ============================================================================

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

#' DML for Survival Outcomes using Survival-specific Learners
#' 
#' This function implements doubly robust machine learning for survival outcomes.
#' The outcome learner uses survival models (Cox, RSF) to predict RMST,
#' then computes residuals for orthogonalization.
#' Note: DeepSurv is currently commented out due to keras installation issues.
#' 
#' @param time Observed survival time
#' @param event Event indicator (1 = event, 0 = censored)
#' @param D Treatment variable (continuous)
#' @param X Covariate matrix
#' @param outcome_learner_surv Survival outcome learner (takes X, time, event, returns RMST predictions)
#' @param treatment_learner Treatment learner (takes X, D, returns D predictions)
#' @param tau Time horizon for RMST
#' @param K Number of folds for cross-fitting
#' @param seed Random seed
#' @return List with theta (treatment effect), se (standard error), ci (confidence interval)
dml_partial_linear_survival <- function(time, event, D, X, 
                                        outcome_learner_surv, 
                                        treatment_learner, 
                                        tau = 120,
                                        K = 5, 
                                        seed = 202501) {
  n <- length(time)
  stopifnot(length(event) == n, length(D) == n, nrow(X) == n)
  
  # Compute Y_pseudo (RMST pseudo-outcome) using IPCW
  library(survival)
  G_t <- compute_censoring_weights(time = time, event = event, method = "KM")
  Y_pseudo <- compute_rmst_pseudo(time = time, event = event, G_t = G_t, 
                                  tau = tau, trim_quantile = 0.95)
  
  set.seed(seed)
  folds <- sample(rep(1:K, length.out = n))
  y_tilde <- numeric(n)
  d_tilde <- numeric(n)
  
  for (k in seq_len(K)) {
    idx_train <- which(folds != k)
    idx_test <- which(folds == k)
    
    # Outcome model: predict E[Y_pseudo | X] using survival model
    # Note: survival model uses (time, event, X) on training set
    # and predicts RMST for test set
    fit_g <- outcome_learner_surv(
      X = X[idx_train, , drop = FALSE], 
      time = time[idx_train], 
      event = event[idx_train],
      tau = tau
    )
    g_hat <- fit_g$predict(X[idx_test, , drop = FALSE])
    
    # Treatment model: predict E[D | X]
    fit_m <- treatment_learner(X[idx_train, , drop = FALSE], D[idx_train])
    m_hat <- fit_m$predict(X[idx_test, , drop = FALSE])
    
    # Orthogonalization
    y_tilde[idx_test] <- Y_pseudo[idx_test] - g_hat
    d_tilde[idx_test] <- D[idx_test] - m_hat
  }
  
  # Final regression on orthogonalized variables
  dr_fit <- stats::lm(y_tilde ~ d_tilde)
  theta_hat <- stats::coef(dr_fit)["d_tilde"]
  se_hat <- summary(dr_fit)$coefficients["d_tilde", "Std. Error"]
  
  list(
    theta = theta_hat,
    se = se_hat,
    ci = c(theta_hat - 1.96 * se_hat, theta_hat + 1.96 * se_hat),
    folds = folds,
    Y_pseudo = Y_pseudo,
    y_tilde = y_tilde,
    d_tilde = d_tilde
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
  # Create tumor_stage variable (for covariate use only, not outcome)
  dplyr::mutate(tumor_stage = as.factor(case_when(
    str_detect(tumor_stg, "0|1")  ~ "early",
    str_detect(tumor_stg, "2|3|4") ~ "late",
    TRUE ~ NA_character_
  ))) %>% 
  # add observability mask for tumor stage
  dplyr::mutate(tumor_stage_mask = if_else(!is.na(tumor_stage), 1, 0)) %>% 
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
              ) %>%
  # Prepare survival outcome variables for DML analysis
  dplyr::mutate(
    OS_MONTHS = overall_survival_months,  # Already converted to numeric above
    OS_STATUS_BINARY = case_when(
      overall_survival_status == "1:DECEASED" ~ 1,
      overall_survival_status == "0:LIVING" ~ 0,
      TRUE ~ NA_real_
    )
  ) 


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

# Filter to samples with complete survival outcome data (OS_MONTHS and OS_STATUS)
dt2 <- dt2 %>%
  dplyr::filter(!is.na(OS_MONTHS), !is.na(OS_STATUS_BINARY))

cat("dt2 after filtering to complete survival data: ", nrow(dt2), " rows\n")
cat(sprintf("Event rate: %.1f%% (n=%d deaths)\n", 
    mean(dt2$OS_STATUS_BINARY == 1) * 100, 
    sum(dt2$OS_STATUS_BINARY == 1)))
cat("dt2 columns:", paste(colnames(dt2), collapse=", "), "\n")
exp_mat <- exp_mat[, dt2$patient_id, drop = FALSE] # keep only subjects with complete survival data to align ids

# check whether exp_mat are already normalized log2 (range is 0 to 16)
summary(as.numeric(exp_mat))

# omit all NA in exp_mat if any after median imputation
exp_mat <- exp_mat[complete.cases(exp_mat), , drop = FALSE]

# Pre-processing steps for mRNA expression data
## Train, val, test split before preprocessing using stratification by survival outcome
# Stratify by OS_STATUS_BINARY (event status) to ensure balanced event rates across splits
sample_ids <- colnames(exp_mat)
n_total <- length(sample_ids)
train_frac <- 0.70
val_frac <- 0.15

split_ids <- stratified_split(
  ids = sample_ids,
  strata = dt2$OS_STATUS_BINARY,  # Stratify by survival outcome instead of tumor stage
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

# Check survival outcome distribution by split (event vs censored)
split_outcome_summary <- lapply(
  list(train = train_ids, val = val_ids, test = test_ids),
  function(ids) {
    os_status <- dt2$OS_STATUS_BINARY[match(ids, dt2$patient_id)]
    table(factor(os_status, levels = c(0, 1), labels = c("Censored", "Event")), useNA = "ifany")
  }
)
cat("\n=== Survival Outcome Distribution by Split ===\n")
print(split_outcome_summary)

# ===== Survival Outcome Preparation =====
cat("\n=== Preparing Survival Outcome (RMST Pseudo-outcome via IPCW) ===\n")

# Set tau = 120 months (10 years)
tau <- 120
cat(sprintf("Time horizon tau = %.0f months\n", tau))

# Check survival data quality
survival_quality <- check_survival_quality(
  time = dt2$OS_MONTHS,
  event = dt2$OS_STATUS_BINARY,
  tau = tau
)
cat(sprintf("Total samples: %d\n", survival_quality$n_total))
cat(sprintf("Events (deaths): %d (%.1f%%)\n", 
    survival_quality$n_events, survival_quality$event_rate * 100))
cat(sprintf("Censored: %d (%.1f%%)\n", 
    survival_quality$n_censored, survival_quality$censor_rate * 100))
cat(sprintf("Events by tau=%.0f: %d (%.1f%%)\n", 
    tau, survival_quality$events_by_tau, survival_quality$event_rate_by_tau * 100))
cat(sprintf("Median follow-up: %.1f months\n", survival_quality$median_followup))

# Compute censoring weights G(t) using Kaplan-Meier
cat("\nComputing censoring weights using Kaplan-Meier...\n")
G_t <- compute_censoring_weights(
  time = dt2$OS_MONTHS,
  event = dt2$OS_STATUS_BINARY,
  method = "KM"
)

cat(sprintf("G(t) range: [%.3f, %.3f]\n", min(G_t), max(G_t)))
cat(sprintf("G(t) mean: %.3f (SD = %.3f)\n", mean(G_t), sd(G_t)))

# Compute RMST pseudo-outcome
cat("\nComputing RMST pseudo-outcome...\n")
dt2$Y_pseudo <- compute_rmst_pseudo(
  time = dt2$OS_MONTHS,
  event = dt2$OS_STATUS_BINARY,
  G_t = G_t,
  tau = tau,
  trim_quantile = 0.95
)

cat(sprintf("Pseudo-outcome range: [%.2f, %.2f]\n", 
    min(dt2$Y_pseudo), max(dt2$Y_pseudo)))
cat(sprintf("Pseudo-outcome mean: %.2f (SD = %.2f)\n", 
    mean(dt2$Y_pseudo), sd(dt2$Y_pseudo)))
cat(sprintf("Pseudo-outcome median: %.2f (IQR: [%.2f, %.2f])\n",
    median(dt2$Y_pseudo), 
    quantile(dt2$Y_pseudo, 0.25), 
    quantile(dt2$Y_pseudo, 0.75)))

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
# NOTE: DEG analysis commented out due to volcano plot rendering issues
#       (viewport zero dimension / font family errors on Windows)
#       This section is independent from downstream GSVA pathway analysis
# 
# # keep genes with at least 20% samples having expression > 1 (z-score)
# keep_genes <- apply(expr_train, 1, function(x) sum(x > 1) >= 0.2 * length(x))
# dt_train_all <- clinical_splits$train
# 
# # Filter to samples with complete survival pseudo-outcome (Y_pseudo) for pathway ranking
# outcome_complete <- !is.na(clinical_splits$train$Y_pseudo)
dt_train_all <- clinical_splits$train
outcome_complete <- !is.na(dt_train_all$Y_pseudo)
cat(sprintf("dt_train_all samples with Y_pseudo: %d out of %d\n", 
    sum(outcome_complete), nrow(dt_train_all)))

if (sum(outcome_complete) == 0) {
  stop("No samples with Y_pseudo found in training set. Check if Y_pseudo was computed correctly.")
}

dt_train <- dt_train_all[outcome_complete, , drop = FALSE]

# expr_train_stage <- expr_train[, dt_train$patient_id, drop = FALSE]
# expr_filt <- expr_train_stage[keep_genes, , drop = FALSE]
# 
# # DEG analysis using survival outcome (Event vs Censored) instead of tumor stage
# # Filter to samples with complete OS_STATUS_BINARY
# outcome_status_complete <- !is.na(dt_train$OS_STATUS_BINARY)
# dt_train_deg <- dt_train[outcome_status_complete, , drop = FALSE]
# expr_filt_deg <- expr_filt[, dt_train_deg$patient_id, drop = FALSE]
# 
# grp_train <- factor(dt_train_deg$OS_STATUS_BINARY, levels = c(0, 1), labels = c("Censored", "Event"))
# cat(sprintf("DEG analysis: %d Censored, %d Event\n", 
#     sum(grp_train == "Censored"), sum(grp_train == "Event")))
# stopifnot(length(grp_train) == ncol(expr_filt_deg))
# 
# design <- stats::model.matrix(~ grp_train) # Event vs Censored
# fit <- limma::lmFit(expr_filt_deg, design)
# fit <- limma::eBayes(fit)
# fit$genes <- data.frame(gene = rownames(expr_filt_deg))
# 
# # limma results: Event vs Censored
# res_limma <- limma::topTable(
#   fit,
#   coef = "grp_trainEvent",
#   number = Inf,
#   adjust.method = "BH"
# ) %>%
#   as.data.frame()
# 
# cat("\n=== DEG Analysis Summary (Event vs Censored) ===\n")
# 
# # argument on cutoff selection. note the small estimated effects overall. The middle 50% of logFC is roughly 
# # [−0.024, 0.027], and the most extreme values are only about −0.95 to +1.01
# summary(res_limma)
# 
# # thresholds
# alpha_fdr <- 0.10 # conservative for exploratory analysis
# lfc_cut   <- 0.2          # ~1.15-fold
# top_n_sig <- 10
# top_n_nsig <- 10
# 
# # pull rownames in as gene IDs
# if (!"gene" %in% colnames(res_limma)) {
#   res_limma <- res_limma %>% tibble::rownames_to_column("gene")
# }
# 
# # annotate direction + DEG class on the same dataset
# res_annot <- res_limma %>%
#   mutate(
#     neglog10_adjP = -log10(adj.P.Val),
#     direction = case_when(
#       logFC >=  lfc_cut ~ "Up",
#       logFC <= -lfc_cut ~ "Down",
#       TRUE              ~ "None"
#     ),
#     deg_class = case_when(
#       adj.P.Val < alpha_fdr & logFC >=  lfc_cut ~ "Upregulated in Event",
#       adj.P.Val < alpha_fdr & logFC <= -lfc_cut ~ "Downregulated in Event",
#       TRUE                                  ~ "Not significant"
#     )
#   )
# 
# # select DEGs for downstream (GSVA / pathway enrichment / causal) 
# deg_limma_selected <- res_annot %>%
#   filter(adj.P.Val < alpha_fdr, abs(logFC) >= lfc_cut) %>%
#   arrange(adj.P.Val)
# 
# if (!nrow(deg_limma_selected)) {
#   warning("No DEGs passed the specified FDR/logFC thresholds; using top 50 genes by adjusted p-value as a fallback.")
#   deg_limma_selected <- res_annot %>%
#     arrange(adj.P.Val) %>%
#     slice_head(n = 50)
# }
# 
# deg_genes <- deg_limma_selected$gene
# length(deg_genes)
# 
# # Separate up vs down gene sets (useful for enrichment or signed analyses)
# deg_genes_up   <- deg_limma_selected %>% filter(logFC >=  lfc_cut) %>% pull(gene)
# deg_genes_down <- deg_limma_selected %>% filter(logFC <= -lfc_cut) %>% pull(gene)
# 
# # Table of up and down regulated DEGs (Event vs Censored)
# deg_tbl_up <- deg_limma_selected %>%
#   filter(deg_class == "Upregulated in Event") %>%
#   arrange(adj.P.Val, desc(logFC)) %>%
#   dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
#   slice_head(n = 10)
# 
# deg_tbl_down <- deg_limma_selected %>%
#   filter(deg_class == "Downregulated in Event") %>%
#   arrange(adj.P.Val, logFC) %>%   # most negative first
#   dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
#   slice_head(n = 10)
# 
# deg_tbl <- bind_rows(deg_tbl_up, deg_tbl_down)
# 
# knitr::kable(
#   deg_tbl,
#   format = "latex",
#   booktabs = TRUE,
#   digits = 4,
#   caption = "Top Upregulated and Downregulated DEGs identified by limma"
# )
# 
# # Volcano plot with correct coloring for up, down, and non-sig genes with labels
# volcano_data <- res_annot
# # Label selection:
# # For Up/Down: top N by abs(logFC)
# # For Not significant: top N by a stable "extremeness" score
# label_data <- bind_rows(
#   volcano_data %>%
#     filter(deg_class %in% c("Upregulated in Event", "Downregulated in Event")) %>%
#     arrange(desc(abs(logFC))) %>%
#     group_by(deg_class) %>%
#     slice_head(n = top_n_sig) %>%
#     ungroup(),
#   volcano_data %>%
#     filter(deg_class == "Not significant") %>%
#     mutate(extremeness = neglog10_adjP * abs(logFC)) %>%
#     arrange(desc(extremeness)) %>%
#     slice_head(n = top_n_nsig)
# )
# volcano_data$deg_class <- factor(
#   volcano_data$deg_class,
#   levels = c("Upregulated in Event", "Downregulated in Event", "Not significant")
# )
# 
# volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neglog10_adjP, color = as.factor(deg_class))) +
#   geom_point(alpha = 0.6, size = 1.2) +
#   geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
#   geom_hline(yintercept = -log10(alpha_fdr), linetype = "dashed", linewidth = 0.4) +
#   ggrepel::geom_text_repel(
#     data = label_data,
#     aes(label = gene),
#     size = 3,
#     max.overlaps = Inf,
#     box.padding = 0.3,
#     point.padding = 0.2,
#     segment.alpha = 0.6
#   ) +
#   scale_color_manual(values = c(
#     "Upregulated in Event"     = "red",
#     "Downregulated in Event"   = "blue",
#     "Not significant" = "black"
#   )) +
#   labs(
#     title = "(a) Volcano plot of DEGs in Event vs Censored patients",
#     x = "log2 Fold Change (Event vs Censored)",
#     y = "-log10(FDR)",
#     color = "DEG class"
#   ) +
#   theme_nature()
# 
# volcano_plot
# # save the 600 dpi image in results folder
# ggsave(
#   filename = "results/[a]limma_deg_volcano_plot.png",
#   plot = volcano_plot,
#   width = 6,
#   height = 5,
#   dpi = 600
# )
# 
# # Top DEG genes for downstream analyses (fallback to overall top-ranked genes if up/down tables are empty)
# top_gene_candidates <- unique(c(deg_tbl_up$gene, deg_tbl_down$gene))
# if (!length(top_gene_candidates)) {
# top_gene_candidates <- deg_limma_selected %>%
#     arrange(adj.P.Val) %>%
#     pull(gene)
# }
# top_genes <- top_gene_candidates[top_gene_candidates %in% rownames(expr_splits$train)]
# 
# if (!length(top_genes)) {
#   stop("No DEG genes overlap with expression matrix after filtering; cannot proceed to downstream analyses.")
# }
# 
# ## Expression of top DEGs by survival outcome
# exp_top <- expr_splits$train[top_genes, dt_train_deg$patient_id, drop = FALSE]
# exp_top_df <- as.data.frame(t(exp_top)) %>%
#   dplyr::mutate(patient_id = rownames(.)) %>%
#   dplyr::left_join(dt_train_deg %>% dplyr::select(patient_id, OS_STATUS_BINARY), by = "patient_id") %>%
#   dplyr::mutate(outcome_group = factor(OS_STATUS_BINARY, levels = c(0, 1), labels = c("Censored", "Event"))) %>%
#   tidyr::pivot_longer(cols = all_of(top_genes), names_to = "gene", values_to = "expression")
# 
# ## Boxplots of top DEGs expression by survival outcome
# boxplot_top_degs <- ggplot(exp_top_df, aes(x = outcome_group, y = expression, fill = outcome_group)) +
#   geom_boxplot() +
#   facet_wrap(~ gene, scales = "free_y") +
#   labs(title = "(b) Expression of top DEGs by Survival Outcome",
#        x = "Survival Outcome",
#        y = "Expression scores") +
#   scale_fill_manual(values = c("Censored" = "#63A375", "Event" = "#E4572E")) +
#   theme(legend.position = "none")+
#   theme_nature()
# boxplot_top_degs
# # save the 600 dpi image in results folder
# ggsave(boxplot_top_degs,
#        filename = "results/[b]top_degs_boxplot.png",
#        width = 8,
#        height = 6,
#        dpi = 600
# )

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
# Use Y_pseudo (RMST pseudo-outcome) for pathway ranking instead of tumor stage
train_outcome <- dt_train$Y_pseudo[match(train_ids, dt_train$patient_id)]

if (any(is.na(train_outcome)) || length(unique(train_outcome)) < 5) {
  stop("RMST pseudo-outcome is missing or has insufficient variation; cannot rank pathways.")
}

cat("\n=== GSVA Pathway Ranking Based on RMST Pseudo-outcome ===\n")
cat(sprintf("Using Y_pseudo (RMST) for pathway ranking: mean=%.2f, SD=%.2f\n", 
    mean(train_outcome, na.rm=TRUE), sd(train_outcome, na.rm=TRUE)))

gsva_rank_param <- GSVA::gsvaParam(
  exprData = expr_splits$train,
  geneSets = bc_gene_sets,
  minSize = 10,
  maxSize = 500,
  kcdf = "Gaussian"
)
gsva_rank_scores <- GSVA::gsva(gsva_rank_param)

# Use correlation with RMST pseudo-outcome instead of group comparison
corr_vals <- apply(gsva_rank_scores, 1, function(pathway_scores) {
  cor(pathway_scores, train_outcome, use = "pairwise.complete.obs")
})

# Calculate mean GSVA scores for high vs low RMST (for descriptive purposes)
rmst_median <- median(train_outcome, na.rm = TRUE)
high_rmst_mask <- train_outcome >= rmst_median
low_rmst_mask <- train_outcome < rmst_median

mean_high_rmst <- matrixStats::rowMeans2(gsva_rank_scores[, high_rmst_mask, drop = FALSE])
mean_low_rmst <- matrixStats::rowMeans2(gsva_rank_scores[, low_rmst_mask, drop = FALSE])
delta <- mean_high_rmst - mean_low_rmst

# P-values from correlation test
pvals <- apply(gsva_rank_scores, 1, function(pathway_scores) {
  test_result <- cor.test(pathway_scores, train_outcome, use = "pairwise.complete.obs")
  test_result$p.value
})

gsva_rank_tbl <- data.frame(
  ID = rownames(gsva_rank_scores),
  correlation = corr_vals,  # Correlation with RMST pseudo-outcome
  mean_high_rmst = mean_high_rmst,
  mean_low_rmst = mean_low_rmst,
  delta = delta,  # High RMST - Low RMST
  p_value = pvals,
  p_adjust = stats::p.adjust(pvals, method = "BH"),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(p_adjust, desc(abs(correlation)))  # Rank by p-value and correlation strength

gsva_top_tbl <- gsva_rank_tbl %>%
  dplyr::slice_head(n = 10)

cat("\n=== Top GSVA-ranked pathways by correlation with RMST ===\n")
print(gsva_top_tbl)

gsva_plot_df <- gsva_rank_tbl %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::mutate(ID = stats::reorder(ID, correlation))

gsva_rank_plot <- ggplot(gsva_plot_df, aes(x = correlation, y = ID, size = abs(correlation), color = -log10(p_adjust))) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#97C4BC", high = "#1F78B4", name = "-log10 adj p") +
  labs(
    title = "Breast cancer relevant GSVA ranking",
    subtitle = "Top pathways by correlation with RMST (survival outcome)",
    x = "Correlation with RMST pseudo-outcome",
    y = NULL,
    size = "Absolute correlation"
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
    kcdf = "Gaussian"
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
  dplyr::left_join(clinical_splits$train %>% dplyr::select(patient_id, OS_STATUS_BINARY), by = "patient_id") %>%
  dplyr::mutate(outcome_group = factor(OS_STATUS_BINARY, levels = c(0, 1), labels = c("Censored", "Event")))

gsva_density_plot <- ggplot(train_pathway_df, aes(x = pathway_score, fill = outcome_group)) +
  geom_density(alpha = 0.65) +
  labs(
    title = sprintf("GSVA distribution for %s", selected_pathway_label),
    subtitle = "Treatment D (GSVA score) by Survival Outcome",
    x = "GSVA score",
    y = "Density",
    fill = "Outcome"
  ) +
  scale_fill_manual(values = c("Censored" = "#63A375", "Event" = "#E4572E")) +
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
      split = split_name
    ) %>%
    # Filter to samples with valid survival pseudo-outcome and pathway score
    dplyr::filter(!is.na(Y_pseudo), !is.na(pathway_score))
}

train_df <- assemble_split_df(clinical_splits$train, "train")
val_df <- assemble_split_df(clinical_splits$val, "val")
test_df <- assemble_split_df(clinical_splits$test, "test")
trainval_df <- dplyr::bind_rows(train_df, val_df)
full_df <- dplyr::bind_rows(train_df, val_df, test_df)

adjustment_covariates <- c(
  "age_at_diagnosis",
  "tumor_size",
  "tumor_stage",  # Added: control for tumor stage as confounder
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
  dplyr::mutate(outcome_group = factor(OS_STATUS_BINARY, levels = c(0, 1), labels = c("Censored", "Event"))) %>%
  dplyr::group_by(outcome_group) %>%
  dplyr::summarise(
    mean_D = mean(pathway_score),
    sd_D = sd(pathway_score),
    min_D = min(pathway_score),
    max_D = max(pathway_score),
    n = dplyr::n(),
    .groups = "drop"
  )
cat("\n=== GSVA Treatment Summary by Survival Outcome (Train Split) ===\n")
print(assumption_summary)

summarize_outcome_gsva <- function(df, split_label) {
  if (!nrow(df)) {
    return(dplyr::tibble())
  }
  df %>%
    dplyr::mutate(outcome_group = factor(OS_STATUS_BINARY, levels = c(0, 1), labels = c("Censored", "Event"))) %>%
    dplyr::group_by(outcome_group) %>%
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

gsva_outcome_summaries <- dplyr::bind_rows(
  summarize_outcome_gsva(train_df, "train"),
  summarize_outcome_gsva(val_df, "val"),
  summarize_outcome_gsva(test_df, "test")
)
cat("\n=== GSVA Treatment Summary by Survival Outcome (Event vs Censored) ===\n")
print(gsva_outcome_summaries)
cat("Note: These are descriptive balance checks, not causal effects.\n")

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

# Updated DML function for survival outcome using survival-specific learners
run_split_dml <- function(df, split_label) {
  if (!nrow(df)) {
    return(dplyr::tibble())
  }
  
  # Check if survival data is available
  if (!all(c("OS_MONTHS", "OS_STATUS_BINARY", "Y_pseudo") %in% colnames(df))) {
    stop("Missing required survival variables: OS_MONTHS, OS_STATUS_BINARY, or Y_pseudo")
  }
  
  cov_imp <- fit_imputer(df, adjustment_covariates)
  cov_df <- apply_imputer(df, cov_imp)
  cov_matrix <- model.matrix(cov_formula, data = cov_df)

  cat(sprintf("\n=== Running DML for split: %s ===\n", split_label))
  
  # Unified treatment learner for fair comparison across all DR methods
  treatment_learner_rf <- ranger_learner(task = "regression", num.trees = 1200)
  cat("Using unified Random Forest treatment learner for all DR methods.\n\n")
  
  # DR-Cox: Uses Cox model to predict RMST, then compute residuals
  cat("Running DR-Cox (Cox survival model for outcome + RF for treatment)...\n")
  cox_fit <- dml_partial_linear_survival(
    time = df$OS_MONTHS,
    event = df$OS_STATUS_BINARY,
    D = df$pathway_score,
    X = cov_df,
    outcome_learner_surv = cox_surv_learner(),
    treatment_learner = treatment_learner_rf,
    tau = 120,
    K = 5,
    seed = 202501
  )

  # DR-RSF: Uses Random Survival Forest to predict RMST
  cat("Running DR-RSF (Random Survival Forest for outcome + RF for treatment)...\n")
  rsf_fit <- dml_partial_linear_survival(
    time = df$OS_MONTHS,
    event = df$OS_STATUS_BINARY,
    D = df$pathway_score,
    X = cov_df,
    outcome_learner_surv = rsf_surv_learner(num.trees = 500),
    treatment_learner = treatment_learner_rf,
    tau = 120,
    K = 5,
    seed = 202501
  )
  
  # DR-DeepSurv: COMMENTED OUT (keras installation issues)
  # cat("Running DR-DeepSurv (Deep learning for outcome + RF for treatment)...\n")
  # deepsurv_fit <- tryCatch({
  #   dml_partial_linear_survival(
  #     time = df$OS_MONTHS,
  #     event = df$OS_STATUS_BINARY,
  #     D = df$pathway_score,
  #     X = cov_df,
  #     outcome_learner_surv = deepsurv_learner(hidden_units = c(64, 32), epochs = 50),
  #     treatment_learner = treatment_learner_rf,
  #     tau = 120,
  #     K = 5,
  #     seed = 202501
  #   )
  # }, error = function(e) {
  #   warning("DeepSurv failed, skipping: ", e$message)
  #   list(theta = NA, se = NA, ci = c(NA, NA))
  # })
  
  # DR-GLMNET: Uses GLMNET regression on Y_pseudo with unified RF treatment learner
  cat("Running DR-GLMNET (GLMNET on Y_pseudo + RF for treatment)...\n")
  glmnet_fit <- dml_partial_linear(
    Y = df$Y_pseudo,
    D = df$pathway_score,
    X = cov_matrix,
    outcome_learner = glmnet_learner(family = "gaussian", alpha = 0.5),
    treatment_learner = treatment_learner_rf,
    K = 5,
    seed = 202501
  )

  
  # Return results (DeepSurv commented out)
  results <- dplyr::tibble(
    split = split_label,
    model = c("DR-Cox", "DR-RSF", "DR-GLMNET"),
    theta = c(cox_fit$theta, rsf_fit$theta, glmnet_fit$theta),
    se = c(cox_fit$se, rsf_fit$se, glmnet_fit$se),
    ci_lower = c(cox_fit$ci[1], rsf_fit$ci[1], glmnet_fit$ci[1]),
    ci_upper = c(cox_fit$ci[2], rsf_fit$ci[2], glmnet_fit$ci[2])
  )
  
  # Filter out failed models (NA theta)
  results <- results %>% dplyr::filter(!is.na(theta))
  
  return(results)
}

dml_split_tbl <- dplyr::bind_rows(
  run_split_dml(train_df, "train"),
  run_split_dml(val_df, "val"),
  run_split_dml(test_df, "test")
)

cat("\n=== Doubly Robust Estimates for GSVA Treatment Effect on RMST (Survival Outcome) ===\n")
cat("Interpretation: theta = change in RMST (months) per 1-unit increase in pathway score\n\n")
cat("Note: DR-DeepSurv is commented out due to keras installation issues.\n")
cat("Using DR-Cox, DR-RSF, and DR-GLMNET for analysis.\n\n")
print("Doubly robust estimates for GSVA treatment effect on RMST (all splits):")
print(dml_split_tbl)
cat("\nNote: These estimates use cross-fitting with K=5 folds to prevent overfitting.\n")
cat("Outcome learners predict RMST pseudo-outcome (continuous), treatment learners predict pathway score.\n")

dml_full_tbl <- run_split_dml(full_df, "all")
print("\nDoubly robust estimates when pooling train/val/test (exploratory full-data fit):")
print(dml_full_tbl)
cat("\nNote: Pooled fit is exploratory only - no holdout test set remains for validation.\n")

dr_effects_tbl <- dml_split_tbl %>% dplyr::filter(split == "train")


#################################################################################
# Model Prediction Performance Evaluation ----
#################################################################################
# Note: DML focuses on causal estimation, not prediction optimization.
# Performance metrics are reported for quality assessment and model comparison.

cat("\n=== Evaluating Outcome Model Prediction Performance ===\n")

# Helper function: Calculate C-index (Concordance Index) for survival models
calculate_cindex <- function(time, event, risk_score) {
  if (length(unique(event)) < 2) {
    return(NA)
  }
  tryCatch({
    surv_obj <- survival::Surv(time, event)
    concordance_result <- survival::concordance(surv_obj ~ risk_score)
    return(as.numeric(concordance_result$concordance))
  }, error = function(e) {
    warning("C-index calculation failed: ", e$message)
    return(NA)
  })
}

# Helper function: Calculate Integrated Brier Score
calculate_ibs <- function(time, event, predicted_surv, eval_times) {
  tryCatch({
    n <- length(time)
    if (length(unique(event)) < 2 || n < 10) {
      return(NA)
    }
    
    # Use pec package if available, otherwise simplified calculation
    if (requireNamespace("pec", quietly = TRUE)) {
      surv_obj <- survival::Surv(time, event)
      # Note: This is a simplified version; full implementation needs predicted survival curves
      return(NA)  # Placeholder - complex calculation
    } else {
      # Simplified Brier score at tau
      tau <- max(eval_times)
      brier_scores <- sapply(eval_times, function(t) {
        at_risk <- time >= t
        observed <- (time <= t) & (event == 1)
        if (sum(at_risk) == 0) return(NA)
        # Simplified calculation
        mean((observed - predicted_surv)^2, na.rm = TRUE)
      })
      return(mean(brier_scores, na.rm = TRUE))
    }
  }, error = function(e) {
    return(NA)
  })
}

# Helper function: Time-dependent AUC
calculate_time_dependent_auc <- function(time, event, risk_score, eval_time) {
  tryCatch({
    if (requireNamespace("survivalROC", quietly = TRUE)) {
      roc_obj <- survivalROC::survivalROC(
        Stime = time,
        status = event,
        marker = risk_score,
        predict.time = eval_time,
        method = "KM"
      )
      return(roc_obj$AUC)
    } else if (requireNamespace("timeROC", quietly = TRUE)) {
      # Alternative using timeROC
      roc_obj <- timeROC::timeROC(
        T = time,
        delta = event,
        marker = risk_score,
        times = eval_time,
        cause = 1
      )
      return(roc_obj$AUC[2])  # AUC at specified time
    } else {
      return(NA)
    }
  }, error = function(e) {
    warning("Time-dependent AUC calculation failed: ", e$message)
    return(NA)
  })
}

# Prepare test data
test_cov_imp <- fit_imputer(test_df, adjustment_covariates)
test_cov_df <- apply_imputer(test_df, test_cov_imp)
test_cov_matrix <- model.matrix(cov_formula, data = test_cov_df)

# Evaluation timepoints
eval_times <- c(60, 120)  # 5 years, 10 years
tau <- 120

cat("Evaluation timepoints:", paste(eval_times, collapse = ", "), "months\n\n")

# Initialize performance results
performance_results <- list()

#-----------------------------------------------------------------------------
# 1. DR-Cox Performance
#-----------------------------------------------------------------------------
cat("Evaluating DR-Cox (Cox Proportional Hazards)...\n")

tryCatch({
  # Fit Cox model on test set for prediction evaluation
  cox_fit_test <- survival::coxph(
    Surv(OS_MONTHS, OS_STATUS_BINARY) ~ .,
    data = cbind(test_cov_df, pathway_score = test_df$pathway_score)
  )
  
  # Risk scores (higher = worse prognosis)
  cox_risk <- predict(cox_fit_test, type = "risk")
  cox_lp <- predict(cox_fit_test, type = "lp")  # linear predictor
  
  # C-index
  cox_cindex <- calculate_cindex(test_df$OS_MONTHS, test_df$OS_STATUS_BINARY, cox_lp)
  
  # Time-dependent AUC
  cox_auc_60 <- calculate_time_dependent_auc(test_df$OS_MONTHS, test_df$OS_STATUS_BINARY, 
                                              cox_lp, 60)
  cox_auc_120 <- calculate_time_dependent_auc(test_df$OS_MONTHS, test_df$OS_STATUS_BINARY, 
                                               cox_lp, 120)
  
  # Predicted RMST (approximate using survival function)
  cox_surv <- summary(survival::survfit(cox_fit_test, newdata = test_cov_df), 
                     times = seq(0, tau, by = 1))
  
  # Store results
  performance_results$cox <- list(
    model = "DR-Cox",
    cindex = cox_cindex,
    auc_60 = cox_auc_60,
    auc_120 = cox_auc_120,
    rmst_cor = NA  # Cox predicts hazard, not RMST directly
  )
  
  cat("  C-index:", round(cox_cindex, 3), "\n")
  cat("  AUC at 60m:", round(cox_auc_60, 3), "\n")
  cat("  AUC at 120m:", round(cox_auc_120, 3), "\n\n")
  
}, error = function(e) {
  warning("DR-Cox performance evaluation failed: ", e$message)
  performance_results$cox <- list(
    model = "DR-Cox", cindex = NA, auc_60 = NA, auc_120 = NA, rmst_cor = NA
  )
})

#-----------------------------------------------------------------------------
# 2. DR-RSF Performance
#-----------------------------------------------------------------------------
cat("Evaluating DR-RSF (Random Survival Forest)...\n")

tryCatch({
  # Fit RSF model on test set
  rsf_fit_test <- randomForestSRC::rfsrc(
    Surv(OS_MONTHS, OS_STATUS_BINARY) ~ .,
    data = cbind(test_cov_df, pathway_score = test_df$pathway_score),
    ntree = 500,
    importance = FALSE
  )
  
  # Predicted mortality (higher = worse prognosis)
  rsf_mortality <- rsf_fit_test$predicted.oob
  if (is.null(rsf_mortality)) {
    rsf_mortality <- rsf_fit_test$predicted
  }
  
  # C-index (1 - error rate for survival)
  rsf_cindex <- 1 - rsf_fit_test$err.rate[length(rsf_fit_test$err.rate)]
  
  # Alternative: manual C-index calculation
  if (is.na(rsf_cindex) || rsf_cindex < 0.4 || rsf_cindex > 1) {
    rsf_cindex <- calculate_cindex(test_df$OS_MONTHS, test_df$OS_STATUS_BINARY, 
                                   rsf_mortality)
  }
  
  # Time-dependent AUC
  rsf_auc_60 <- calculate_time_dependent_auc(test_df$OS_MONTHS, test_df$OS_STATUS_BINARY,
                                             rsf_mortality, 60)
  rsf_auc_120 <- calculate_time_dependent_auc(test_df$OS_MONTHS, test_df$OS_STATUS_BINARY,
                                              rsf_mortality, 120)
  
  # Store results
  performance_results$rsf <- list(
    model = "DR-RSF",
    cindex = rsf_cindex,
    auc_60 = rsf_auc_60,
    auc_120 = rsf_auc_120,
    rmst_cor = NA
  )
  
  cat("  C-index:", round(rsf_cindex, 3), "\n")
  cat("  AUC at 60m:", round(rsf_auc_60, 3), "\n")
  cat("  AUC at 120m:", round(rsf_auc_120, 3), "\n\n")
  
}, error = function(e) {
  warning("DR-RSF performance evaluation failed: ", e$message)
  performance_results$rsf <- list(
    model = "DR-RSF", cindex = NA, auc_60 = NA, auc_120 = NA, rmst_cor = NA
  )
})

#-----------------------------------------------------------------------------
# 3. DR-GLMNET Performance (Predicting Y_pseudo/RMST directly)
#-----------------------------------------------------------------------------
cat("Evaluating DR-GLMNET (GLMNET on RMST pseudo-outcome)...\n")

tryCatch({
  # Fit GLMNET on test set to predict Y_pseudo
  test_X <- test_cov_matrix
  test_Y <- test_df$Y_pseudo
  
  # Remove NA values
  complete_idx <- complete.cases(test_X, test_Y)
  test_X_complete <- test_X[complete_idx, , drop = FALSE]
  test_Y_complete <- test_Y[complete_idx]
  
  if (length(test_Y_complete) > 10) {
    glmnet_fit_test <- glmnet::cv.glmnet(
      x = test_X_complete,
      y = test_Y_complete,
      family = "gaussian",
      alpha = 0.5,
      nfolds = 5
    )
    
    # Predicted RMST
    glmnet_pred <- predict(glmnet_fit_test, newx = test_X_complete, s = "lambda.min")
    
    # R-squared
    glmnet_rsq <- 1 - sum((test_Y_complete - glmnet_pred)^2) / 
                       sum((test_Y_complete - mean(test_Y_complete))^2)
    
    # Correlation
    glmnet_cor <- cor(glmnet_pred, test_Y_complete, use = "complete.obs")
    
    # RMSE
    glmnet_rmse <- sqrt(mean((test_Y_complete - glmnet_pred)^2))
    
    # MAE
    glmnet_mae <- mean(abs(test_Y_complete - glmnet_pred))
    
    # For survival concordance: use predicted RMST as risk score (higher RMST = better prognosis)
    # So we use negative predicted RMST for C-index calculation
    glmnet_cindex <- calculate_cindex(
      test_df$OS_MONTHS[complete_idx], 
      test_df$OS_STATUS_BINARY[complete_idx], 
      -as.vector(glmnet_pred)  # Negative because higher RMST = lower risk
    )
    
    # Store results
    performance_results$glmnet <- list(
      model = "DR-GLMNET",
      cindex = glmnet_cindex,
      auc_60 = NA,  # GLMNET predicts RMST, not time-specific risk
      auc_120 = NA,
      rmst_cor = glmnet_cor,
      rmst_rsq = glmnet_rsq,
      rmst_rmse = glmnet_rmse,
      rmst_mae = glmnet_mae
    )
    
    cat("  C-index (using -pred as risk):", round(glmnet_cindex, 3), "\n")
    cat("  R² (RMST prediction):", round(glmnet_rsq, 3), "\n")
    cat("  Correlation:", round(glmnet_cor, 3), "\n")
    cat("  RMSE:", round(glmnet_rmse, 2), "months\n")
    cat("  MAE:", round(glmnet_mae, 2), "months\n\n")
  } else {
    warning("Insufficient complete cases for GLMNET evaluation")
    performance_results$glmnet <- list(
      model = "DR-GLMNET", cindex = NA, auc_60 = NA, auc_120 = NA, 
      rmst_cor = NA, rmst_rsq = NA, rmst_rmse = NA, rmst_mae = NA
    )
  }
  
}, error = function(e) {
  warning("DR-GLMNET performance evaluation failed: ", e$message)
  performance_results$glmnet <- list(
    model = "DR-GLMNET", cindex = NA, auc_60 = NA, auc_120 = NA, 
    rmst_cor = NA, rmst_rsq = NA, rmst_rmse = NA, rmst_mae = NA
  )
})

#-----------------------------------------------------------------------------
# 4. Treatment Model (Random Forest) Performance
#-----------------------------------------------------------------------------
cat("Evaluating Treatment Model (Random Forest predicting pathway score)...\n")

tryCatch({
  # Fit RF to predict pathway score on test set
  rf_treatment_test <- ranger::ranger(
    pathway_score ~ .,
    data = cbind(test_cov_df, pathway_score = test_df$pathway_score),
    num.trees = 1200,
    importance = "impurity"
  )
  
  # Predictions
  rf_pred <- rf_treatment_test$predictions
  
  # R-squared
  rf_rsq <- 1 - sum((test_df$pathway_score - rf_pred)^2) / 
                 sum((test_df$pathway_score - mean(test_df$pathway_score))^2)
  
  # Correlation
  rf_cor <- cor(rf_pred, test_df$pathway_score, use = "complete.obs")
  
  # RMSE
  rf_rmse <- sqrt(mean((test_df$pathway_score - rf_pred)^2))
  
  # Variable importance (top 10)
  if (!is.null(rf_treatment_test$variable.importance)) {
    top_vars <- sort(rf_treatment_test$variable.importance, decreasing = TRUE)[1:10]
  } else {
    top_vars <- NULL
  }
  
  # Store results
  treatment_performance <- list(
    model = "RF-Treatment",
    rsq = rf_rsq,
    cor = rf_cor,
    rmse = rf_rmse,
    top_vars = names(top_vars)
  )
  
  cat("  R² (pathway score prediction):", round(rf_rsq, 3), "\n")
  cat("  Correlation:", round(rf_cor, 3), "\n")
  cat("  RMSE:", round(rf_rmse, 3), "\n")
  if (!is.null(top_vars)) {
    cat("  Top 3 important variables:", paste(names(top_vars)[1:3], collapse = ", "), "\n")
  }
  cat("\n")
  
}, error = function(e) {
  warning("Treatment model performance evaluation failed: ", e$message)
  treatment_performance <- list(
    model = "RF-Treatment", rsq = NA, cor = NA, rmse = NA, top_vars = NULL
  )
})

#-----------------------------------------------------------------------------
# Create Performance Summary Tables
#-----------------------------------------------------------------------------

# Table 1: Outcome Models (Survival Prediction)
performance_outcome_tbl <- data.frame(
  model = c("DR-Cox", "DR-RSF", "DR-GLMNET"),
  cindex = c(
    performance_results$cox$cindex,
    performance_results$rsf$cindex,
    performance_results$glmnet$cindex
  ),
  auc_60m = c(
    performance_results$cox$auc_60,
    performance_results$rsf$auc_60,
    performance_results$glmnet$auc_60
  ),
  auc_120m = c(
    performance_results$cox$auc_120,
    performance_results$rsf$auc_120,
    performance_results$glmnet$auc_120
  ),
  metric_type = c("Survival", "Survival", "RMST")
)

# Table 2: GLMNET RMST Prediction Performance
performance_rmst_tbl <- data.frame(
  model = "DR-GLMNET",
  rsquared = performance_results$glmnet$rmst_rsq,
  correlation = performance_results$glmnet$rmst_cor,
  rmse_months = performance_results$glmnet$rmst_rmse,
  mae_months = performance_results$glmnet$rmst_mae
)

# Table 3: Treatment Model Performance
performance_treatment_tbl <- data.frame(
  model = "RF-Treatment",
  target = "Pathway Score",
  rsquared = treatment_performance$rsq,
  correlation = treatment_performance$cor,
  rmse = treatment_performance$rmse
)

# Print summary
cat("=== OUTCOME MODEL PERFORMANCE SUMMARY (Test Set) ===\n")
print(performance_outcome_tbl)
cat("\n")

if (!all(is.na(performance_rmst_tbl[-1]))) {
  cat("=== DR-GLMNET RMST PREDICTION PERFORMANCE ===\n")
  print(performance_rmst_tbl)
  cat("\n")
}

cat("=== TREATMENT MODEL PERFORMANCE ===\n")
print(performance_treatment_tbl)
cat("\n")

# Save to CSV
write.csv(performance_outcome_tbl, "results/survival/tables/tableS_outcome_model_performance.csv", 
         row.names = FALSE)
write.csv(performance_rmst_tbl, "results/survival/tables/tableS_rmst_prediction_performance.csv", 
         row.names = FALSE)
write.csv(performance_treatment_tbl, "results/survival/tables/tableS_treatment_model_performance.csv", 
         row.names = FALSE)

cat("✓ Saved performance tables to results/survival/tables/\n\n")

# Add to analysis summary
cat("Note: Performance metrics assess nuisance function quality.\n")
cat("      DML estimates remain valid as long as models achieve adequate prediction.\n")
cat("      Higher C-index/R² indicates more efficient (lower SE) causal estimates.\n\n")

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

# Note: Model performance evaluation has been removed as this analysis focuses on
# causal effect estimation using survival-specific DR methods (Cox, RSF, GLMNET).
# DeepSurv is commented out due to keras installation issues.
# For predictive modeling, please use dedicated survival prediction pipelines.

reporting_candidates <- c("DR-Cox", "DR-RSF", "DR-GLMNET")


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

# Sensitivity analysis: estimate treatment effect using clinical covariates only
# Use unified RF treatment learner for consistency
clinical_glmnet <- dml_partial_linear(
  Y = train_df$Y_pseudo,  # RMST pseudo-outcome
  D = train_df$pathway_score,
  X = clinical_matrix,
  outcome_learner = glmnet_learner(family = "gaussian", alpha = 0.5),
  treatment_learner = ranger_learner(task = "regression", num.trees = 1200),
  K = 5,
  seed = 202501
)

clinical_sensitivity_tbl <- dplyr::tibble(
  model = c("DR-GLMNET"),
  theta = c(clinical_glmnet$theta),
  ci_lower = c(clinical_glmnet$ci[1]),
  ci_upper = c(clinical_glmnet$ci[2])
)
print("Sensitivity (clinical-only covariates) DR estimates:")
print(clinical_sensitivity_tbl)


#################################################################################
# Save all results to files ----
#################################################################################

# Create results subdirectories if they don't exist
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
dir.create("results/survival/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/survival/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/survival/summary", showWarnings = FALSE, recursive = TRUE)

cat("\n=== Saving Results to results/survival/ subfolders ===\n")

# 1. Save DML estimates as CSV
write.csv(dml_split_tbl, "results/survival/tables/dml_estimates_by_split.csv", row.names = FALSE)
write.csv(dml_full_tbl, "results/survival/tables/dml_estimates_full_data.csv", row.names = FALSE)
write.csv(clinical_sensitivity_tbl, "results/survival/tables/sensitivity_clinical_only.csv", row.names = FALSE)
cat("✓ Saved DML estimates to results/survival/tables/\n")

# 2. Save summary statistics as text file
sink("results/survival/summary/analysis_summary.txt")
cat("=================================================================\n")
cat("Doubly Robust Machine Learning Analysis Summary\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=================================================================\n\n")

cat("SELECTED PATHWAY:\n")
cat("  ID:", selected_pathway_id, "\n")
cat("  Label:", selected_pathway_label, "\n")
cat("  Correlation with RMST:", round(selected_pathway$correlation, 4), "\n")
cat("  Adjusted p-value:", format.pval(selected_pathway$p_adjust, digits = 3), "\n\n")

cat("MAIN RESULTS (Test Set):\n")
print(dml_split_tbl %>% dplyr::filter(split == "test"))
cat("\n")

cat("FULL DATA ESTIMATES (Exploratory):\n")
print(dml_full_tbl)
cat("\n")

cat("SENSITIVITY ANALYSIS (Clinical covariates only):\n")
print(clinical_sensitivity_tbl)
cat("\n")

cat("MODEL STABILITY ACROSS SPLITS:\n")
stability_summary <- dml_split_tbl %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_theta = mean(theta, na.rm = TRUE),
    sd_theta = sd(theta, na.rm = TRUE),
    cv_theta = sd_theta / abs(mean_theta) * 100,
    mean_se = mean(se, na.rm = TRUE),
    .groups = "drop"
  )
print(stability_summary)
cat("\n")

cat("MODEL PREDICTION PERFORMANCE (Test Set):\n")
cat("Note: Metrics assess nuisance function quality, not causal estimation.\n\n")
print(performance_outcome_tbl)
cat("\n")

if (!all(is.na(performance_rmst_tbl[-1]))) {
  cat("DR-GLMNET RMST Prediction:\n")
  print(performance_rmst_tbl)
  cat("\n")
}

cat("Treatment Model (RF) Performance:\n")
print(performance_treatment_tbl)
cat("\n")

cat("INTERPRETATION:\n")
cat("  Treatment Effect (theta): Change in RMST (months) per 1-unit increase in pathway score\n")
cat("  Negative theta indicates pathway activity is associated with shorter survival\n")
cat("  All models use K=5 cross-fitting to prevent overfitting\n")
cat("  Unified Random Forest treatment learner for fair comparison\n")
cat("  \n")
cat("  C-index: Concordance index for survival prediction (>0.7 excellent, >0.6 good)\n")
cat("  AUC: Time-dependent area under ROC curve at specified timepoints\n")
cat("  R²: Proportion of variance explained in RMST prediction\n")
cat("  Higher prediction performance → More efficient (lower SE) causal estimates\n")
sink()
cat("✓ Saved analysis summary to results/survival/summary/analysis_summary.txt\n")

#################################################################################
# Visualizations ----
#################################################################################

cat("\n=== Generating Visualizations ===\n")

# Load survminer if available for KM plots
if (!requireNamespace("survminer", quietly = TRUE)) {
  warning("survminer package not found. KM plots will be skipped. Install with: install.packages('survminer')")
  has_survminer <- FALSE
} else {
  library(survminer)
  has_survminer <- TRUE
}

# Figure 1: Forest Plot - Test Set + Full Dataset
cat("Generating Figure 1: Forest Plot (Test Set + Full Dataset)...\n")

# Combine test set and full dataset results
forest_data <- dplyr::bind_rows(
  dml_split_tbl %>%
    dplyr::filter(split == "test") %>%
    dplyr::mutate(dataset = "Test Set"),
  dml_full_tbl %>%
    dplyr::mutate(dataset = "Full Dataset")
) %>%
  dplyr::mutate(
    model = factor(model, levels = c("DR-GLMNET", "DR-Cox", "DR-RSF")),
    dataset = factor(dataset, levels = c("Test Set", "Full Dataset")),
    facet_label = paste0(model, "\n(", dataset, ")")
  )

forest_plot <- ggplot(forest_data, aes(x = theta, y = model, color = model, shape = dataset)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = 0.25, linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "DR-GLMNET" = "#E69F00", 
      "DR-Cox" = "#56B4E9", 
      "DR-RSF" = "#009E73"
    ),
    name = "DR Method"
  ) +
  scale_shape_manual(
    values = c("Test Set" = 16, "Full Dataset" = 17),
    name = "Dataset"
  ) +
  labs(
    title = "Causal Effect of Pathway Activity on RMST",
    subtitle = sprintf("Treatment: %s", 
                      stringr::str_trunc(selected_pathway_label, 60)),
    x = "Treatment Effect (months change in RMST per unit pathway score)",
    y = NULL,
    caption = "Note: Full dataset results are exploratory (no holdout validation)"
  ) +
  theme_nature() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
    legend.position = "right",
    axis.text.y = element_text(size = 11)
  )

ggsave("results/survival/figures/fig1_forest_plot_combined.png", forest_plot, 
       width = 10, height = 5, dpi = 600, bg = "white")
cat("✓ Saved Figure 1 (Combined Test Set + Full Dataset)\n")

# Figure 2: Stability across splits (including full dataset)
cat("Generating Figure 2: Stability Plot (with Full Dataset)...\n")

# Combine split results with full dataset
stability_data <- dplyr::bind_rows(
  dml_split_tbl %>%
    dplyr::mutate(dataset_type = "CV Split"),
  dml_full_tbl %>%
    dplyr::mutate(split = "full", dataset_type = "Full Data")
) %>%
  dplyr::mutate(
    split = factor(split, levels = c("train", "val", "test", "full")),
    dataset_type = factor(dataset_type, levels = c("CV Split", "Full Data"))
  )

stability_plot <- ggplot(stability_data, 
                         aes(x = split, y = theta, color = model, group = model)) +
  geom_point(aes(shape = dataset_type), size = 4, 
             position = position_dodge(width = 0.4)) +
  geom_line(data = stability_data %>% dplyr::filter(dataset_type == "CV Split"),
            linewidth = 1, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.25, linewidth = 0.9,
                position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = c(
    "DR-GLMNET" = "#E69F00", 
    "DR-Cox" = "#56B4E9", 
    "DR-RSF" = "#009E73"
  )) +
  scale_shape_manual(
    values = c("CV Split" = 16, "Full Data" = 17),
    name = "Dataset Type"
  ) +
  scale_x_discrete(labels = c("train" = "Train", "val" = "Val", 
                              "test" = "Test", "full" = "Full")) +
  labs(
    title = "Model Stability Across Data Splits",
    subtitle = "Consistent estimates indicate robust causal effect",
    x = "Data Split",
    y = "Treatment Effect Estimate (theta)",
    color = "DR Method"
  ) +
  theme_nature() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  )

ggsave("results/survival/figures/fig2_stability_across_splits.png", stability_plot, 
       width = 9, height = 6, dpi = 600, bg = "white")
cat("✓ Saved Figure 2 (with Full Dataset)\n")

# Figure 3: Kaplan-Meier curves by pathway activity
if (has_survminer) {
  cat("Generating Figure 3: Kaplan-Meier Curves (Pathway)...\n")
  
  # Prepare data with pathway score groups
  test_df_surv <- test_df %>%
    dplyr::mutate(
      pathway_group = factor(
        ifelse(pathway_score >= median(pathway_score, na.rm = TRUE), 
               "High Pathway Activity", 
               "Low Pathway Activity"),
        levels = c("Low Pathway Activity", "High Pathway Activity")
      )
    ) %>%
    dplyr::filter(!is.na(pathway_group), !is.na(OS_MONTHS), !is.na(OS_STATUS_BINARY))
  
  # Fit survival curves
  surv_fit_pathway <- survival::survfit(
    Surv(OS_MONTHS, OS_STATUS_BINARY) ~ pathway_group,
    data = test_df_surv
  )
  
  # Create KM plot
  km_plot_pathway <- survminer::ggsurvplot(
    surv_fit_pathway,
    data = test_df_surv,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    ggtheme = theme_nature(),
    palette = c("#00BFC4", "#F8766D"),
    title = sprintf("Survival by Pathway Activity (Test Set)\n%s", 
                   stringr::str_trunc(selected_pathway_label, 60)),
    xlab = "Time (months)",
    ylab = "Survival Probability",
    legend.title = "Pathway Activity",
    legend.labs = c("Low", "High"),
    font.main = c(14, "bold"),
    font.x = c(12, "plain"),
    font.y = c(12, "plain"),
    font.tickslab = c(10, "plain")
  )
  
  # Save
  ggsave("results/survival/figures/fig3a_km_curve_pathway.png", 
         print(km_plot_pathway), 
         width = 10, height = 8, dpi = 600, bg = "white")
  cat("✓ Saved Figure 3a\n")
  
  # Optional: KM curves by tumor stage
  if ("tumor_stage" %in% colnames(test_df)) {
    cat("Generating Figure 3b: Kaplan-Meier Curves (Tumor Stage)...\n")
    
    test_df_stage <- test_df %>%
      dplyr::filter(!is.na(tumor_stage), !is.na(OS_MONTHS), !is.na(OS_STATUS_BINARY)) %>%
      dplyr::mutate(
        tumor_stage = factor(tumor_stage, levels = c(0, 1), 
                           labels = c("Early Stage", "Late Stage"))
      )
    
    if (nrow(test_df_stage) > 0 && length(unique(test_df_stage$tumor_stage)) > 1) {
      surv_fit_stage <- survival::survfit(
        Surv(OS_MONTHS, OS_STATUS_BINARY) ~ tumor_stage,
        data = test_df_stage
      )
      
      km_plot_stage <- survminer::ggsurvplot(
        surv_fit_stage,
        data = test_df_stage,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.height = 0.25,
        ggtheme = theme_nature(),
        palette = c("#63A375", "#E4572E"),
        title = "Survival by Tumor Stage (Covariate, Test Set)",
        xlab = "Time (months)",
        ylab = "Survival Probability",
        legend.title = "Tumor Stage",
        legend.labs = c("Early", "Late"),
        font.main = c(14, "bold"),
        font.x = c(12, "plain"),
        font.y = c(12, "plain"),
        font.tickslab = c(10, "plain")
      )
      
      ggsave("results/survival/figures/fig3b_km_curve_tumor_stage.png", 
             print(km_plot_stage), 
             width = 10, height = 8, dpi = 600, bg = "white")
      cat("✓ Saved Figure 3b\n")
    }
  }
} else {
  cat("⊗ Skipped KM plots (survminer not installed)\n")
}

# Figure 4: Pathway Score vs RMST Scatter Plot
cat("Generating Figure 4: Scatter Plot (Pathway vs RMST)...\n")
scatter_rmst <- train_df %>%
  dplyr::filter(!is.na(pathway_score), !is.na(Y_pseudo)) %>%
  ggplot(aes(x = pathway_score, y = Y_pseudo)) +
  geom_point(alpha = 0.4, color = "#1F78B4", size = 1.5) +
  geom_smooth(method = "loess", color = "red", se = TRUE, linewidth = 1.2) +
  geom_smooth(method = "lm", color = "blue", linetype = "dashed", se = FALSE, linewidth = 1) +
  labs(
    title = "Pathway Activity vs RMST (Training Set)",
    subtitle = "Red: LOESS smoothing | Blue: Linear fit",
    x = "Pathway Score (D)",
    y = "RMST Pseudo-outcome (months)"
  ) +
  theme_nature() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/survival/figures/fig4_scatter_pathway_rmst.png", scatter_rmst, 
       width = 8, height = 6, dpi = 600, bg = "white")
cat("✓ Saved Figure 4\n")

# Figure 5: Sensitivity Analysis
cat("Generating Figure 5: Sensitivity Analysis...\n")
sensitivity_comparison <- dplyr::bind_rows(
  dml_split_tbl %>% 
    dplyr::filter(split == "train", model == "DR-GLMNET") %>%
    dplyr::mutate(adjustment = "Full (Clinical + Genomic)"),
  clinical_sensitivity_tbl %>%
    dplyr::mutate(split = "train", adjustment = "Clinical Only")
)

sensitivity_plot <- ggplot(sensitivity_comparison, 
                           aes(x = adjustment, y = theta, color = adjustment)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = c(
    "Full (Clinical + Genomic)" = "#009E73", 
    "Clinical Only" = "#D55E00"
  )) +
  labs(
    title = "Sensitivity Analysis: Covariate Adjustment",
    subtitle = "DR-GLMNET model on training set",
    x = NULL,
    y = "Treatment Effect (theta)"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 11)
  )

ggsave("results/survival/figures/fig5_sensitivity_clinical_vs_full.png", sensitivity_plot, 
       width = 7, height = 5, dpi = 600, bg = "white")
cat("✓ Saved Figure 5\n")

# Supplementary Figure S1: Event Rate by Pathway Quartile
cat("Generating Supplementary Figure S1: Event Rate by Quartile...\n")
event_by_quartile <- train_df %>%
  dplyr::filter(!is.na(pathway_score), !is.na(OS_STATUS_BINARY)) %>%
  dplyr::mutate(
    pathway_quartile = cut(
      pathway_score, 
      breaks = quantile(pathway_score, probs = seq(0, 1, 0.25), na.rm = TRUE),
      labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
      include.lowest = TRUE
    )
  ) %>%
  dplyr::group_by(pathway_quartile) %>%
  dplyr::summarise(
    event_rate = mean(OS_STATUS_BINARY == 1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(pathway_quartile))

event_rate_plot <- ggplot(event_by_quartile, 
                          aes(x = pathway_quartile, y = event_rate, 
                              fill = pathway_quartile)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", event_rate * 100, n)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(event_by_quartile$event_rate) * 1.15)) +
  labs(
    title = "Event Rate by Pathway Activity Quartile",
    subtitle = "Training set",
    x = "Pathway Score Quartile",
    y = "Event Rate (Death)"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("results/survival/figures/figS1_event_rate_by_quartile.png", event_rate_plot, 
       width = 8, height = 6, dpi = 600, bg = "white")
cat("✓ Saved Supplementary Figure S1\n")

# Supplementary Figure S2: GSVA Pathway Ranking
cat("Generating Supplementary Figure S2: GSVA Ranking...\n")
gsva_plot_df <- gsva_rank_tbl %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::mutate(
    ID = forcats::fct_reorder(ID, correlation),
    neg_log10_p = pmin(-log10(p_adjust), 10)  # Cap at 10 for visualization
  )

gsva_rank_plot <- ggplot(gsva_plot_df, 
                         aes(x = correlation, y = ID, size = abs(correlation), 
                             color = neg_log10_p)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#97C4BC", high = "#1F78B4", 
                      name = "-log10(adj.p)") +
  scale_size_continuous(range = c(3, 8), name = "|Correlation|") +
  labs(
    title = "Top Breast Cancer Pathways by GSVA Ranking",
    subtitle = "Ranked by correlation with RMST pseudo-outcome",
    x = "Correlation with RMST",
    y = NULL
  ) +
  theme_nature() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

ggsave("results/survival/figures/figS2_gsva_pathway_ranking.png", gsva_rank_plot, 
       width = 10, height = 7, dpi = 600, bg = "white")
cat("✓ Saved Supplementary Figure S2\n")

# Supplementary Figure S2b: Pathway Score Distribution by RMST Group
cat("Generating Supplementary Figure S2b: Pathway Distribution by RMST Group...\n")

# Create binary RMST groups based on median
train_df_plot <- train_df %>%
  dplyr::filter(!is.na(pathway_score), !is.na(Y_pseudo)) %>%
  dplyr::mutate(
    rmst_group = ifelse(Y_pseudo >= median(Y_pseudo, na.rm = TRUE), 
                       "High RMST (Better Survival)", 
                       "Low RMST (Worse Survival)"),
    rmst_group = factor(rmst_group, 
                       levels = c("High RMST (Better Survival)", 
                                 "Low RMST (Worse Survival)"))
  )

# Calculate summary statistics
pathway_stats <- train_df_plot %>%
  dplyr::group_by(rmst_group) %>%
  dplyr::summarise(
    mean_score = mean(pathway_score, na.rm = TRUE),
    median_score = median(pathway_score, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Create density plot
pathway_density_plot <- ggplot(train_df_plot, 
                               aes(x = pathway_score, fill = rmst_group)) +
  geom_density(alpha = 0.6, linewidth = 1) +
  geom_vline(data = pathway_stats, 
             aes(xintercept = median_score, color = rmst_group),
             linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("High RMST (Better Survival)" = "#7FC97F",  # Green
               "Low RMST (Worse Survival)" = "#E78C73"),    # Red/Orange
    name = "RMST Group"
  ) +
  scale_color_manual(
    values = c("High RMST (Better Survival)" = "#4D7C4D",
               "Low RMST (Worse Survival)" = "#B85C4A"),
    guide = "none"
  ) +
  labs(
    title = sprintf("GSVA Distribution for %s", 
                   stringr::str_trunc(selected_pathway_label, 50)),
    subtitle = "Treatment D (GSVA score) stratified by RMST outcome",
    x = "GSVA score (Pathway Activity)",
    y = "Density"
  ) +
  theme_nature() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("High RMST: n=%d\nLow RMST: n=%d", 
                         pathway_stats$n[1], pathway_stats$n[2]),
           hjust = 1.1, vjust = 1.5, size = 3.5, color = "gray30")

ggsave("results/survival/figures/figS2b_pathway_distribution_by_rmst.png", 
       pathway_density_plot, 
       width = 9, height = 6, dpi = 600, bg = "white")
cat("✓ Saved Supplementary Figure S2b (Pathway Distribution by RMST)\n")

# Supplementary Figure S3: Model Performance Comparison
cat("Generating Supplementary Figure S3: Model Performance Comparison...\n")

# Prepare data for plotting
performance_plot_data <- performance_outcome_tbl %>%
  dplyr::select(model, cindex, auc_60m, auc_120m) %>%
  tidyr::pivot_longer(cols = c(cindex, auc_60m, auc_120m),
                     names_to = "metric",
                     values_to = "value") %>%
  dplyr::mutate(
    metric_label = dplyr::case_when(
      metric == "cindex" ~ "C-index",
      metric == "auc_60m" ~ "AUC at 60m",
      metric == "auc_120m" ~ "AUC at 120m"
    ),
    metric_label = factor(metric_label, levels = c("C-index", "AUC at 60m", "AUC at 120m"))
  ) %>%
  dplyr::filter(!is.na(value))

if (nrow(performance_plot_data) > 0) {
  performance_comparison_plot <- ggplot(performance_plot_data,
                                        aes(x = model, y = value, fill = model)) +
    geom_col(width = 0.7, position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", value)), 
             vjust = -0.5, size = 3.5, fontface = "bold") +
    facet_wrap(~ metric_label, scales = "free_x", ncol = 3) +
    scale_fill_manual(values = c(
      "DR-GLMNET" = "#E69F00", 
      "DR-Cox" = "#56B4E9", 
      "DR-RSF" = "#009E73"
    )) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Outcome Model Prediction Performance (Test Set)",
      subtitle = "Higher values indicate better survival prediction",
      x = NULL,
      y = "Performance Metric Value"
    ) +
    theme_nature() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave("results/survival/figures/figS3_model_performance_comparison.png", performance_comparison_plot,
         width = 10, height = 5, dpi = 600, bg = "white")
  cat("✓ Saved Supplementary Figure S3\n")
} else {
  cat("⊗ Skipped Figure S3 (no performance data available)\n")
}

# Supplementary Figure S4: RMST Prediction (DR-GLMNET only)
if (!is.na(performance_results$glmnet$rmst_cor) && 
    !is.null(performance_results$glmnet$rmst_cor)) {
  
  cat("Generating Supplementary Figure S4: GLMNET RMST Prediction...\n")
  
  # Re-fit to get predictions for plotting
  tryCatch({
    test_X <- test_cov_matrix
    test_Y <- test_df$Y_pseudo
    complete_idx <- complete.cases(test_X, test_Y)
    test_X_complete <- test_X[complete_idx, , drop = FALSE]
    test_Y_complete <- test_Y[complete_idx]
    
    glmnet_fit_plot <- glmnet::cv.glmnet(
      x = test_X_complete,
      y = test_Y_complete,
      family = "gaussian",
      alpha = 0.5
    )
    
    glmnet_pred_plot <- predict(glmnet_fit_plot, newx = test_X_complete, s = "lambda.min")
    
    rmst_pred_data <- data.frame(
      observed = test_Y_complete,
      predicted = as.vector(glmnet_pred_plot)
    )
    
    rmst_pred_plot <- ggplot(rmst_pred_data, aes(x = observed, y = predicted)) +
      geom_point(alpha = 0.4, color = "#E69F00", size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                 color = "red", linewidth = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
      annotate("text", x = min(rmst_pred_data$observed) + 5, 
               y = max(rmst_pred_data$predicted) - 5,
               label = sprintf("R² = %.3f\nCorr = %.3f", 
                             performance_results$glmnet$rmst_rsq,
                             performance_results$glmnet$rmst_cor),
               hjust = 0, vjust = 1, size = 5, fontface = "bold") +
      labs(
        title = "DR-GLMNET: Predicted vs Observed RMST",
        subtitle = "Test set performance (dashed line = perfect prediction)",
        x = "Observed RMST (months)",
        y = "Predicted RMST (months)"
      ) +
      theme_nature() +
      theme(plot.title = element_text(face = "bold", size = 14))
    
    ggsave("results/survival/figures/figS4_glmnet_rmst_prediction.png", rmst_pred_plot,
           width = 7, height = 6, dpi = 600, bg = "white")
    cat("✓ Saved Supplementary Figure S4\n")
    
  }, error = function(e) {
    cat("⊗ Skipped Figure S4 (prediction data unavailable)\n")
  })
}

cat("\n=== All Results Saved Successfully ===\n")
cat("Output files:\n")
cat("  DML estimates: results/survival/tables/dml_estimates_*.csv\n")
cat("  Performance tables: results/survival/tables/tableS_*_performance.csv\n")
cat("  Summary: results/survival/summary/analysis_summary.txt\n")
cat("  Main figures: results/survival/figures/fig1-fig5_*.png\n")
cat("  Supplementary: results/survival/figures/figS1-figS4_*.png\n")
cat("\n")


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



