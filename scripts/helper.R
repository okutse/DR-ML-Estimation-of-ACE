# Helpers

## Helper functions for preprocessing
compute_qn_reference <- function(mat) {
  sorted <- apply(mat, 2, function(col) sort(col, na.last = TRUE))
  rowMeans(sorted, na.rm = TRUE)
}

apply_qn <- function(mat, reference) {
  ranked <- apply(mat, 2, function(col) rank(col, ties.method = "min", na.last = "keep"))
  normalized <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat),
                       dimnames = list(rownames(mat), colnames(mat)))
  for (j in seq_len(ncol(mat))) {
    idx <- which(!is.na(ranked[, j]))
    if (length(idx) > 0) {
      normalized[idx, j] <- reference[ranked[idx, j]]
    }
  }
  normalized
}

safe_log2 <- function(mat) log2(pmax(mat, 0) + 1)

zscore_with_stats <- function(mat, means, sds) {
  sweep(sweep(mat, 1, means, FUN = "-"), 1, sds, FUN = "/")
}

align_clinical_split <- function(dt, split_ids) {
  dt %>%
    dplyr::filter(patient_id %in% split_ids) %>%
    dplyr::arrange(match(patient_id, split_ids))
}


# PLR
## fit and predict with glmnet
# fit_glmnet_cv <- function(X, y, family = "gaussian", alpha = 1) {
#   X_mat <- as.matrix(X)
#   cvfit <- cv.glmnet(X_mat, y, family = family, alpha = alpha)
#   list(
#     model = cvfit,
#     predict = function(newX) {
#       predict(cvfit, newx = as.matrix(newX), s = "lambda.min")[,1]
#     }
#   )
# }
