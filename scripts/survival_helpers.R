# Title: Survival outcome helper functions for IPCW-based RMST pseudo-outcome construction
# Author: Amos Okutse, Man-Fang Liang
# Date: Jan 2025

#' Compute censoring survival function G(t) using Kaplan-Meier or Cox model
#' 
#' @param time Observed time (OS_MONTHS)
#' @param event Event indicator (1 = event/death, 0 = censored)
#' @param X Optional covariate matrix for Cox censoring model
#' @param method "KM" for Kaplan-Meier or "Cox" for covariate-adjusted
#' @return Vector of G(t) values (censoring survival probability at each observed time)
compute_censoring_weights <- function(time, event, X = NULL, method = "KM") {
  library(survival)
  
  # Reverse event indicator: treat censoring as the "event"
  cens_indicator <- 1 - event
  
  if (method == "KM") {
    # Kaplan-Meier censoring model (marginal)
    surv_cens <- Surv(time = time, event = cens_indicator)
    km_fit <- survfit(surv_cens ~ 1)
    
    # Extract G(t) at each observed time
    # Use summary with extend=TRUE to handle times beyond last event
    km_summary <- summary(km_fit, times = time, extend = TRUE)
    G_t <- km_summary$surv
    
  } else if (method == "Cox" && !is.null(X)) {
    # Cox censoring model (covariate-adjusted)
    surv_cens <- Surv(time = time, event = cens_indicator)
    cox_data <- as.data.frame(X)
    cox_fit <- coxph(surv_cens ~ ., data = cox_data)
    
    # Predict survival function for each subject at their observed time
    surv_obj <- survfit(cox_fit, newdata = cox_data)
    
    # Extract G(t_i) for subject i at time t_i
    G_t <- numeric(length(time))
    for (i in seq_along(time)) {
      time_idx <- which.min(abs(surv_obj$time - time[i]))
      G_t[i] <- surv_obj$surv[time_idx, i]
    }
    
  } else {
    stop("method must be 'KM' or 'Cox'. For Cox, provide covariate matrix X.")
  }
  
  # Stabilization: set lower bound to prevent extreme weights
  G_t[G_t < 0.01] <- 0.01
  G_t[is.na(G_t)] <- 0.01
  
  return(G_t)
}


#' Compute RMST pseudo-outcome using IPCW
#' 
#' Formula: Y_pseudo = (Delta * min(T, tau) / G(T)) + ((1 - Delta / G(T)) * tau)
#' where Delta = event indicator, G(T) = censoring survival at observed time T
#' 
#' @param time Observed time (OS_MONTHS)
#' @param event Event indicator (1 = event/death, 0 = censored)
#' @param G_t Censoring survival function at each observed time
#' @param tau Time horizon for RMST
#' @param trim_quantile Quantile for trimming extreme weights (default 0.95)
#' @return Vector of RMST pseudo-outcomes
compute_rmst_pseudo <- function(time, event, G_t, tau, trim_quantile = 0.95) {
  
  # Restrict time to tau
  T_tilde <- pmin(time, tau)
  
  # IPCW weights: Delta / G(t)
  weights <- event / G_t
  
  # Trim extreme weights (only among events)
  if (sum(event == 1) > 0) {
    weight_threshold <- quantile(weights[event == 1], trim_quantile, na.rm = TRUE)
    weights[weights > weight_threshold] <- weight_threshold
  }
  
  # RMST pseudo-outcome formula
  # For events: contribute weighted survival time
  # For censored: contribute tau, downweighted
  Y_pseudo <- (weights * T_tilde) + ((1 - weights) * tau)
  
  # Ensure values are in valid range [0, tau]
  Y_pseudo[Y_pseudo < 0] <- 0
  Y_pseudo[Y_pseudo > tau] <- tau
  
  return(Y_pseudo)
}


#' Check survival data quality and censoring rate
#' 
#' @param time Observed time
#' @param event Event indicator
#' @param tau Time horizon
#' @return List with diagnostic statistics
check_survival_quality <- function(time, event, tau) {
  n <- length(time)
  n_events <- sum(event == 1, na.rm = TRUE)
  n_censored <- sum(event == 0, na.rm = TRUE)
  
  event_rate <- n_events / n
  censor_rate <- n_censored / n
  
  # Events by tau
  events_by_tau <- sum(event == 1 & time <= tau, na.rm = TRUE)
  event_rate_by_tau <- events_by_tau / n
  
  # Median follow-up (reverse KM)
  surv_cens <- Surv(time = time, event = 1 - event)
  km_cens <- survfit(surv_cens ~ 1)
  median_followup <- summary(km_cens)$table["median"]
  
  list(
    n_total = n,
    n_events = n_events,
    n_censored = n_censored,
    event_rate = event_rate,
    censor_rate = censor_rate,
    events_by_tau = events_by_tau,
    event_rate_by_tau = event_rate_by_tau,
    median_followup = median_followup,
    tau = tau
  )
}
