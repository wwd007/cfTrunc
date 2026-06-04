integrate_survival_curve <- function(curve, tau = NA_real_) {
  survival <- as.matrix(curve$surv)
  if (!is.na(tau)) {
    keep <- curve$time < tau
    times <- c(0, curve$time[keep], tau)
    survival <- rbind(1, survival[keep, , drop = FALSE])
  } else {
    times <- c(0, curve$time)
    survival <- rbind(1, survival[-nrow(survival), , drop = FALSE])
  }
  unname(colSums(diff(times) * survival))
}

predict_cox_mean <- function(model, newdata, tau = NA_real_) {
  curve <- survival::survfit(model, newdata = newdata)
  integrate_survival_curve(curve, tau)
}

predict_aft_semipar <- function(model, newdata, tau = NA_real_, data,
                                trunc_type = "left", L = "L", R = "R") {
  if (inherits(model, "lm")) {
    exp_resid <- exp(model$residuals)
    n <- length(exp_resid)
    if (trunc_type == "left") {
      exp_resid_L <- exp(log(data[[L]]) - model$fitted.values)
      exp_resid_surv <- survival::survfit(
        survival::Surv(exp_resid_L, exp_resid, rep(1, n)) ~ 1,
        se.fit = FALSE
      )
    } else if (trunc_type == "right") {
      exp_resid_R <- exp(log(data[[R]]) - model$fitted.values)
      exp_resid_surv <- survival::survfit(
        survival::Surv(-exp_resid_R, -exp_resid, rep(1, n)) ~ 1,
        se.fit = FALSE
      )
      exp_resid_surv$time <- rev(-exp_resid_surv$time)
      exp_resid_surv$surv <- rev(1 - exp_resid_surv$surv)
    } else {
      exp_resid_L <- exp(log(data[[L]]) - model$fitted.values)
      exp_resid_R <- exp(log(data[[R]]) - model$fitted.values)
      npmle_fit <- cdfDT(exp_resid, exp_resid_L, exp_resid_R, display = FALSE)
      exp_resid_surv <- list(time = npmle_fit$time, surv = npmle_fit$Survival)
    }

    if (!is.na(tau)) {
      exp_nlp <- exp(-stats::predict(model, newdata))
      tau_scaled <- tau * exp_nlp
      pred <- vapply(seq_len(nrow(newdata)), function(i) {
        keep <- exp_resid_surv$time < tau_scaled[i]
        times <- c(0, exp_resid_surv$time[keep], tau_scaled[i])
        surv_probs <- c(1, exp_resid_surv$surv[keep])
        sum(diff(times) * surv_probs)
      }, numeric(1))
      return(pred / exp_nlp)
    }

    times <- c(0, exp_resid_surv$time)
    surv_probs <- c(1, exp_resid_surv$surv[-length(exp_resid_surv$surv)])
    area <- sum(diff(times) * surv_probs)
    return(unname(exp(stats::predict(model, newdata)) * area))
  }

  linpred <- cbind(1, model$data$x) %*% model$coefficients[, 2]
  exp_resid <- model$data$y * exp(-linpred)
  exp_resid_L <- exp(log(data[[L]]) - linpred)
  exp_resid_surv <- survival::survfit(
    survival::Surv(exp_resid_L, exp_resid, model$data$d) ~ 1,
    se.fit = FALSE
  )

  varnames <- dimnames(model$coefficients)[[1]][-1]
  Xmat <- as.matrix(newdata[, varnames, drop = FALSE])
  linpred_new <- as.vector(cbind(1, Xmat) %*% model$coefficients[, 2])
  exp_nlp <- exp(-linpred_new)
  tau_scaled <- tau * exp_nlp
  pred <- vapply(seq_len(nrow(newdata)), function(i) {
    keep <- exp_resid_surv$time < tau_scaled[i]
    times <- c(0, exp_resid_surv$time[keep], tau_scaled[i])
    surv_probs <- c(1, exp_resid_surv$surv[keep])
    sum(diff(times) * surv_probs)
  }, numeric(1))
  pred / exp_nlp
}

predict_rmst_ltrcforest <- function(object, newdata, tau,
                                    start_col = "L", stop_col = "X",
                                    status_col = "delta") {
  time_points <- seq(0, tau, length.out = 50)
  x_var <- names(object$xvar)
  newdata <- add_ltrcforest_dummy(newdata, x_var)
  nd <- newdata[, x_var, drop = FALSE]
  nd[[start_col]] <- 0
  nd[[stop_col]] <- tau
  nd[[status_col]] <- 1

  pred <- LTRCforests::predictProb(
    object = object, newdata = nd, time.eval = time_points
  )
  survival <- forest_survival_matrix(
    pred$survival.probs, nrow(newdata), length(pred$survival.times)
  )
  interval_length <- diff(pred$survival.times)
  colSums(survival[-nrow(survival), , drop = FALSE] * interval_length)
}

predict_surv_at_times_ltrcforest <- function(object, newdata, times,
                                             start_col, stop_col, status_col,
                                             start_value = 0) {
  times <- as.numeric(times)
  times_eval <- sort(unique(times))
  x_var <- names(object$xvar)
  newdata <- add_ltrcforest_dummy(newdata, x_var)
  nd <- newdata[, x_var, drop = FALSE]
  nd[[start_col]] <- start_value
  nd[[stop_col]] <- max(times_eval)
  nd[[status_col]] <- 1

  pred <- LTRCforests::predictProb(
    object = object, newdata = nd, time.eval = times_eval
  )
  survival <- forest_survival_matrix(
    pred$survival.probs, nrow(newdata), length(times_eval)
  )
  index <- match(times, times_eval)
  vapply(seq_len(nrow(newdata)), function(i) survival[index[i], i], numeric(1))
}

forest_survival_matrix <- function(survival, n_observations, n_times) {
  if (is.null(dim(survival))) {
    return(matrix(survival, nrow = n_times, ncol = n_observations))
  }
  survival <- as.matrix(survival)
  if (nrow(survival) == n_times && ncol(survival) == n_observations) {
    return(survival)
  }
  if (nrow(survival) == n_observations && ncol(survival) == n_times) {
    return(t(survival))
  }
  stop("Unexpected survival-probability dimensions returned by LTRCforests.",
       call. = FALSE)
}

add_ltrcforest_dummy <- function(data, predictors) {
  if (".cf_rf_duplicate" %in% predictors &&
      !".cf_rf_duplicate" %in% names(data)) {
    source_predictor <- setdiff(predictors, ".cf_rf_duplicate")[1]
    data$.cf_rf_duplicate <- data[[source_predictor]]
  }
  data
}
