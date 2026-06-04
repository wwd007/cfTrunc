fit_left_truncation_model <- function(data_tr, X, Z, L, model,
                                      b0, seed, mtry, ntree) {
  reverse_data <- data_tr
  reverse_data$.X_tilde <- b0 - reverse_data[[X]]
  reverse_data$.L_tilde <- b0 - reverse_data[[L]]
  reverse_data$.delta_tilde <- 1L

  if (model == "marginal") {
    fit <- survival::survfit(
      survival::Surv(.X_tilde, .L_tilde, .delta_tilde) ~ 1,
      data = reverse_data, timefix = FALSE
    )
  } else if (model == "reversed-cox") {
    fmla <- make_surv_formula(".L_tilde", Z, ".X_tilde", ".delta_tilde")
    fit <- survival::coxph(fmla, data = reverse_data, model = TRUE, x = TRUE)
  } else {
    reverse_data <- prepare_ltrcforest_data(reverse_data, Z)
    fmla <- make_surv_formula(
      ".L_tilde", ltrcforest_predictors(Z), ".X_tilde", ".delta_tilde"
    )
    set.seed(seed + 1)
    fit <- LTRCforests::ltrcrrf(
      formula = fmla, data = reverse_data, mtry = mtry, ntree = ntree
    )
  }
  structure(list(model = model, fit = fit), class = "cf_left_truncation_model")
}

predict_left_inclusion_probability <- function(object, newdata, X, L, b0, eps) {
  reverse_times <- b0 - newdata[[X]]
  if (object$model == "marginal") {
    probability <- predict_marginal_survival_at_times(object$fit, reverse_times)
  } else if (object$model == "reversed-cox") {
    probability <- predict_cox_survival_at_times(object$fit, newdata, reverse_times)
  } else {
    probability <- predict_surv_at_times_ltrcforest(
      object$fit, newdata, reverse_times,
      start_col = ".X_tilde", stop_col = ".L_tilde",
      status_col = ".delta_tilde"
    )
  }
  clamp_probability(probability, eps)
}

fit_left_outcome_model <- function(data_tr, X, Z, L, delta, model,
                                   inclusion_probability, seed, mtry, ntree) {
  if (model == "cox") {
    fmla <- make_surv_formula(X, Z, L, delta)
    fit <- survival::coxph(fmla, data = data_tr, model = TRUE, x = TRUE)
  } else if (model == "aft") {
    fmla <- make_surv_formula(X, Z, status = delta)
    data_tr$.outcome_weight <- 1 / inclusion_probability
    invisible(utils::capture.output(
      fit <- aftgee::aftgee(
        fmla, data = data_tr, weights = .outcome_weight,
        B = 0, binit = "lm"
      )
    ))
  } else {
    data_tr <- prepare_ltrcforest_data(data_tr, Z)
    fmla <- make_surv_formula(X, ltrcforest_predictors(Z), L, delta)
    set.seed(seed)
    fit <- LTRCforests::ltrcrrf(
      formula = fmla, data = data_tr, mtry = mtry, ntree = ntree
    )
  }
  structure(list(model = model, fit = fit), class = "cf_left_outcome_model")
}

predict_left_outcome_mean <- function(object, data_tr, newdata,
                                      L, X, delta, tau) {
  if (object$model == "cox") {
    return(predict_cox_mean(object$fit, newdata, tau))
  }
  if (object$model == "aft") {
    return(predict_aft_semipar(
      object$fit, newdata, tau, data = data_tr, trunc_type = "left", L = L
    ))
  }
  predict_rmst_ltrcforest(
    object$fit, newdata, tau, start_col = L, stop_col = X, status_col = delta
  )
}

fit_left_censoring_model <- function(data_tr, X, Z, L, delta, model,
                                     seed, mtry, ntree) {
  if (model == "none" || all(data_tr[[delta]] == 1)) {
    return(structure(list(model = "constant"), class = "cf_left_censoring_model"))
  }
  data_cens <- data_tr
  data_cens$.delta_cens <- 1 - data_cens[[delta]]
  if (model == "marginal") {
    fit <- survival::survfit(
      survival::Surv(data_cens[[L]], data_cens[[X]], data_cens$.delta_cens) ~ 1,
      se.fit = FALSE, stype = 2
    )
  } else if (model == "cox") {
    fmla <- make_surv_formula(X, Z, L, ".delta_cens")
    fit <- survival::coxph(fmla, data = data_cens, model = TRUE, x = TRUE)
  } else {
    data_cens <- prepare_ltrcforest_data(data_cens, Z)
    fmla <- make_surv_formula(X, ltrcforest_predictors(Z), L, ".delta_cens")
    set.seed(seed + 2)
    fit <- LTRCforests::ltrcrrf(
      formula = fmla, data = data_cens, mtry = mtry, ntree = ntree
    )
  }
  structure(list(model = model, fit = fit), class = "cf_left_censoring_model")
}

predict_left_censoring_survival <- function(object, newdata, times, L, X, eps) {
  if (object$model == "constant") {
    return(rep(1, nrow(newdata)))
  }
  if (object$model == "marginal") {
    probability <- predict_marginal_survival_at_times(object$fit, times)
  } else if (object$model == "cox") {
    probability <- predict_cox_survival_at_times(object$fit, newdata, times)
  } else {
    probability <- predict_surv_at_times_ltrcforest(
      object$fit, newdata, times,
      start_col = L, stop_col = X, status_col = ".delta_cens"
    )
  }
  clamp_probability(probability, eps)
}

ltrcforest_predictors <- function(Z) {
  if (length(Z) == 1L) c(Z, ".cf_rf_duplicate") else Z
}

prepare_ltrcforest_data <- function(data, Z) {
  if (length(Z) == 1L) {
    data$.cf_rf_duplicate <- data[[Z]]
  }
  data
}
