validate_seq_ordering <- function(data, X, R, Rp, argument) {
  if (any(!(data[[X]] < data[[R]] & data[[R]] < data[[Rp]]))) {
    stop(paste0(argument, " must satisfy ", X, " < ", R, " < ", Rp, "."),
         call. = FALSE)
  }
}

seq_trunc_second_survival <- function(data, R, Rp, times = data[[R]]) {
  fit <- survival::survfit(
    survival::Surv(data[[R]], data[[Rp]], rep(1, nrow(data))) ~ 1,
    timefix = FALSE
  )
  predict_marginal_survival_at_times(fit, times)
}

fit_seq_trunc_cox <- function(data, X, Z, R, Rp) {
  second_survival <- seq_trunc_second_survival(data, R, Rp)
  po_weight <- 1 / pmax(second_survival, 1e-4)
  thresholds <- unname(stats::quantile(data[[X]], seq(0.1, 0.8, length.out = 5)))
  long_data <- do.call(rbind, lapply(thresholds, function(threshold) {
    out <- data[, Z, drop = FALSE]
    out$.indicator <- as.integer(data[[X]] <= threshold)
    out$.threshold <- threshold
    out$.po_weight <- po_weight
    out$.id <- seq_len(nrow(data))
    out
  }))
  fmla <- stats::reformulate(c("factor(.threshold)", Z), response = ".indicator")
  fit <- geepack::geese(
    fmla, data = long_data, id = long_data$.id,
    weights = long_data$.po_weight,
    scale.fix = TRUE, family = stats::binomial(), jack = TRUE,
    mean.link = "cloglog", corstr = "independence"
  )
  structure(list(beta = fit$beta[Z]), class = "cf_seq_trunc_cox")
}

predict_seq_trunc_cox_mst <- function(object, data_tr, newdata, X, Z, R, Rp) {
  ordering <- order(data_tr[[X]])
  data_tr <- data_tr[ordering, , drop = FALSE]
  event_time <- data_tr[[X]]
  second_survival <- pmax(seq_trunc_second_survival(data_tr, R, Rp), 1e-4)
  risk <- as.vector(exp(as.matrix(data_tr[, Z, drop = FALSE]) %*% object$beta))
  times <- sort(unique(event_time))
  hazard <- vapply(times, function(time) {
    sum(1 / second_survival[event_time == time]) /
      sum(risk[event_time >= time] / second_survival[event_time >= time])
  }, numeric(1))
  baseline_survival <- exp(-cumsum(hazard))
  new_risk <- as.vector(exp(as.matrix(newdata[, Z, drop = FALSE]) %*% object$beta))
  curve <- list(time = times, surv = outer(baseline_survival, new_risk, "^"))
  integrate_survival_curve(curve)
}

fit_seq_trunc_aft <- function(data, X, Z, R, Rp) {
  second_survival <- seq_trunc_second_survival(data, R, Rp)
  data$.po_weight <- 1 / pmax(second_survival, 1e-4)
  fmla <- stats::reformulate(Z, response = paste0("log(", X, ")"))
  fit <- geepack::geese(
    fmla, data = data, id = seq_len(nrow(data)),
    weights = data$.po_weight,
    scale.fix = TRUE, family = stats::gaussian(), jack = TRUE,
    mean.link = "identity", corstr = "independence"
  )
  structure(
    list(beta = fit$beta[c("(Intercept)", Z)]),
    class = "cf_seq_trunc_aft"
  )
}

seq_trunc_npmle <- function(event_time, first_truncation, second_truncation) {
  n <- length(event_time)
  fit <- survival::survfit(
    survival::Surv(first_truncation, second_truncation, rep(1, n)) ~ 1,
    timefix = FALSE
  )
  second_survival <- pmax(
    predict_marginal_survival_at_times(fit, first_truncation), 1e-4
  )
  evaluation_points <- sort(event_time)
  cdf <- vapply(evaluation_points, function(point) {
    sum((event_time <= point) / second_survival) / sum(1 / second_survival)
  }, numeric(1))
  list(time = evaluation_points, cdf = cdf)
}

predict_seq_trunc_aft_mst <- function(object, data_tr, newdata, X, Z, R, Rp) {
  design_tr <- cbind("(Intercept)" = 1, as.matrix(data_tr[, Z, drop = FALSE]))
  linpred_tr <- as.vector(design_tr %*% object$beta)
  event_resid <- exp(log(data_tr[[X]]) - linpred_tr)
  first_resid <- exp(log(data_tr[[R]]) - linpred_tr)
  second_resid <- exp(log(data_tr[[Rp]]) - linpred_tr)
  npmle <- seq_trunc_npmle(event_resid, first_resid, second_resid)
  times <- c(0, npmle$time)
  survival <- c(1, (1 - npmle$cdf)[-length(npmle$cdf)])
  area <- sum(diff(times) * survival)

  design_new <- cbind("(Intercept)" = 1, as.matrix(newdata[, Z, drop = FALSE]))
  unname(exp(design_new %*% object$beta) * area)
}

conformal_pred_seq <- function(data_tr, data_ca, data_te,
                               X, Z, R, Rp, outcome_model, eps, alpha) {
  require_columns(data_tr, c(X, Z, R, Rp), "data_tr")
  require_columns(data_ca, c(X, Z, R, Rp), "data_ca")
  require_columns(data_te, Z, "data_te")
  validate_seq_ordering(data_tr, X, R, Rp, "data_tr")
  validate_seq_ordering(data_ca, X, R, Rp, "data_ca")

  if (outcome_model == "cox") {
    outcome_fit <- fit_seq_trunc_cox(data_tr, X, Z, R, Rp)
    mu_ca <- predict_seq_trunc_cox_mst(outcome_fit, data_tr, data_ca, X, Z, R, Rp)
    mu_te <- predict_seq_trunc_cox_mst(outcome_fit, data_tr, data_te, X, Z, R, Rp)
  } else {
    outcome_fit <- fit_seq_trunc_aft(data_tr, X, Z, R, Rp)
    mu_ca <- predict_seq_trunc_aft_mst(outcome_fit, data_tr, data_ca, X, Z, R, Rp)
    mu_te <- predict_seq_trunc_aft_mst(outcome_fit, data_tr, data_te, X, Z, R, Rp)
  }

  H_ca <- clamp_probability(seq_trunc_second_survival(data_tr, R, Rp, data_ca[[R]]), eps)
  scores <- abs(data_ca[[X]] - mu_ca)
  prediction_interval(mu_te, scores, 1 / H_ca, alpha)
}
