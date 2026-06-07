match_choice <- function(value, choices, argument) {
  if (length(value) != 1L || !value %in% choices) {
    stop(
      paste0(argument, " must be one of: ", paste(choices, collapse = ", "), "."),
      call. = FALSE
    )
  }
  value
}

validate_scalar_probability <- function(value, argument) {
  if (length(value) != 1L || !is.finite(value) || value <= 0 || value >= 1) {
    stop(paste0(argument, " must be strictly between 0 and 1."), call. = FALSE)
  }
}

validate_scalar_logical <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(paste0(argument, " must be TRUE or FALSE."), call. = FALSE)
  }
}

require_columns <- function(data, columns, argument) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      paste0(argument, " is missing columns: ",
             paste(missing_columns, collapse = ", "), "."),
      call. = FALSE
    )
  }
}

make_surv_formula <- function(stop, Z, start = NULL, status = NULL) {
  surv_arguments <- c(start, stop, status)
  response <- paste0("survival::Surv(", paste(surv_arguments, collapse = ", "), ")")
  stats::reformulate(Z, response = response)
}

clamp_probability <- function(probability, eps) {
  pmax(pmin(as.numeric(probability), 1), eps)
}

step_survival_at_times <- function(event_times, survival, times) {
  index <- findInterval(times, event_times, left.open = TRUE)
  out <- rep(1, length(times))
  has_event <- index > 0L
  out[has_event] <- survival[index[has_event]]
  out
}

predict_marginal_survival_at_times <- function(fit, times) {
  step_survival_at_times(fit$time, fit$surv, times)
}

predict_cox_survival_at_times <- function(fit, newdata, times) {
  curve <- survival::survfit(fit, newdata = newdata)
  survival <- as.matrix(curve$surv)
  if (ncol(survival) == 1L && length(times) > 1L) {
    survival <- matrix(survival, nrow = nrow(survival), ncol = length(times))
  }
  vapply(seq_along(times), function(i) {
    step_survival_at_times(curve$time, survival[, i], times[i])
  }, numeric(1))
}

weighted_conformal_quantile <- function(scores, weights, alpha) {
  if (length(scores) != length(weights) || length(scores) == 0L) {
    stop("Calibration scores and weights must have the same positive length.",
         call. = FALSE)
  }
  if (any(!is.finite(scores)) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("Calibration scores and weights must be finite with positive total weight.",
         call. = FALSE)
  }
  ordering <- order(scores)
  standardized_weights <- weights[ordering] / sum(weights)
  index <- which(cumsum(standardized_weights) >= 1 - alpha)[1]
  scores[ordering][index]
}

validate_reverse_time_bound <- function(data_tr, data_ca, X, L, b0) {
  if (length(b0) != 1L || !is.finite(b0)) {
    stop("b0 must be a finite number.", call. = FALSE)
  }
  relevant_times <- c(data_tr[[X]], data_tr[[L]], data_ca[[X]], data_ca[[L]])
  if (any(!is.finite(relevant_times)) || any(relevant_times >= b0)) {
    stop("Training and calibration X and L values must be strictly below b0.",
         call. = FALSE)
  }
}

calc_winkler_score <- function(y, lo, hi, alpha) {
  ans <- hi - lo
  hi_ind <- y > hi
  lo_ind <- y < lo
  ans[hi_ind] <- ans[hi_ind] + 2 / alpha * (y[hi_ind] - hi[hi_ind])
  ans[lo_ind] <- ans[lo_ind] + 2 / alpha * (lo[lo_ind] - y[lo_ind])
  mean(ans)
}

conformal_result <- function(prediction, main_model, data_te, Z,
                             return_main_model, return_Zbhat) {
  if (!return_main_model && !return_Zbhat) {
    return(prediction)
  }

  result <- list(prediction = prediction)
  if (return_main_model) {
    result$main_model <- fitted_main_model(main_model)
  }
  if (return_Zbhat) {
    result$Zbhat <- predict_Zbhat(main_model, data_te, Z)
  }
  result
}

fitted_main_model <- function(main_model) {
  if (inherits(main_model, "cf_left_outcome_model")) {
    return(main_model$fit)
  }
  main_model
}

predict_Zbhat <- function(main_model, newdata, Z) {
  if (inherits(main_model, "cf_left_outcome_model")) {
    if (main_model$model == "rf") {
      stop("Zbhat is not available for random-forest outcome models.",
           call. = FALSE)
    }
    return(predict_Zbhat(main_model$fit, newdata, Z))
  }

  if (inherits(main_model, "cf_seq_trunc_cox") ||
      inherits(main_model, "cf_seq_trunc_aft")) {
    return(linear_predictor_from_named_beta(main_model$beta, newdata))
  }

  if (inherits(main_model, "lm") || inherits(main_model, "coxph")) {
    beta <- stats::coef(main_model)
    return(linear_predictor_from_terms(main_model, newdata, beta))
  }

  if (!is.null(main_model$coefficients) &&
      is.matrix(main_model$coefficients)) {
    beta <- main_model$coefficients[, 2]
    return(linear_predictor_from_named_beta(beta, newdata))
  }

  stop("Zbhat is not available for this outcome model.", call. = FALSE)
}

linear_predictor_from_terms <- function(model, newdata, beta) {
  beta <- beta[!is.na(beta)]
  if (is.null(names(beta)) || any(names(beta) == "")) {
    stop("The fitted main model does not have named coefficients.",
         call. = FALSE)
  }

  terms_object <- stats::delete.response(stats::terms(model))
  design <- stats::model.matrix(terms_object, data = newdata)
  missing_terms <- setdiff(names(beta), colnames(design))
  if (length(missing_terms) > 0L) {
    stop(
      paste0("data_te could not be expanded for terms: ",
             paste(missing_terms, collapse = ", "), "."),
      call. = FALSE
    )
  }
  as.numeric(design[, names(beta), drop = FALSE] %*% beta)
}

linear_predictor_from_named_beta <- function(beta, newdata) {
  beta <- beta[!is.na(beta)]
  if (is.null(names(beta)) || any(names(beta) == "")) {
    stop("The fitted main model does not have named coefficients.",
         call. = FALSE)
  }

  intercept_name <- intersect(c("(Intercept)", "Intercept"), names(beta))
  out <- rep(0, nrow(newdata))
  if (length(intercept_name) > 0L) {
    out <- out + beta[intercept_name[1L]]
  }

  term_names <- setdiff(names(beta), c("(Intercept)", "Intercept"))
  if (length(term_names) == 0L) {
    return(as.numeric(out))
  }
  require_columns(newdata, term_names, "data_te")
  numeric_terms <- vapply(newdata[, term_names, drop = FALSE], is.numeric,
                          logical(1))
  if (any(!numeric_terms)) {
    stop("Zbhat requires numeric covariates for this outcome model.",
         call. = FALSE)
  }
  design <- as.matrix(newdata[, term_names, drop = FALSE])
  as.numeric(out + design %*% beta[term_names])
}
