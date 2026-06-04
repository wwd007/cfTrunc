#' Conformal prediction for truncated time-to-event data
#'
#' @param data_tr Training data set.
#' @param data_ca Calibration data set.
#' @param data_te Test data set.
#' @param X Variable name for the observed event time.
#' @param Z Vector of variable names for covariates.
#' @param L Variable name for the left-truncation time.
#' @param R Variable name for the right-truncation or first sequential
#'   truncation time.
#' @param Rp Variable name for the second sequential truncation time.
#' @param delta Variable name for the event indicator.
#' @param trunc_type Truncation type: \code{"left"}, \code{"right"},
#'   \code{"double"}, or \code{"seq"}.
#' @param censoring Censoring type: \code{"right"} or \code{"none"}.
#' @param target Prediction target. Left truncation supports \code{"RMST"};
#'   other truncation types support \code{"MST"}.
#' @param tau Restriction time used for RMST. By default, the 90th percentile
#'   of observed training times.
#' @param outcome_model Outcome model: \code{"cox"}, \code{"aft"}, or
#'   \code{"rf"}. RF is supported for left truncation only.
#' @param truncation_model Truncation model: \code{"marginal"},
#'   \code{"reversed-cox"}, or \code{"rf"}. Dependent models are supported for
#'   left truncation only.
#' @param censoring_model Censoring model: \code{"none"}, \code{"marginal"},
#'   \code{"cox"}, or \code{"rf"}. Dependent models are supported for
#'   left-truncated RMST only.
#' @param b0 Reverse-time origin used by left-truncation models.
#' @param seed Random seed used by RF models.
#' @param mtry Number of variables sampled at each RF split. By default, the
#'   smaller of four and the number of covariates.
#' @param ntree Number of trees used by RF models.
#' @param eps Lower bound applied to estimated inclusion and censoring
#'   probabilities.
#' @param alpha Prespecified uncertainty level.
#'
#' @return A data frame with columns \code{y_pred}, \code{y_pred_hi}, and
#'   \code{y_pred_lo}.
#' @importFrom survival Surv
#' @export
#'
#' @examples
#' set.seed(42)
#' N <- 1500
#' Z <- runif(N, -1, 1)
#' T <- exp(0.2 + 0.5 * Z + rnorm(N, sd = 0.25))
#' C <- rexp(N, rate = 0.1)
#' L <- runif(N, 0, 0.5)
#' dat <- data.frame(X = pmin(T, C), Z, L, delta = as.integer(T <= C))
#' dat <- dat[dat$L < dat$X, ]
#' dat_tr <- dat[1:(nrow(dat) %/% 2), ]
#' dat_ca <- dat[(nrow(dat) %/% 2 + 1):nrow(dat), ]
#' tau <- unname(quantile(dat_tr$X, 0.9))
#'
#' Z_te <- runif(100, -1, 1)
#' T_te <- exp(0.2 + 0.5 * Z_te + rnorm(100, sd = 0.25))
#' dat_te <- data.frame(Z = Z_te)
#' pred <- conformal_pred(
#'   dat_tr, dat_ca, dat_te,
#'   X = "X", Z = "Z", L = "L", delta = "delta",
#'   trunc_type = "left", censoring = "right", target = "RMST",
#'   tau = tau, outcome_model = "cox", alpha = 0.1
#' )
#' mean(pred$y_pred_hi > pmin(T_te, tau) &
#'        pred$y_pred_lo < pmin(T_te, tau))
conformal_pred <- function(data_tr, data_ca, data_te,
                           X = "X", Z = "Z", L = "L", R = "R", Rp = "Rp",
                           delta = "delta",
                           trunc_type = "left",
                           censoring = "right",
                           target = "RMST",
                           tau = NA_real_,
                           outcome_model = "cox",
                           truncation_model = "marginal",
                           censoring_model = "marginal",
                           b0 = 100,
                           seed = 1,
                           mtry = NULL,
                           ntree = 100L,
                           eps = 1e-3,
                           alpha = 0.1) {
  trunc_type <- match_choice(trunc_type, c("left", "right", "double", "seq"),
                             "trunc_type")
  censoring <- match_choice(censoring, c("right", "none"), "censoring")
  target <- match_choice(target, c("RMST", "MST"), "target")
  outcome_model <- match_choice(outcome_model, c("cox", "aft", "rf"),
                                "outcome_model")
  truncation_model <- match_choice(
    truncation_model, c("marginal", "reversed-cox", "rf"),
    "truncation_model"
  )
  censoring_model <- match_choice(
    censoring_model, c("none", "marginal", "cox", "rf"),
    "censoring_model"
  )
  validate_scalar_probability(alpha, "alpha")
  validate_scalar_probability(eps, "eps")

  if (is.null(mtry)) {
    mtry <- min(4L, length(Z))
  }
  if (length(mtry) != 1L || !is.finite(mtry) || mtry < 1) {
    stop("mtry must be a positive number.", call. = FALSE)
  }
  if (length(ntree) != 1L || !is.finite(ntree) || ntree < 1) {
    stop("ntree must be a positive number.", call. = FALSE)
  }

  if (trunc_type == "left") {
    validate_left_configuration(censoring, target, censoring_model)
    return(conformal_pred_left(
      data_tr = data_tr, data_ca = data_ca, data_te = data_te,
      X = X, Z = Z, L = L, delta = delta,
      censoring = censoring, tau = tau,
      outcome_model = outcome_model,
      truncation_model = truncation_model,
      censoring_model = censoring_model,
      b0 = b0, seed = seed, mtry = mtry, ntree = ntree,
      eps = eps, alpha = alpha
    ))
  }

  validate_non_left_configuration(
    trunc_type, censoring, target, outcome_model,
    truncation_model, censoring_model
  )

  if (trunc_type == "right") {
    return(conformal_pred_right(
      data_tr, data_ca, data_te, X = X, Z = Z, R = R,
      outcome_model = outcome_model, eps = eps, alpha = alpha
    ))
  }
  if (trunc_type == "double") {
    return(conformal_pred_double(
      data_tr, data_ca, data_te, X = X, Z = Z, L = L, R = R,
      outcome_model = outcome_model, alpha = alpha
    ))
  }
  conformal_pred_seq(
    data_tr, data_ca, data_te, X = X, Z = Z, R = R, Rp = Rp,
    outcome_model = outcome_model, eps = eps, alpha = alpha
  )
}

validate_left_configuration <- function(censoring, target, censoring_model) {
  if (target != "RMST") {
    stop("Left truncation supports target = \"RMST\" only.", call. = FALSE)
  }
  if (censoring == "none" && censoring_model != "none") {
    stop("censoring_model must be \"none\" when censoring = \"none\".",
         call. = FALSE)
  }
  if (censoring == "right" && censoring_model == "none") {
    stop("censoring_model cannot be \"none\" when censoring = \"right\".",
         call. = FALSE)
  }
}

validate_non_left_configuration <- function(trunc_type, censoring, target,
                                            outcome_model, truncation_model,
                                            censoring_model) {
  if (censoring != "none" || censoring_model != "none") {
    stop(paste(trunc_type, "truncation requires censoring = \"none\" and",
               "censoring_model = \"none\"."), call. = FALSE)
  }
  if (target != "MST") {
    stop(paste(trunc_type, "truncation supports target = \"MST\" only."),
         call. = FALSE)
  }
  if (outcome_model == "rf") {
    stop("outcome_model = \"rf\" is supported for left truncation only.",
         call. = FALSE)
  }
  if (truncation_model != "marginal") {
    stop("Dependent truncation models are supported for left truncation only.",
         call. = FALSE)
  }
}

conformal_pred_left <- function(data_tr, data_ca, data_te,
                                X, Z, L, delta, censoring, tau,
                                outcome_model, truncation_model,
                                censoring_model, b0, seed, mtry, ntree,
                                eps, alpha) {
  require_columns(data_tr, c(X, Z, L), "data_tr")
  require_columns(data_ca, c(X, Z, L), "data_ca")
  require_columns(data_te, Z, "data_te")

  if (censoring == "right") {
    require_columns(data_tr, delta, "data_tr")
    require_columns(data_ca, delta, "data_ca")
  } else {
    data_tr[[delta]] <- 1L
    data_ca[[delta]] <- 1L
  }

  if (is.na(tau)) {
    tau <- unname(stats::quantile(data_tr[[X]], 0.9))
  }
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be a positive number.", call. = FALSE)
  }
  if (truncation_model %in% c("reversed-cox", "rf")) {
    validate_reverse_time_bound(data_tr, data_ca, X, L, b0)
  }

  truncation_fit <- fit_left_truncation_model(
    data_tr, X, Z, L, truncation_model, b0, seed, mtry, ntree
  )
  H_tr <- predict_left_inclusion_probability(
    truncation_fit, data_tr, X, L, b0, eps
  )
  H_ca <- predict_left_inclusion_probability(
    truncation_fit, data_ca, X, L, b0, eps
  )

  outcome_fit <- fit_left_outcome_model(
    data_tr, X, Z, L, delta, outcome_model, H_tr, seed, mtry, ntree
  )
  mu_ca <- predict_left_outcome_mean(
    outcome_fit, data_tr, data_ca, L, X, delta, tau
  )
  mu_te <- predict_left_outcome_mean(
    outcome_fit, data_tr, data_te, L, X, delta, tau
  )

  censoring_fit <- fit_left_censoring_model(
    data_tr, X, Z, L, delta, censoring_model, seed, mtry, ntree
  )
  S_cens_X <- predict_left_censoring_survival(
    censoring_fit, data_ca, data_ca[[X]], L, X, eps
  )
  S_cens_tau <- predict_left_censoring_survival(
    censoring_fit, data_ca, rep(tau, nrow(data_ca)), L, X, eps
  )

  censoring_weight <- data_ca[[delta]] / S_cens_X * (data_ca[[X]] <= tau) +
    1 / S_cens_tau * (data_ca[[X]] > tau)
  weights <- censoring_weight / H_ca
  scores <- abs(pmin(data_ca[[X]], tau) - mu_ca)
  prediction_interval(mu_te, scores, weights, alpha)
}

conformal_pred_right <- function(data_tr, data_ca, data_te,
                                 X, Z, R, outcome_model, eps, alpha) {
  require_columns(data_tr, c(X, Z, R), "data_tr")
  require_columns(data_ca, c(X, Z, R), "data_ca")
  require_columns(data_te, Z, "data_te")

  truncation_fit <- survival::survfit(
    survival::Surv(data_tr[[X]], data_tr[[R]], rep(1, nrow(data_tr))) ~ 1,
    timefix = FALSE
  )
  H_tr <- clamp_probability(
    predict_marginal_survival_at_times(truncation_fit, data_tr[[X]]), eps
  )
  H_ca <- clamp_probability(
    predict_marginal_survival_at_times(truncation_fit, data_ca[[X]]), eps
  )

  if (outcome_model == "aft") {
    fmla <- stats::reformulate(Z, response = paste0("log(", X, ")"))
    outcome_fit <- stats::lm(fmla, data = data_tr, weights = 1 / H_tr)
    mu_ca <- predict_aft_semipar(
      outcome_fit, data_ca, data = data_tr, trunc_type = "right", R = R
    )
    mu_te <- predict_aft_semipar(
      outcome_fit, data_te, data = data_tr, trunc_type = "right", R = R
    )
  } else {
    right_truncation_fmla <- stats::reformulate(Z, response = X)
    right_truncation_data <- data_tr
    right_truncation_data$.right_truncation <- right_truncation_data[[R]]
    right_truncation_fit <- coxrt::coxph.RT(
      right_truncation_fmla,
      right = .right_truncation,
      data = right_truncation_data
    )
    if (is.null(right_truncation_fit)) {
      stop("The right-truncated Cox model could not be estimated.", call. = FALSE)
    }
    beta <- right_truncation_fit$coef
    fmla <- make_surv_formula(stop = X, Z = Z)
    outcome_fit <- survival::coxph(fmla, data = data_tr, model = TRUE, x = TRUE)
    outcome_fit$coefficients <- beta
    mu_ca <- predict_cox_mean(outcome_fit, data_ca)
    mu_te <- predict_cox_mean(outcome_fit, data_te)
  }

  scores <- abs(data_ca[[X]] - mu_ca)
  prediction_interval(mu_te, scores, 1 / H_ca, alpha)
}

conformal_pred_double <- function(data_tr, data_ca, data_te,
                                  X, Z, L, R, outcome_model, alpha) {
  require_columns(data_tr, c(X, Z, L, R), "data_tr")
  require_columns(data_ca, c(X, Z, L, R), "data_ca")
  require_columns(data_te, Z, "data_te")

  npmle_tr <- cdfDT(data_tr[[X]], data_tr[[L]], data_tr[[R]], display = FALSE)
  npmle_ca <- cdfDT(data_ca[[X]], data_ca[[L]], data_ca[[R]], display = FALSE)
  data_tr$.outcome_weight <- 1 / npmle_tr$P.K
  if (outcome_model == "cox") {
    fmla <- make_surv_formula(stop = X, Z = Z)
    outcome_fit <- survival::coxph(
      fmla, data = data_tr, weights = .outcome_weight,
      model = TRUE, x = TRUE
    )
    mu_ca <- predict_cox_mean(outcome_fit, data_ca)
    mu_te <- predict_cox_mean(outcome_fit, data_te)
  } else {
    fmla <- stats::reformulate(Z, response = paste0("log(", X, ")"))
    outcome_fit <- stats::lm(fmla, data = data_tr, weights = .outcome_weight)
    mu_ca <- predict_aft_semipar(
      outcome_fit, data_ca, data = data_tr, trunc_type = "double", L = L, R = R
    )
    mu_te <- predict_aft_semipar(
      outcome_fit, data_te, data = data_tr, trunc_type = "double", L = L, R = R
    )
  }

  scores <- abs(data_ca[[X]] - mu_ca)
  prediction_interval(mu_te, scores, 1 / npmle_ca$P.K, alpha)
}

prediction_interval <- function(prediction, scores, weights, alpha) {
  q_star <- weighted_conformal_quantile(scores, weights, alpha)
  data.frame(
    y_pred = as.numeric(prediction),
    y_pred_hi = as.numeric(prediction + q_star),
    y_pred_lo = as.numeric(prediction - q_star)
  )
}
