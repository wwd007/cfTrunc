#' Conformal prediction for truncated time-to-event data
#'
#' @param data_tr training data set
#' @param data_ca calibration data set
#' @param data_te test data set
#' @param X variable name for observed event time
#' @param Z vector of variable names for covariates
#' @param L variable name for left truncation time
#' @param R variable name for right truncation time
#' @param delta variable name for event indicator
#' @param trunc_type truncation type: left, right, or double
#' @param cencoring censoring type: right
#' @param target MST or RMST
#' @param tau used when targeting for RMST
#' @param model main model, cox or aft
#' @param lin_pred whether to output linear predictors from the main model
#' @param alpha prespecified uncertainty level, default is 0.1
#'
#' @return a data frame with 3 columns: y_pred, y_pred_hi, y_pred_lo
#' @import survival eha aftgee coxrt
#' @importFrom survival Surv
#'
#' @export
#'
#' @examples
#'# generate data
#'set.seed(42)
#'N <- 1000
#'Z <- runif(N)
#'X <- exp(Z)
#'C <- rexp(N)
#'L <- rgamma(N, 0.25, 0.25)
#'dat <- data.frame(X, Z, C, L)
#'dat$delta <- as.integer(dat$X > dat$C)
#'dat$X <- pmin(dat$X, dat$C)
#'dat <- dat[dat$L < dat$X, ]
#'dat_tr <- dat[1:(nrow(dat) %/% 2), ]
#'dat_ca <- dat[(nrow(dat) %/% 2 + 1):nrow(dat), ]
#'tau <- quantile(c(dat_tr$X), 0.9)
#'
#'N_te <- 500
#'Z <- runif(N_te)
#'X <- exp(Z)
#'C <- rexp(N_te)
#'L <- rgamma(N_te, 1, 1)
#'dat_te <- data.frame(X, Z, C, L)
#'dat_te$delta <- as.integer(dat_te$X > dat_te$C)
#'# conformal prediction
#'pred <- conformal_pred(dat_tr, dat_ca, dat_te,
#'                       X = "X", Z = c("Z"), L = "L", delta = "delta",
#'                       trunc_type = "left", target = "RMST",
#'                       tau = tau, model = "aft", alpha = 0.1)
#'# calculate coverage
#'coverage <- mean(pred$y_pred_hi > pmin(dat_te$X, tau) &
#'                   pred$y_pred_lo < pmin(dat_te$X, tau))
#'print(paste0("The coverage rate is ",coverage*100,"%."))
#'
conformal_pred <- function(data_tr, data_ca, data_te,
                           X = "X", Z = "Z", L = "L", R = "R", delta = "delta",
                           trunc_type = "left", cencoring = "right",
                           target = "RMST", tau = NA, model = "cox",
                           lin_pred = FALSE, alpha = 0.1) {
  if (target == "RMST") {
    if (is.na(tau)) tau <- unname(quantile(c(data_tr$X),0.9))
  }
  if (trunc_type == "left") {
    if (model == "cox") {
      # Surv(L, X, delta) ~ Z1 + Z2 + Z3
      fmla <- as.formula(paste("Surv(",L,",", X,",", delta,") ~ ",
                               paste(Z,collapse = "+")))
      model <- survival::coxph(fmla, data=data_tr)
      Zbhat <- as.vector(as.matrix(data_te[,Z]) %*% as.matrix(model$coefficients))
      mu_hat_tau_n1 <- predict_rmst(model, tau, newdata = data_ca)
      pred_data_te_mu <- predict_rmst(model, tau, newdata = data_te)
    } else if (model == "aft") {
      b0 <- 100 # an arbitary large number
      L_tilde <- b0-data_tr[[L]]
      X_tilde <- b0-data_tr[[X]]
      model.trunc.fit <- survival::survfit(Surv(X_tilde, L_tilde, rep(1,nrow(data_tr))) ~ 1,
                                 data=data_tr, timefix = FALSE) # see line 60
      trunc_time <- rev(b0-model.trunc.fit$time)
      trunc_prob <- c(1-rev(model.trunc.fit$surv), 0)

      t <- data_tr[[X]]
      # X_hat for truncation
      X_hat <- rep(NA, length(t))
      for (i in 1:length(t)) {
        temp_ind <- which(trunc_time<t[i])
        if (length(temp_ind)==0) {
          X_hat[i] <- 0.999 ## use a number close to 1, instead of 1, to avoid divided by 0.
        } else {
          X_hat[i] <- trunc_prob[max(temp_ind)+1]
          if(X_hat[i]==1) X_hat[i]<-0.999
        }
      }
      # Surv(X, delta) ~ Z1 + Z2 + Z3
      fmla <- as.formula(paste("Surv(", X,",", delta,") ~ ",
                               paste(Z,collapse = "+")))
      # model <- aftgee(fmla, data=data_tr, weights=1/(1-X_hat))
      invisible(capture.output(
        model <- aftgee::aftgee(fmla, data=data_tr, weights=1/(1-X_hat), B=0, binit="lm") ))
      Zbhat <- cbind(1, as.matrix(data_te[,Z])) %*% model$coefficients[,2]
      mu_hat_tau_n1 <- predict_aft_semipar(model, newdata = data_ca, tau=tau)
      pred_data_te_mu <- predict_aft_semipar(model, newdata = data_te, tau=tau)
      ### old
      # model <- aftreg(fmla, data=data_tr, dist=aft_dist)
      # mu_hat_tau_n1 <- predict.aftreg(model, newdata = data_ca[,Z], dist = aft_dist)
      # pred_data_te_mu <- predict.aftreg(model, newdata = data_te[,Z], dist = aft_dist)
    }
    ## model to estimate truncation time, fitted by training data
    # b0 <- max(c(data_tr$time,data_tr$L,data_ca$time,data_ca$L,data_te$time,data_te$L))+1
    b0 <- 100 # an arbitary large number
    L_tilde <- b0-data_tr[[L]]
    X_tilde <- b0-data_tr[[X]]
    model.trunc.fit <- survival::survfit(Surv(X_tilde, L_tilde, rep(1,nrow(data_tr))) ~ 1,
                               data=data_tr, timefix = FALSE)
    # Error in aeqSurv(Y) :
    #   aeqSurv exception, an interval has effective length 0
    # Probably the cause is the aeqSurv routine that treats time values such
    #  that tiny differences are treated as a tie. This is actually useful and
    # the error is potentially pointing an issue with the data.
    # However, if we need to force a solution you can use the coxph.options.
    # Just setting timefix = FALSE in the call to coxph should make the trick!
    #   https://stackoverflow.com/questions/46988426/error-in-aeqsurvy-aeqsurv-exception-an-interval-has-effective-length-0
    #   Source: https://rdrr.io/cran/survival/src/R/aeqSurv.R
    trunc_time <- rev(b0-model.trunc.fit$time)
    trunc_prob <- c(1-rev(model.trunc.fit$surv), 0)

    t <- data_ca[[X]]
    # X_hat for truncation
    X_hat <- rep(NA, length(t))
    for (i in 1:length(t)) {
      temp_ind <- which(trunc_time<t[i])
      if (length(temp_ind)==0) {
        X_hat[i] <- 0.999 ## use a number close to 1, instead of 1, to avoid divided by 0.
      } else {
        X_hat[i] <- trunc_prob[max(temp_ind)+1]
        if(X_hat[i]==1) X_hat[i]<-0.999
      }
    }
    # G_hat for censoring
    G <- survival::survfit(Surv(data_tr[[L]], data_tr[[X]], 1-data_tr[[delta]]) ~ 1,
                 se.fit=F, stype=2) # Breslow
    G_hat_T <- B(G, t=data_ca[[X]], newdata=data_ca)
    # w <- 1/ (1-X_hat) * (0 + data_ca[[delta]]/ ((1-G_hat_T)))
    w <- 1/ (1-X_hat) * ## 1/Pr(L<x)
      (data_ca[[delta]] / (1-G_hat_T)*(data_ca[[X]]<=tau) +
         1 / (1-G_hat_T)*(data_ca[[X]]>tau))
  }
  else if (trunc_type == "right") {
    model.trunc.fit <- survival::survfit(Surv(data_tr[[X]], data_tr[[R]], rep(1, nrow(data_tr)))~1)
    trunc_time <- rev(model.trunc.fit$time)
    trunc_prob <- c(1-rev(model.trunc.fit$surv), 0)
    if (model == "aft") {
      t <- data_tr[[X]]
      G_hat <- rep(NA, length(t))
      for (i in 1:length(t)) {
        # right truncation: find where trunc_time>t[i]
        temp_ind <- which(trunc_time>t[i])
        if (length(temp_ind)==0) {
          G_hat[i] <- 0.999 ## use a number close to 1, instead of 1, to avoid divided by 0.
        } else {
          G_hat[i] <- trunc_prob[max(temp_ind)+1]
          if(G_hat[i]==1) G_hat[i]<-0.999
        }
      }
      fmla <- as.formula(paste("log(",X,") ~", paste(Z, collapse = "+")))
      model <- lm(fmla, data=data_tr, weights=1/(1-G_hat))
      Zbhat <- cbind(1, as.matrix(data_te[,Z])) %*% model$coefficients
      mu_hat_tau_n1 <- predict_aft_semipar(model, newdata = data_ca)
      pred_data_te_mu <- predict_aft_semipar(model, newdata = data_te)
    } else if (model=="cox") {
      coxrt_beta <- coxrt:::.get_est(data_tr[[X]], data_tr[[R]], as.matrix(data_tr[,Z]),
                                     rep(1,nrow(data_tr)))[["est"]]
      fmla <- as.formula(paste("Surv(",X,") ~", paste(Z, collapse = "+")))
      model <- survival::coxph(fmla, data = data_tr)
      # hack
      model$coefficients <- coxrt_beta
      Zbhat <- as.vector(as.matrix(data_te[,Z]) %*% as.matrix(model$coefficients))
      # model <- coxph(Surv(L, time, event_1)~Z1, D_n1)
      mu_hat_tau_n1 <- predict_rmst(model, newdata = data_ca)
      pred_data_te_mu <- predict_rmst(model, newdata = data_te)
    }

    t <- data_ca[[X]]
    G_hat <- rep(NA, length(t))
    for (i in 1:length(t)) {
      # right truncation: find where trunc_time>t[i]
      temp_ind <- which(trunc_time>t[i])
      if (length(temp_ind)==0) {
        G_hat[i] <- 0.999 ## use a number close to 1, instead of 1, to avoid divided by 0.
      } else {
        G_hat[i] <- trunc_prob[max(temp_ind)+1]
        if(G_hat[i]==1) G_hat[i]<-0.999
      }
    }
    w <- 1/(1-G_hat) ## 1/Pr(R>x)
  }
  else if (trunc_type == "double") {
    model.NPMLE <- cdfDT(y=data_tr[[X]], l=data_tr[[L]], r=data_tr[[R]], display=F)
    if (model == "cox") {
      fmla <- as.formula(paste("Surv(",X,") ~", paste(Z, collapse = "+")))
      model <- survival::coxph(fmla, data = data_tr, weights = 1/model.NPMLE$P.K)
      Zbhat <- as.vector(as.matrix(data_te[,Z]) %*% as.matrix(model$coefficients))
      mu_hat_tau_n1 <- predict_rmst(model, newdata = data_ca)
      pred_data_te_mu <- predict_rmst(model, newdata = data_te)
    } else if (model == "aft") {
      fmla <- as.formula(paste("log(",X,") ~", paste(Z, collapse = "+")))
      model <- lm(fmla, data=data_tr, weights=1/model.NPMLE$P.K)
      Zbhat <- cbind(1, as.matrix(data_te[,Z])) %*% model$coefficients
      mu_hat_tau_n1 <- predict_aft_semipar(model, newdata = data_ca)
      pred_data_te_mu <- predict_aft_semipar(model, newdata = data_te)
    }
    model.n2.NPMLE <- cdfDT(y=data_ca[[X]], l=data_ca[[L]], r=data_ca[[R]], display=F)
    w <- 1/model.n2.NPMLE$P.K
  }

  if (target == "RMST") {
    R_star <- pmin(data_ca[[X]], tau)-mu_hat_tau_n1
  } else {
    R_star <- data_ca[[X]]-mu_hat_tau_n1
  }
  w_order_std <- w[order(R_star)] / sum(w)
  ind <- which.max(cumsum(w_order_std) >= 1 - alpha)
  q_star <- sort(R_star)[ind]

  if (target == "RMST") {
    pred_data_te <- data.frame(y_pred = pred_data_te_mu)
  } else {
    pred_data_te <- data.frame(y_pred = pmin(pred_data_te_mu, tau))
  }
  pred_data_te$y_pred_hi <- pred_data_te$y_pred + q_star
  pred_data_te$y_pred_lo <- pred_data_te$y_pred - q_star
  if (lin_pred) pred_data_te$Zbhat <- Zbhat
  return (pred_data_te)
}

# Surv(L, X, delta) ~ Z1+Z2+Z3
#
# Z=rnorm(100)
# X=exp(Z)
# fmla = as.formula(paste("Surv(X)~Z"))
# library(survival)
# model=survreg(fmla)
