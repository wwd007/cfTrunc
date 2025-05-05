predict_aft_semipar <- function(model, newdata, tau = NA) {
  if (class(model) == "lm") {
    exp_resid <- exp(model$residuals)
    exp_resid_surv <- survival::survfit(Surv(exp_resid) ~ 1, se.fit=F)
    if (!is.na(tau)) { # RMST
      exp_nlp <- exp(-predict(model, newdata))
      tau2 <- tau * exp_nlp
      pred <- rep(NA, nrow(newdata))
      for (i in 1:nrow(newdata)) {
        times <- c(0, exp_resid_surv$time[exp_resid_surv$time < tau2[i]], tau2[i])
        surv_probs <- c(1, exp_resid_surv$surv[exp_resid_surv$time < tau2[i]])
        pred[i] <- sum(diff(times) * surv_probs)
      }
      pred <- pred / exp_nlp
    } else { # MST
      times <- c(0, exp_resid_surv$time)
      surv_probs <- c(1, exp_resid_surv$surv[-length(exp_resid_surv$surv)])
      area <- sum(diff(times) * surv_probs)
      pred <- unname(exp(predict(model, newdata))*area)
    }
  } else if (class(model) == "aftgee") {
    linpred <- cbind(1, model$data$x) %*% model$coefficients[,2]

    time <- model$data$y * exp(-linpred)
    delta <- model$data$d
    exp_resid_surv <- survival::survfit(Surv(time, delta) ~ 1, se.fit=F)
    if (!is.na(tau)) { # RMST
      varnames <- dimnames(model$coefficients)[[1]][-1]
      Xmat <- as.matrix(newdata[,varnames,drop=F])
      linpred_new <- as.vector(cbind(1, Xmat) %*% model$coefficients[,2])
      exp_nlp <- exp(-linpred_new)
      tau2 <- tau * exp_nlp
      pred <- rep(NA, nrow(newdata))
      for (i in 1:nrow(newdata)) {
        times <- c(0, exp_resid_surv$time[exp_resid_surv$time < tau2[i]], tau2[i])
        surv_probs <- c(1, exp_resid_surv$surv[exp_resid_surv$time < tau2[i]])
        pred[i] <- sum(diff(times) * surv_probs)
      }
      pred <- pred / exp_nlp
    } else { # MST
      warning("Trying to calculate MST for censored data. MST may be biased.")
      times <- c(0, exp_resid_surv$time)
      surv_probs <- c(1, exp_resid_surv$surv[-length(exp_resid_surv$surv)])
      area <- sum(diff(times) * surv_probs)
      pred <- as.vector(exp(linpred)*area)
    }
  }
  return (pred)
}

predict_rmst <- function(cox_model, tau = NA, newdata) {
  surv_curve <- survival::survfit(cox_model, newdata = newdata)
  if (!is.na(tau)) {
    times <- c(0, surv_curve$time[surv_curve$time < tau], tau)
    surv_probs <- rbind(1, surv_curve$surv[surv_curve$time < tau,])
  } else {
    times <- c(0, surv_curve$time)
    surv_probs <- rbind(1, surv_curve$surv[-nrow(surv_curve$surv),])
  }
  area <- colSums(diff(times) * surv_probs)
  return (unname(area))
}


predict.aft.semipar.seqTrun <- function(model, newdata) {
  if (class(model) == "lm") {
    exp_resid <- exp(model$residuals)
    exp_resid_surv <- survfit(Surv(exp_resid) ~ 1, se.fit=F)
    if (!is.na(tau)) { # RMST
      exp_nlp <- exp(-predict(model, newdata))
      tau2 <- tau * exp_nlp
      pred <- rep(NA, nrow(newdata))
      for (i in 1:nrow(newdata)) {
        times <- c(0, exp_resid_surv$time[exp_resid_surv$time < tau2[i]], tau2[i])
        surv_probs <- c(1, exp_resid_surv$surv[exp_resid_surv$time < tau2[i]])
        pred[i] <- sum(diff(times) * surv_probs)
      }
      pred <- pred / exp_nlp
    } else { # MST
      times <- c(0, exp_resid_surv$time)
      surv_probs <- c(1, exp_resid_surv$surv[-length(exp_resid_surv$surv)])
      area <- sum(diff(times) * surv_probs)
      pred <- unname(exp(predict(model, newdata))*area)
    }
  }
  return (pred)
}


survfit.seqTrun <- function(cox_model, data, newdata) {
  # data: X, R, RR, Zj
  N <- nrow(data)
  beta <- cox_model$`Coefficient estimate`
  data <- data[order(data[["X"]]), ]
  km.fit <- survfit(Surv(R, RR, rep(1, nrow(data)))~1, data, se.fit=F)
  times <- data[["X"]]
  times_ind <- findInterval(times, km.fit$time)
  times_ind <- times_ind + 1
  surv_prob <- c(1, km.fit$surv)
  S_RR_X <- surv_prob[times_ind]
  exp_lin_pred <- exp(data[["Z1"]] * beta)
  denom <- exp_lin_pred / S_RR_X

  times_uniq <- unique(times)
  times_uniq_ind <- findInterval(times_uniq, times)
  times_uniq_ind <- c(0, times_uniq_ind)
  lambda_0 <- rep(NA, length(times_uniq))
  for (i in 1:length(times_uniq)) {
    lambda_0[i] <- sum(1/S_RR_X[(times_uniq_ind[i]+1):times_uniq_ind[i+1]]/
                         sum(denom[(times_uniq_ind[i]+1):N]))
  }
  Lambda_0 <- cumsum(lambda_0)
  surv_baseline <- exp(-Lambda_0)
  ans <- list(
    time = times_uniq,
    surv = outer(surv_baseline, exp(newdata$Z * beta), "^")
  )
  return (ans)
}

predict_rmst_seqTrun <- function(cox_model, tau = NA, data, newdata) {
  surv_curve <- survfit.seqTrun(cox_model, data, newdata)
  if (!is.na(tau)) {
    times <- c(0, surv_curve$time[surv_curve$time < tau], tau)
    surv_probs <- rbind(1, surv_curve$surv[surv_curve$time < tau,])
  } else {
    times <- c(0, surv_curve$time)
    surv_probs <- rbind(1, surv_curve$surv[-nrow(surv_curve$surv),])
  }
  area <- colSums(diff(times) * surv_probs)
  return (unname(area))
}

predict.aft.semipar.seqTrun <- function(model, data, newdata) {
  # model fitted by seqTrun.POReg.AFT
  # prediction by using NPMLE to estimate distribution of log(X)
  beta_hat <- model$mean[, 1]
  X_hat <- cbind(1, data$Z1) %*% beta_hat
  exp_res_logX <- exp(log(data$X) - X_hat)
  exp_res_logR <- exp(log(data$R) - X_hat)
  exp_res_logRR <- exp(log(data$RR) - X_hat)
  model_npmle <- seqTrunNPMLE_simp(
    data = cbind(exp_res_logX, exp_res_logR, exp_res_logRR),
    a = sort(exp_res_logX[,1]),
    t.a = max(exp_res_logX))
  # exp_resid_surv <- survfit(Surv(exp_resid) ~ 1, se.fit=F)
  times <- c(0, sort(exp_res_logX))
  surv_probs <- 1-model_npmle$`Point estimate`
  surv_probs <- c(1, surv_probs[-length(surv_probs)])
  area <- sum(diff(times) * surv_probs)
  pred <- unname(exp(cbind(1, newdata$Z1) %*% beta_hat)*area)
  return (pred)
}
