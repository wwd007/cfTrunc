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
