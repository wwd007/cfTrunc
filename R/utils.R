B <- function(cox_model, t, newdata) {
  N <- nrow(newdata)
  ans <- rep(NA, N)
  if ("coxph" %in% class(cox_model)) {
    surv_curve <- survival::survfit(cox_model, newdata = newdata, type="breslow")
    for (i in 1:N) {
      temp <- which(surv_curve$time<t[i])
      if (length(temp)==0) {
        ans[i] <- 0
      } else {
        ans[i] <- 1-surv_curve$surv[max(temp), i]
      }
    }
  } else if ("survfit" %in% class(cox_model)) {
    surv_curve <- cox_model
    for (i in 1:N) {
      temp <- which(surv_curve$time<t[i])
      if (length(temp)==0) {
        ans[i] <- 0
      } else {
        ans[i] <- 1-surv_curve$surv[max(temp)]
      }
    }
  }
  return (ans)
}

calc_winkler_score <- function(y, lo, hi, alpha) {
  ans <- hi - lo
  hi_ind <- y > hi
  lo_ind <- y < lo
  ans[hi_ind] <- ans[hi_ind] + 2/alpha*(y[hi_ind] - hi[hi_ind])
  ans[lo_ind] <- ans[lo_ind] + 2/alpha*(lo[lo_ind] - y[lo_ind])
  return (mean(ans))
}
