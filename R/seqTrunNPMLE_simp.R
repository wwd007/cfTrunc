seqTrunNPMLE_simp <- function (data, a, t.a, B, V) 
{
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[, 1]
  x <- temp.dat[, 2]
  r <- temp.dat[, 3]
  w <- rep(1, length(x))
  # npmle.temp <- matrix(0, B, length(a))
  # beta <- seqTrun:::fn_est(1, l, x, r, w)[[2]]
  # point_est <- sapply(a, seqTrun:::fun_test, l = l, x = x, r = r, 
  #                     w = w)
  
  n <- length(x)
  fit <- survival::survfit(survival::Surv(x, r, event = rep(1, n)) ~ 1, 
                           weights = w, timefix = FALSE)
  s.r.x <- summary(fit, time = x)$surv
  s.r.x <- s.r.x[match(x, sort(x))]
  point_est <- 1/sum(w/s.r.x) * sapply(a, function(a) sum((l <= a) * w/s.r.x))
  # if (V == TRUE) {
  #   F.l.cov.temp <- sapply(sort(a), fun_covar_est, var2 = t.a, 
  #                          l = l, x = x, r = r)
  #   F.l.variance.temp <- sapply(sort(a), fun_variance_est, 
  #                               var2 = t.a, l = l, x = x, r = r, w = w)
  #   F.l.variance.temp <- ifelse(F.l.variance.temp > 0, F.l.variance.temp, 
  #                               0)
  #   F.l.se.temp.new <- sqrt(F.l.variance.temp)
  # }
  # if (B > 1) {
  #   for (j in 1:B) {
  #     set.seed(j)
  #     w <- 4 * rbeta(length(x), 0.5, 1.5)
  #     w <- w/sum(w)
  #     F.npmle.est.bp <- sapply(sort(a), fun_test, l = l, 
  #                              x = x, r = r, w = w)
  #     F.npmle.est.bp <- F.npmle.est.bp/F.npmle.est.bp[length(a)]
  #     npmle.temp[j, ] <- F.npmle.est.bp
  #   }
  #   F.np.se.temp <- apply(npmle.temp, 2, sd)
  # }
  # out <- list(a, point_est, F.np.se.temp, F.l.se.temp.new, 
  #             beta)
  # names(out) <- c("Evaluation point", "Point estimate", "Bootstrap SE", 
  #                 "Analytical SE", "Truncation probability")
  out <- list(a, point_est)
  names(out) <- c("Evaluation point", "Point estimate")
  return(out)
}
