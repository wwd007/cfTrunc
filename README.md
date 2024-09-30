
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cfTrunc

<!-- badges: start -->
<!-- badges: end -->

The goal of `cfTrunc` is to generate robust prediction intervals using
conformal prediction approach for time-to-event data subject to
truncation with or without censoring. Currently, `cfTrunc` supports
left, right, and double truncation. More truncation types will be added
later.

## Installation

You can install the development version of cfTrunc like so:

``` r
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("wwd007/cfTrunc")
```

## Example

This is a basic example which shows you how to generate prediction
intervals for left truncated and right censored data:

``` r
# load the package
library(cfTrunc)

# generate data
set.seed(42)
N <- 1000
Z <- runif(N)
X <- exp(Z)
C <- rexp(N)
L <- rgamma(N, 0.25, 0.25)
dat <- data.frame(X, Z, C, L)
dat$delta <- as.integer(dat$X > dat$C)
dat$X <- pmin(dat$X, dat$C)
dat <- dat[dat$L < dat$X, ]
dat_tr <- dat[1:(nrow(dat) %/% 2), ]
dat_ca <- dat[(nrow(dat) %/% 2 + 1):nrow(dat), ]
tau <- quantile(c(dat_tr$X), 0.9)

N_te <- 500
Z <- runif(N_te)
X <- exp(Z)
C <- rexp(N_te)
L <- rgamma(N_te, 1, 1)
dat_te <- data.frame(X, Z, C, L)
dat_te$delta <- as.integer(dat_te$X > dat_te$C)

# conformal prediction
pred <- conformal_pred(dat_tr, dat_ca, dat_te,
                       X = "X", Z = c("Z"), L = "L", delta = "delta",
                       trunc_type = "left", target = "RMST",
                       tau = tau, model = "aft", alpha = 0.1)

# calculate coverage
coverage <- mean(pred$y_pred_hi > pmin(dat_te$X, tau) &
                   pred$y_pred_lo < pmin(dat_te$X, tau))
print(paste0("The coverage rate is ", coverage*100, "%."))
```

## Acknowledgment

The implementation of the distribution function estimation under double
truncation is adapted from `cdfDT()` function from `SurvTrunc` R
package.
