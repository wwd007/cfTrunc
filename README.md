
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cfTrunc

<!-- badges: start -->
<!-- badges: end -->

The goal of `cfTrunc` is to generate robust prediction intervals using
a conformal prediction approach for time-to-event data subject to
truncation with or without censoring. The package supports left, right,
double, and sequential truncation.

For left-truncated RMST prediction, outcome models can be Cox, AFT, or
random forests. Left-truncation weights can use marginal, reversed-time
Cox, or random-forest models. Right-censoring weights can use marginal,
Cox, or random-forest models. These choices can be combined modularly.
Right, double, and sequential truncation currently support uncensored
MST prediction with Cox or AFT outcome models.

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
N <- 1500
Z <- runif(N, -1, 1)
T <- exp(0.2 + 0.5 * Z + rnorm(N, sd = 0.25))
C <- rexp(N, rate = 0.1)
L <- runif(N, 0, 0.5)
dat <- data.frame(X = pmin(T, C), Z, C, L)
dat$delta <- as.integer(T <= C)
dat <- dat[dat$L < dat$X, ]
dat_tr <- dat[1:(nrow(dat) %/% 2), ]
dat_ca <- dat[(nrow(dat) %/% 2 + 1):nrow(dat), ]
tau <- quantile(c(dat_tr$X), 0.9)

N_te <- 500
Z_te <- runif(N_te, -1, 1)
T_te <- exp(0.2 + 0.5 * Z_te + rnorm(N_te, sd = 0.25))
dat_te <- data.frame(Z = Z_te)

# conformal prediction
pred <- conformal_pred(dat_tr, dat_ca, dat_te,
                       X = "X", Z = c("Z"), L = "L", delta = "delta",
                       trunc_type = "left", censoring = "right",
                       target = "RMST", tau = tau,
                       outcome_model = "aft", alpha = 0.1)

# calculate coverage
coverage <- mean(pred$y_pred_hi > pmin(T_te, tau) &
                   pred$y_pred_lo < pmin(T_te, tau))
print(paste0("The coverage rate is ", coverage*100, "%."))
```

## Acknowledgment

The implementation of the distribution function estimation under double
truncation is adapted from `cdfDT()` function from `SurvTrunc` R
package.
