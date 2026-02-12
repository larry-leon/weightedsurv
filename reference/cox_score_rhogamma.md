# Cox score with rho-gamma weights

Calculates the Cox score statistic for weighted log-rank tests.

## Usage

``` r
cox_score_rhogamma(
  beta,
  time,
  delta,
  z,
  w_hat = rep(1, length(time)),
  wt_rg = rep(1, length(time)),
  score_only = TRUE
)
```

## Arguments

- beta:

  Log hazard ratio parameter.

- time:

  Numeric vector of event times.

- delta:

  Numeric vector of event indicators.

- z:

  Numeric vector of group indicators.

- w_hat:

  Numeric vector of weights.

- wt_rg:

  Numeric vector of rho-gamma weights.

- score_only:

  Logical; if TRUE, returns only the score.

## Value

Numeric value of the score or a list with additional results.
