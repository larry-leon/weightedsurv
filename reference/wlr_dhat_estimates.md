# Weighted Log-Rank and Difference Estimate at a Specified Time

Computes the weighted log-rank statistic, its variance, the difference
in survival at a specified time (`tzero`), the variance of the
difference, their covariance, and correlation, using flexible
time-dependent weights. The weighting scheme is selected via the
`scheme` argument and is calculated using `wt.rg.S`.

## Usage

``` r
wlr_dhat_estimates(
  dfcounting,
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 0),
  tzero = NULL
)
```

## Arguments

- dfcounting:

  List output from `df_counting` containing risk sets, event counts, and
  survival estimates.

- scheme:

  Character string specifying weighting scheme. One of: "fh"
  (Fleming-Harrington), "schemper", "XO", "MB", "custom_time", or
  "custom_code".

- scheme_params:

  Named list with numeric weighting parameters `rho` and `gamma` (used
  for "fh" and "custom_code" schemes).

- tzero:

  Time point at which to evaluate the difference in survival (default:
  24).

## Value

A list with elements:

- lr:

  Weighted log-rank test statistic.

- sig2_lr:

  Variance of the log-rank statistic.

- dhat:

  Difference in survival at `tzero`.

- cov_wlr_dhat:

  Covariance between log-rank and difference at `tzero`.

- sig2_dhat:

  Variance of the difference at `tzero`.

- cor_wlr_dhat:

  Correlation between log-rank and difference at `tzero`.

## Details

The weighting scheme is selected via the `scheme` argument and
calculated using `wt.rg.S`. Supports standard Fleming-Harrington,
Schemper, XO (Xu & O'Quigley), MB (Maggir-Burman), custom time-based,
and custom code weights.
