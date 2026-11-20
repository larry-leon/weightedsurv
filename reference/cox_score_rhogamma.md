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

Numeric value of the score or, with `score_only = FALSE`, a list with
`score`, `sig2_score` (variance of the score at `beta`), `sig2_beta_asy`
(`sig2_score / i_bhat^2`), `K_wt_rg` and `i_bhat` (the information).

## Details

With `w_hat` absent or all equal, `sig2_score` and `sig2_beta_asy` are
the formulas used before
([`score_calculation`](https://larry-leon.github.io/weightedsurv/reference/score_calculation.md))
and are unchanged. When `w_hat` is not all equal, `sig2_score` is the
sum of squared per-subject score contributions \\\sum_i
\psi_i(\beta)^2\\ described below and `sig2_beta_asy` is that sum over
`i_bhat^2`. The score, `K_wt_rg` and `i_bhat` are unchanged. `time` is
taken to be sorted, as
[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
returns it; the time weight of a tied time is the one at its first
observation.

## Variance of the log-rank and Cox score with subject weights

With no subject weights, or with all subject weights equal, the weighted
log-rank test and the weighted Cox estimate are computed as before, and
their results are unchanged. When the subject weights are not all equal,
the variance of the weighted score is the sum of squared per-subject
score contributions, \$\$\widehat{Var}\\U(\beta)\\ = \sum_i
\psi_i(\beta)^2,\$\$ \$\$\psi_i(\beta) = a_i \Big( D_i K(U_i) \\Z_i -
\bar z(U_i, \beta)\\ - e^{\beta Z_i} \sum\_{s \le U_i} K(s) \\Z_i - \bar
z(s, \beta)\\ d\hat\Lambda(s, \beta) \Big),\$\$ where \\U_i\\, \\D_i\\,
\\Z_i\\ and \\a_i\\ are subject \\i\\'s time, event indicator, arm (0/1)
and subject weight, \\K(s)\\ is the time weight of the scheme, \\S_0(s,
\beta) = \sum_j a_j 1(U_j \ge s) e^{\beta Z_j}\\, \\\bar z(s, \beta) =
\sum_j a_j 1(U_j \ge s) Z_j e^{\beta Z_j} / S_0(s, \beta)\\ and
\\d\hat\Lambda(s, \beta) = \sum_j a_j D_j 1(U_j = s) / S_0(s, \beta)\\.
Sums over \\s\\ run over distinct times (events at a tied time are added
up), and a term whose \\S_0\\ is 0 is 0. The score itself, \\-\sum_i
\psi_i(\beta)\\, is the one computed before and does not change, and
neither does any point estimate. The log-rank z-score is the score at
\\\beta = 0\\ over \\\sqrt{\sum_i \psi_i(0)^2}\\; the variance of the
log hazard ratio estimate is \\\sum_i \psi_i(\hat\beta)^2 /
i(\hat\beta)^2\\ with \\i(\hat\beta)\\ the information (unchanged); in a
stratified analysis the contributions are computed within each stratum
and the score and its variance are sums over strata. For
Fleming-Harrington (0, 0) weights this is the robust variance of
`survival::coxph(weights = w, robust = TRUE, ties = "breslow")`: the
same robust score test and the same robust standard error. Results no
longer depend on how the weights are scaled (weights `w`, `10 * w` and
`w / sum(w)` give the same z-score and standard error). By default the
weights are treated as fixed and known. For weights estimated by a
logistic propensity model, `weight_variance = "estimated"` carries their
estimation into the test and the intervals (section "Estimated
weights"). When the weights are extreme (low effective sample size) the
intervals can be too short under either option.
