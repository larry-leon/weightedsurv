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
  NULL, no difference computed). The curves and the pooled survival are
  read at the last grid point `<= tzero`; before the first grid point
  the difference is 0 with variance 0, and at or after the last grid
  point the values at the last grid point are used.

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

With no subject weights or all subject weights equal the three variances
are the formulas used before, and these results are unchanged. When
`dfcounting` is a
[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
result whose subject weights are not all equal, the log-rank statistic,
the survival difference and their covariance rest on the same
per-subject contributions: `sig2.score` is \\\sum_i \psi_i(0)^2\\
(below; the same as `df_counting()$sig2_lr`), `sig2_dhat` is \\\sum_i
d_i^2\\ and `cov_wlr_dhat` is \\\sum_i \\-\psi_i(0)\\ d_i\\, where
\\d_i\\ is subject \\i\\'s contribution to the survival difference at
the grid point read for `tzero`: \\d_i = -\hat S_1 a_i \psi^{KM}\_i\\ in
the treatment arm and \\d_i = \hat S_0 a_i \psi^{KM}\_i\\ in the control
arm, with \\\psi^{KM}\_i\\ the Kaplan-Meier influence contribution of
that arm (see
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
section "Variance of the Kaplan-Meier estimate"). `sig2_dhat` then
equals `KM_diff(draws = 0)$sig2_dhat` at that grid point. `score` and
`dhat` are unchanged.

With a
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
result made with `weight_variance = "estimated"`, `sig2.score`,
`z.score`, `sig2_dhat`, `cov_wlr_dhat` and `cor_wlr_dhat` carry the
estimation of the weights (section "Estimated weights"); `sig2_dhat`
then equals `KM_diff(weight_variance = "estimated")$sig2_dhat` at that
grid point.

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

## Estimated weights

Propensity-score weights are estimated from the data. `weight_variance`
says how standard errors, intervals, bands and tests treat them.

- `"known"` (default):

  The weights are fixed, known numbers: the formulas of the two sections
  above. `weight_model` is ignored.

- `"estimated"`:

  The weights were estimated by a logistic regression for the treatment,
  given as `weight_model` (the object returned by
  [`ps_weights`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)),
  and that estimation is carried into the inference.

Estimates never change: survival curves, their difference, RMST, the
log-rank score, hazard ratios, counts and risk tables are identical
under the two options. Only standard errors, intervals, bands and tests
differ.

Notation. There are \\n\\ subjects, the rows analysed. Subject \\i\\ has
arm \\A_i\\, row \\x_i\\ of the model matrix of the propensity model
(\\p\\ columns), fitted probability \\e_i\\ and weight \\w_i\\ (written
\\a_i\\ in the score formulas above). The influence of subject \\i\\ on
the fitted coefficients and the derivative of its weight in the
coefficients are \$\$IF_i = \Big\\\sum_j x_j x_j' e_j (1 -
e_j)\Big\\^{-1} x_i (A_i - e_i), \qquad dw_i = -\frac{1 - e_i}{e_i} x_i
\\ (A_i = 1), \quad dw_i = \frac{e_i}{1 - e_i} x_i \\ (A_i = 0)\$\$ (for
stabilised weights \\w_i\\ and \\dw_i\\ are multiplied by the arm
proportion).

Kaplan-Meier estimate of a group \\g\\ (an arm; an arm within a stratum;
all subjects). With \\\psi_i(t)\\ the influence value of the section
"Variance of the Kaplan-Meier estimate", computed from the group's own
weighted counts, \$\$D_g(t) = \sum\_{k \in g} \psi_k(t)\\ dw_k, \qquad
\Phi_i(t) = 1(i \in g)\\ w_i \psi_i(t) + IF_i' D_g(t) \quad (i = 1,
\ldots, n),\$\$ \$\$\widehat{Var}\[\hat S_g(t)\] = \hat S_g(t)^2
\sum\_{i=1}^n \Phi_i(t)^2.\$\$ Every subject contributes, also those
outside the group, through the model. The variance of the difference of
the two arms is \$\$\widehat{Var}\[\hat S_1(t) - \hat S_0(t)\] =
\sum\_{i=1}^n \\\hat S_1(t) \Phi^1_i(t) - \hat S_0(t)
\Phi^0_i(t)\\^2,\$\$ which is **not** the sum of the two arm variances:
the arms share the model's coefficients. A resampled curve of arm \\a\\
is \\-\hat S_a(t) \sum\_{i=1}^n G_i \Phi^a_i(t)\\ with the
standard-normal multipliers \\G_i\\ that
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
generates in any case; the resampled difference curves, the simultaneous
band and the RMST variance and band of
[`cumulative_rmst_bands`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
follow from them unchanged.

Score of the log-rank test and of the weighted Cox estimate. With
\\\psi_i(\beta)\\ the contribution of the section "Variance of the
log-rank and Cox score with subject weights", computed within each
stratum, \$\$B(\beta) = \sum\_{j=1}^n \frac{\psi_j(\beta)}{a_j}\\ dw_j,
\qquad \psi^{ew}\_i(\beta) = \psi_i(\beta) + IF_i' B(\beta), \qquad
\widehat{Var}\\U(\beta)\\ = \sum\_{i=1}^n \psi^{ew}\_i(\beta)^2.\$\$ The
log-rank z-score, the variance of the log hazard ratio estimate (the sum
at \\\hat\beta\\ over the squared information) and a resampled score
(\\-\sum_i G_i \psi^{ew}\_i(\beta)\\) are formed as before with
\\\psi^{ew}\\ in place of \\\psi\\. In a stratified analysis the
variance is this single sum over all subjects, not a sum of per-stratum
variances. The cumulative variance of
[`wlr_cumulative`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
is \\\sum_i \\c_i(t) + IF_i' B(t)\\^2\\ with \\c_i(t)\\ the contribution
of subject \\i\\ to the score up to \\t\\ and \\B(t) = \sum_j \\c_j(t) /
a_j\\ dw_j\\; in `wlr_dhat_estimates` the contribution to the survival
difference is \\d_i + IF_i' (\hat S_0 D_0 - \hat S_1 D_1)\\ at the grid
point read for `tzero`, and the covariance uses both adjusted
contributions. The outputs of
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
that come from the `survival` package (`logrank_results`,
`cox_results$score`, `$se`, `$hr_ci`, `$pval`) carry the same
adjustment, applied to
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)'s weighted
score and dfbeta residuals \\r\\: \\r_i + IF_i' \sum_j (r_j / w_j)
dw_j\\.

How to use it. Fit the propensity model on exactly the rows to be
analysed with
[`ps_weights`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md),
put its weights in a column of `df`, and pass
`weight_variance = "estimated"` and `weight_model` to
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`cumulative_rmst_bands`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
or
[`plotKM.band_subgroups`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md).
[`cox_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`wlr_cumulative`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
and `wlr_dhat_estimates` take the option from the
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
result they are given (its elements `weight_variance` and `ew_parts`).
The rows of `df`, the treatment column and the weight column must be
those of `weight_model`; anything else is an error.

Scope and limits of this version.

- A logistic `glm` for a binary treatment, fitted on exactly the rows
  analysed. For weights from any other model use `"known"`.

- A subgroup analysis that reuses weights from a model fitted on the
  full sample must use `"known"`: the model's rows are then not the rows
  analysed. (The subgroup curves of
  [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)
  are point estimates and are not affected.)

- `draws > 0` (the resampling-based pointwise variance) is not available
  with `"estimated"`; `draws.band` and the resampling of
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  are.

- The time weights of the scheme (Fleming-Harrington and the others)
  depend on the weighted pooled curve; that dependence on the estimated
  weights is ignored.

- With stabilised weights the arm proportion is treated as a fixed
  number. The Kaplan-Meier quantities do not depend on that: they are
  unchanged when an arm's weights are multiplied by a constant.

- With `"estimated"` the influence-based formulas always run; whether
  the weights are all equal is not tested.

- Carrying the estimation through does not repair extreme weights: when
  the effective sample size (`df_counting()$ess`) is a small fraction of
  an arm, standard errors, intervals and tests can be too optimistic
  under either option.

The device of carrying the estimation of a nuisance quantity into the
inference by perturbing its estimate is that of Leon, Cai and Wei
(2009); it is used here in linearised form.
