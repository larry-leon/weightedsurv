# Cumulative RMST bands for survival curves

Calculates cumulative Restricted Mean Survival Time (RMST) and
confidence bands for survival curves using resampling. Optionally plots
the cumulative RMST curve, pointwise confidence intervals, and
simultaneous confidence bands.

## Usage

``` r
cumulative_rmst_bands(
  df,
  fit,
  tte.name,
  event.name,
  treat.name,
  weight.name = NULL,
  draws_sb = 1000,
  xlab = "months",
  ylim_pad = 0.5,
  rmst_max_legend = "left",
  rmst_max_cex = 0.7,
  plot = TRUE,
  mean_rule = c("trapezoid", "rectangle"),
  weight_variance = c("known", "estimated"),
  weight_model = NULL
)
```

## Arguments

- df:

  Data frame containing survival data.

- fit:

  Survival fit object (output from KM_diff).

- tte.name:

  Name of time-to-event variable in `df`.

- event.name:

  Name of event indicator variable in `df`.

- treat.name:

  Name of treatment group variable in `df`.

- weight.name:

  Optional name of weights variable in `df`.

- draws_sb:

  Number of resampling draws for simultaneous bands (default: 1000).

- xlab:

  Label for x-axis (default: "months").

- ylim_pad:

  Padding for y-axis limits (default: 0.5).

- rmst_max_legend:

  Position for RMST legend (default: "left").

- rmst_max_cex:

  Text size for RMST legend (default: 0.7).

- plot:

  Logical; if TRUE, plot the results. Default is TRUE.

- mean_rule:

  Character; how the survival difference is integrated over the time
  grid of `fit` (`fit$at_points`). `"trapezoid"` (default) averages the
  difference at consecutive grid points: \\R(t_j) = \sum\_{k \le j}
  \\d(t\_{k-1}) + d(t_k)\\/2 \\ (t_k - t\_{k-1})\\. `"rectangle"` is the
  left-rectangle sum \\R(t_j) = \sum\_{k \le j} d(t\_{k-1}) (t_k -
  t\_{k-1})\\, the exact area under the step function. In both, \\R =
  0\\ at the first grid point. `"rectangle"` equals the difference of
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)'s
  restricted means (`summary(fit, rmean = tau)`) when the time grid
  holds every event time, which it does for the default
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  grid; it does not when the fit was made with `modify_tau = TRUE`. The
  rule applies to the estimate and to both sets of resampled curves;
  nothing else changes.

- weight_variance:

  How the weights are treated in standard errors, intervals, bands and
  tests. `"known"` (default): as fixed, known numbers. `"estimated"`: as
  estimated by the logistic propensity model given in `weight_model`,
  whose estimation is carried into the inference. Estimates never
  change. See the section "Estimated weights".

- weight_model:

  With `weight_variance = "estimated"`, the object returned by
  [`ps_weights`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
  for the rows of `df`, in the same order; its weights must be the
  column `weight.name`. Ignored with `"known"`.

## Value

A list with elements:

- at_points:

  Time points used for RMST calculation

- rmst_time:

  Cumulative RMST estimates

- sig2_rmst_time:

  Variance of RMST estimates

- rmst_time_lower:

  Pointwise lower confidence interval

- rmst_time_upper:

  Pointwise upper confidence interval

- rmst_maxtau_ci:

  RMST and CI at maximum time

- rmst_text:

  Text summary for legend

- c_alpha_band:

  Critical value for simultaneous band

- rmst_time_sb_lower:

  Simultaneous band lower bound

- rmst_time_sb_upper:

  Simultaneous band upper bound

- mean_rule:

  The integration rule used

- weight_variance:

  `"known"` or `"estimated"`: how the weights were treated

## Details

The pointwise variance of the cumulative RMST difference and the
simultaneous band come from the resampled survival-difference curves of
`fit` (`fit$dhat_star`) and of a second
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
call with `seedstart = 99999`. With weights that are not all equal,
those resampled curves are built from the per-subject contributions of
the robust variance (see
[`resampling_survival`](https://larry-leon.github.io/weightedsurv/reference/resampling_survival.md)),
so the RMST variance and bands are the robust ones; with no weights or
all weights equal they are unchanged. By default the weights are treated
as fixed and known. For weights estimated by a logistic propensity
model, `weight_variance = "estimated"` (with the same option in the
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
call that made `fit`; a `fit` made with the other option is an error)
carries their estimation into both sets of resampled curves, and so into
the variance, the intervals and the band (section "Estimated weights").
With extreme weights (low effective sample size) the intervals can be
too narrow under either option.

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
band and the RMST variance and band of `cumulative_rmst_bands` follow
from them unchanged.

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
a_j\\ dw_j\\; in
[`wlr_dhat_estimates`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
the contribution to the survival difference is \\d_i + IF_i' (\hat S_0
D_0 - \hat S_1 D_1)\\ at the grid point read for `tzero`, and the
covariance uses both adjusted contributions. The outputs of
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
`cumulative_rmst_bands` or
[`plotKM.band_subgroups`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md).
[`cox_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`wlr_cumulative`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
and
[`wlr_dhat_estimates`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
take the option from the
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
