# Plot Kaplan-Meier Survival Difference Curves with Subgroups and Confidence Bands

Plots the difference in Kaplan-Meier survival curves between two groups,
optionally including simultaneous confidence bands and subgroup curves.
Also displays risk tables for the overall population and specified
subgroups.

## Usage

``` r
plotKM.band_subgroups(
  df,
  tte.name,
  event.name,
  treat.name,
  weight.name = NULL,
  sg_labels = NULL,
  ltype = "s",
  lty = 1,
  draws = 0,
  lwd = 2,
  sg_colors = NULL,
  color = "lightgrey",
  ymax.pad = 0.01,
  ymin.pad = -0.01,
  taus = c(-Inf, Inf),
  yseq_length = 5,
  cex_Yaxis = 0.8,
  risk_cex = 0.8,
  by.risk = 6,
  risk.add = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  ymin.del = 0.035,
  y.risk1 = NULL,
  y.risk2 = NULL,
  ymin2 = NULL,
  risk_offset = NULL,
  risk.pad = 0.01,
  risk_delta = 0.0275,
  tau_add = NULL,
  time.zero.pad = 0,
  time.zero.label = 0,
  xlabel = NULL,
  ylabel = NULL,
  Maxtau = NULL,
  seedstart = 8316951,
  ylim = NULL,
  draws.band = 1000,
  qtau = 0.025,
  show_resamples = FALSE,
  modify_tau = FALSE,
  weight_variance = c("known", "estimated"),
  weight_model = NULL
)
```

## Arguments

- df:

  Data frame containing survival data.

- tte.name:

  Name of the time-to-event column.

- event.name:

  Name of the event indicator column (0/1; logical and 1/2 coding as in
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) are
  also accepted).

- treat.name:

  Name of the treatment group column (0/1).

- weight.name:

  Optional name of the weights column.

- sg_labels:

  Character vector of subgroup definitions (as logical expressions).

- ltype:

  Line type for curves (default: "s").

- lty:

  Line style for curves (default: 1).

- draws:

  Number of draws for the resampling-based pointwise variance, passed to
  [`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  (default: 0). With the default the pointwise standard error, the
  pointwise limits and the standardisation of the simultaneous band use
  the analytic variance (Greenwood's formula; the robust variance when
  the weights are not all equal). `draws > 0` only has an effect when
  the weights are not all equal: the variance is then the variance over
  that many resampling draws. The default was 20 before version 0.2.0.

- lwd:

  Line width for curves (default: 2).

- sg_colors:

  Colors for subgroup curves.

- color:

  Color for confidence band polygon (default: "lightgrey").

- ymax.pad:

  Padding for y-axis limits.

- ymin.pad:

  Padding for y-axis limits.

- taus:

  Vector for time truncation (default: c(-Inf, Inf)).

- yseq_length:

  Number of y-axis ticks (default: 5).

- cex_Yaxis:

  Text size for axis.

- risk_cex:

  Text size for risk table.

- by.risk:

  Interval for risk table time points (default: 6).

- risk.add:

  Additional time points for risk table.

- xmax:

  Additional graphical and calculation parameters.

- ymin:

  Additional graphical and calculation parameters.

- ymax:

  Additional graphical and calculation parameters.

- ymin.del:

  Additional graphical and calculation parameters.

- y.risk1:

  Additional graphical and calculation parameters.

- y.risk2:

  Additional graphical and calculation parameters.

- ymin2:

  Additional graphical and calculation parameters.

- risk_offset:

  Additional graphical and calculation parameters.

- risk.pad:

  Additional graphical and calculation parameters.

- risk_delta:

  Additional graphical and calculation parameters.

- tau_add:

  Additional graphical and calculation parameters.

- time.zero.pad:

  Additional graphical and calculation parameters.

- time.zero.label:

  Additional graphical and calculation parameters.

- xlabel:

  Additional graphical and calculation parameters.

- ylabel:

  Additional graphical and calculation parameters.

- Maxtau:

  Upper end of the time range. By default the smaller of the two arms'
  last event times; it must be supplied when an arm has no events.

- seedstart:

  Additional graphical and calculation parameters.

- ylim:

  Additional graphical and calculation parameters.

- draws.band:

  Number of draws for the simultaneous confidence band (default: 1000);
  `0` draws the pointwise limits only. The band's critical value is a
  quantile of the suprema of that many resampled curves. The default was
  20 before version 0.2.0, which left the critical value to 20 draws; a
  call that supplies `draws.band` is unaffected.

- qtau:

  Quantile for time range in simultaneous bands (default: 0.025).

- show_resamples:

  Logical; whether to plot resampled curves (default: FALSE).

- modify_tau:

  Logical; restrict time range for bands.

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

(Invisible) list containing KM_diff results, time points, subgroup
curves, risk tables, and confidence intervals.

## Details

The pointwise standard error of the difference, the pointwise limits and
the standardisation of the simultaneous band come from
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
with the `draws` and `draws.band` given here. Since version 0.2.0 the
default is `draws = 0`: the pointwise variance is the analytic one
(Greenwood's formula, or the robust variance when the weights are not
all equal). Before, the default `draws = 20` made the pointwise standard
error of a weighted plot, and with it the band, a 20-draw resampling
variance. Plots without weights, or with all weights equal, are
numerically unaffected by `draws`. Since version 0.2.0 the default of
`draws.band` is 1000 (it was 20), so the critical value of the
simultaneous band of a call that does not pass `draws.band` rests on
1000 resampled curves. Zero-weight rows contribute nothing and do not
affect which variance formula runs.

`weight_variance` and `weight_model` are passed to the overall
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
call (pointwise limits and band). The subgroup curves are point
estimates and are computed as before.

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
[`cumulative_rmst_bands`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
or `plotKM.band_subgroups`.
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
  analysed. (The subgroup curves of `plotKM.band_subgroups()` are point
  estimates and are not affected.)

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
