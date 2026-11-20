# Weighted and Stratified Survival Analysis

Performs comprehensive weighted and/or stratified survival analysis,
including Cox proportional hazards model, logrank/Fleming-Harrington
tests, and calculation of risk/event sets, Kaplan-Meier curves,
quantiles, and variance estimates.

## Usage

``` r
df_counting(
  df,
  tte.name = "tte",
  event.name = "event",
  treat.name = "treat",
  weight.name = NULL,
  strata.name = NULL,
  arms = c("treat", "control"),
  time.zero = 0,
  tpoints.add = c(0),
  by.risk = 6,
  time.zero.label = 0,
  risk.add = NULL,
  get.cox = TRUE,
  cox.digits = 2,
  lr.digits = 2,
  cox.eps = 0.001,
  lr.eps = 0.001,
  verbose = FALSE,
  qprob = 0.5,
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 0),
  conf_level = 0.95,
  check.KM = TRUE,
  check.seKM = FALSE,
  draws = 0,
  seedstart = 8316951,
  stop.onerror = FALSE,
  censoring_allmarks = TRUE,
  weight_variance = c("known", "estimated"),
  weight_model = NULL
)
```

## Arguments

- df:

  Data frame containing survival data.

- tte.name:

  Character; name of the time-to-event column in `df`.

- event.name:

  Character; name of the event indicator column in `df` (1=event,
  0=censored). Logical (TRUE = event) and 1/2 coding (2 = event, as
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) reads
  it) are also accepted.

- treat.name:

  Character; name of the treatment/group column in `df` (must be coded
  as 0=control, 1=experimental).

- weight.name:

  Character or NULL; name of the weights column in `df`. If NULL, equal
  weights are used.

- strata.name:

  Character or NULL; name of the strata column in `df` for stratified
  analysis.

- arms:

  Character vector of length 2; group labels. Default:
  `c("treat","control")`.

- time.zero:

  Numeric; time value to use as zero. Default: 0.

- tpoints.add:

  Numeric vector; additional time points to include in calculations.
  Default: `c(0)`.

- by.risk:

  Numeric; interval for risk set time points. Default: 6.

- time.zero.label:

  Numeric; label for time zero in output. Default: 0.0.

- risk.add:

  Numeric vector or NULL; additional specific risk points to include.

- get.cox:

  Logical; whether to fit Cox proportional hazards model. Default: TRUE.

- cox.digits:

  Integer; number of decimal places for Cox output formatting. Default:
  2.

- lr.digits:

  Integer; number of decimal places for logrank output formatting.
  Default: 2.

- cox.eps:

  Numeric; threshold for Cox p-value formatting (values below shown as
  "\<eps"). Default: 0.001.

- lr.eps:

  Numeric; threshold for logrank p-value formatting. Default: 0.001.

- verbose:

  Logical; whether to print warnings and diagnostic messages. Default:
  FALSE.

- qprob:

  Numeric in (0,1); survival level for the KM quantile table: the
  quantile reported is the time at which survival drops to `qprob`, i.e.
  `quantile(survfit(...), probs = 1 - qprob)`. Default: 0.5 (median).

- scheme:

  Character; weighting scheme for logrank/Fleming-Harrington test.
  Options: "fh", "schemper", "XO", "MB", "custom_time", "fh_exp1",
  "fh_exp2". Default: "fh".

- scheme_params:

  List; parameters for the selected weighting scheme. Default:
  `list(rho = 0, gamma = 0)`.

  - For "fh": `rho` and `gamma` (Fleming-Harrington parameters)

  - For "MB": `mb_tstar` (cutoff time)

  - For "custom_time": `t.tau`, `w0.tau`, `w1.tau`

- conf_level:

  Numeric in (0,1); confidence level for quantile intervals. Default:
  0.95.

- check.KM:

  Logical; whether to check KM curve validity against
  [`survival::survfit`](https://rdrr.io/pkg/survival/man/survfit.html).
  Default: TRUE.

- check.seKM:

  Logical; whether to check KM standard error estimates. Default: FALSE.

- draws:

  Integer; number of draws for resampling-based variance estimation.
  Default: 0 (no resampling).

- seedstart:

  Integer; random seed for reproducible resampling. Default: 8316951.

- stop.onerror:

  Logical; whether to stop execution on errors (TRUE) or issue warnings
  (FALSE). Default: FALSE.

- censoring_allmarks:

  Logical; if FALSE, removes events from censored time points. Default:
  TRUE.

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

A list with the following components:

- cox_results:

  List with the results of
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html):
  `cox_fit`, `hr`, `hr_ci`, `pval`, `score`, `cox_text` and `se`, the
  standard error of the log hazard ratio behind `hr_ci` and `pval` (the
  robust one of `coxph(robust = TRUE)`; `NA` when the model could not be
  fitted). `score` is the score test of the model: with no weights, or
  with all positive weights equal, `summary(cox_fit)$sctest`; with
  subject weights that are not all equal, the robust score test
  `summary(cox_fit)$robscore` (the model-based one reads the weights as
  counts).

- logrank_results:

  List with `chisq`, `pval`, `logrank_text` and `method`. With no
  weights, or with all positive weights equal, the test of
  [`survival::survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html)
  (`method = "survdiff"`; `rho` from the scheme when it is `"fh"`). With
  subject weights that are not all equal
  [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) has no
  counterpart, and the result is the robust score test at 0 of the
  weighted Cox model
  (`coxph(weights = , init = 0, iter.max = 0, robust = TRUE, ties = "breslow")`
  on the positive-weight rows, strata included;
  `method = "coxph_robust_score"`): the `rho = 0` test whatever the
  scheme. Its `chisq` equals `z.score^2` (`z.score_stratified^2`) for
  Fleming-Harrington (0, 0) weights, so
  [`check_results`](https://larry-leon.github.io/weightedsurv/reference/check_results.md)
  agrees on weighted data.

- z.score:

  Standardized weighted log-rank test statistic (with subject weights
  that are not all equal, standardized by the score-residual variance;
  see the section on subject weights)

- at_points:

  Vector of all time points used in calculations

- surv0, surv1:

  Kaplan-Meier survival estimates for control and treatment groups

- sig2_surv0, sig2_surv1:

  Variance estimates for survival curves

- survP:

  Pooled survival estimates

- survG:

  Censoring distribution estimates

- quantile_results:

  Data frame with median survival and confidence intervals by group

- lr, sig2_lr:

  Weighted log-rank statistic and its variance (`get_score` holds the
  same as returned by
  [`wlr_cumulative`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md);
  `zlogrank_text` is the one-sided p-value of `z.score`)

- riskpoints0, riskpoints1:

  Risk set counts at specified time points

- z.score_stratified:

  Stratified z-score (if stratified analysis), with `lr_stratified` and
  `sig2_lr_stratified`

- ess:

  Data frame with one row per arm, in the order of `quantile_results`
  (`arms[1]`, the treated arm, then `arms[2]`): `group`; `n`, the rows
  of the arm with positive weight; `ess`, the effective sample size
  \\(\sum w)^2 / \sum w^2\\ over those rows; and `ratio = ess / n`.
  Without weights, and with constant weights, `ess = n`. A ratio near 1
  means the weights are even; a low ratio means that a few subjects
  carry most of the weight, and standard errors, intervals and tests can
  then be too optimistic, whether the weights are treated as known or as
  estimated. No warning or message depends on it.

- weight_variance:

  `"known"` or `"estimated"`: how the weights were treated.

- ew_parts:

  Only with `weight_variance = "estimated"`: list with `IF` and `dw` of
  `weight_model`, rows in the order of `time`; read by
  [`cox_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
  [`wlr_cumulative`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  and
  [`wlr_dhat_estimates`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md).

## Details

This function implements a comprehensive survival analysis framework
supporting:

- Weighted observations via `weight.name`

- Stratified analysis via `strata.name`

- Multiple weighting schemes for log-rank tests

- Resampling-based variance estimation

- Automatic validation against `survival` package results

The function performs time-fixing using
[`survival::aeqSurv`](https://rdrr.io/pkg/survival/man/aeqSurv.html) to
handle tied event times.

Input checks: missing values (in any named column, including
`strata.name`), non-finite or negative times, negative weights and an
empty arm are errors; rows are never dropped silently. Zero weights are
allowed: those rows contribute nothing, and the Cox model and log-rank
test are fitted on the positive-weight rows. Zero-weight rows do not
affect which variance formula runs either: whether the weights are "all
equal" is decided on the positive weights, so a 0/1 weight column gives
the results of the positive-weight rows analysed without weights. A
stratum in which one arm is empty contributes nothing to the stratified
statistics (that arm's stratum curve is `NA`). On degenerate data (an
arm without events, no events at all, one subject per arm) a result is
still returned; quantities that cannot be computed are `NA`, with one
warning that says why. For stratified analyses, stratum-specific
estimates are computed and combined appropriately.

## Variance of the Kaplan-Meier estimate

With no weights, or with all weights equal, the variance of a
Kaplan-Meier estimate is Greenwood's formula, and these results are
unchanged. When the weights are not all equal it is the robust
(infinitesimal-jackknife) variance \$\$\widehat{Var}\[\hat S(t)\] = \hat
S(t)^2 \sum_i w_i^2 \psi_i(t)^2,\$\$ \$\$\psi_i(t) = \frac{D_i 1(U_i \le
t)}{Y_w(U_i) - dN_w(U_i)} - \sum\_{s \le \min(t, U_i)}
\frac{dN_w(s)}{Y_w(s)\\Y_w(s) - dN_w(s)\\},\$\$ where \\U_i\\, \\D_i\\
and \\w_i\\ are subject \\i\\'s time, event indicator and weight,
\\Y_w(s)\\ is the weighted number at risk, \\dN_w(s)\\ the weighted
number of events at \\s\\, the sum runs over event times, and a term is
0 where \\Y_w(s) - dN_w(s) \le 10^{-10} Y_w(s)\\. It equals the standard
error of `survival::survfit(weights = w, robust = TRUE)` (on the
time-fixed data), reduces to Greenwood's formula at equal weights, and
does not depend on the scale of the weights. With `draws > 0` and
unequal weights the variance is instead the variance over resampling
draws of \\\sum_i G_i w_i \psi_i(t)\\ (standard-normal multipliers
\\G_i\\), which estimates the same quantity; resampled curves for
weighted data are built from the same per-subject contributions (see
[`resampling_survival`](https://larry-leon.github.io/weightedsurv/reference/resampling_survival.md)).
By default the weights are treated as fixed and known. For weights
estimated by a logistic propensity model,
`weight_variance = "estimated"` carries their estimation into the
standard error (section "Estimated weights" of `df_counting`). When the
weights are extreme (low effective sample size) the standard error can
be too small under either option.

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
a_j\\ dw_j\\; in
[`wlr_dhat_estimates`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
the contribution to the survival difference is \\d_i + IF_i' (\hat S_0
D_0 - \hat S_1 D_1)\\ at the grid point read for `tzero`, and the
covariance uses both adjusted contributions. The outputs of
`df_counting()` that come from the `survival` package
(`logrank_results`, `cox_results$score`, `$se`, `$hr_ci`, `$pval`) carry
the same adjustment, applied to
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)'s weighted
score and dfbeta residuals \\r\\: \\r_i + IF_i' \sum_j (r_j / w_j)
dw_j\\.

How to use it. Fit the propensity model on exactly the rows to be
analysed with
[`ps_weights`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md),
put its weights in a column of `df`, and pass
`weight_variance = "estimated"` and `weight_model` to `df_counting()`,
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`cumulative_rmst_bands`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
or
[`plotKM.band_subgroups`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md).
[`cox_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`wlr_cumulative`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
and
[`wlr_dhat_estimates`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
take the option from the `df_counting()` result they are given (its
elements `weight_variance` and `ew_parts`). The rows of `df`, the
treatment column and the weight column must be those of `weight_model`;
anything else is an error.

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

## Weighting Schemes

- fh:

  Fleming-Harrington: w(t) = S(t)^rho \* (1-S(t))^gamma

- MB:

  Magirr-Burman: w(t) = 1/max(S(t), S(t\*))

- schemper:

  Schemper: w(t) = S(t)/G(t) where G is the censoring distribution

- XO:

  Xu-O'Quigley: w(t) = S(t)/Y(t) where Y is risk set size

## References

Fleming, T. R. and Harrington, D. P. (1991). Counting Processes and
Survival Analysis. Wiley.

Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests.
Statistics in Medicine, 38(20), 3782-3790.

Leon, L. F., Cai, T. and Wei, L. J. (2009). Robust inferences for
covariate effects on survival time with censored linear regression
models. Statistics in Biosciences, 1, 50-64.
[doi:10.1007/s12561-009-9002-5](https://doi.org/10.1007/s12561-009-9002-5)

## See also

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html),
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html)
[`cox_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
for weighted Cox models
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
for Kaplan-Meier differences

Other survival_analysis:
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md),
[`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)

Other weighted_tests:
[`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)

## Examples

``` r
# Basic survival analysis
library(survival)
str(veteran)
#> 'data.frame':    137 obs. of  8 variables:
#>  $ trt     : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ celltype: Factor w/ 4 levels "squamous","smallcell",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ time    : num  72 411 228 126 118 10 82 110 314 100 ...
#>  $ status  : num  1 1 1 1 1 1 1 1 1 0 ...
#>  $ karno   : num  60 70 60 60 70 20 40 80 50 70 ...
#>  $ diagtime: num  7 5 3 9 11 5 10 29 18 6 ...
#>  $ age     : num  69 64 38 63 65 49 69 68 43 70 ...
#>  $ prior   : num  0 10 0 10 10 0 10 0 0 0 ...
veteran$treat <- as.numeric(veteran$trt) - 1

result <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  arms = c("Treatment", "Control")
)

# Print results
print(result$cox_results$cox_text)
#> [1] "HR = 1.02 (0.72, 1.44)"
print(result$zlogrank_text)
#> [1] "logrank (1-sided) p = 0.54"

# Fleming-Harrington (0,1) weights (emphasizing late differences)
result_fh <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 1)
)

# Stratified analysis
result_strat <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  strata.name = "celltype"
)
```
