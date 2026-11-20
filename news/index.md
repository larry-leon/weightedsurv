# Changelog

## weightedsurv 0.2.0

### Changes in behaviour

- Input checks
  ([`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
  [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)):
  missing values in the time, event, treatment, weight or strata column,
  non-finite or negative times, negative weights, a missing column
  (including `strata.name` and `weight.name`) and an empty treatment arm
  are now errors that name the column or arm (and the subgroup in
  [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)).
  Before, some of these dropped into cryptic errors and some gave wrong
  results without a warning (`NA` survival for an `NA` weight, curves
  from negative weights or negative times). Rows are never dropped
  silently. The event column may be 0/1, logical, or 1/2 as
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) reads
  it, in all three functions (before,
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  rejected 1/2). Zero weights are allowed: in
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  those rows no longer cause a false “Discrepancy in KM curve fit”
  warning, and the Cox model and log-rank test are fitted on the
  positive-weight rows (before,
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) stopped on
  zero weights). On degenerate data (an arm without events, no events at
  all, one subject per arm)
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  now returns a result; quantities that cannot be computed (for example
  the hazard ratio when an arm has no events) are `NA`, with one warning
  that says why. In a stratified analysis a stratum in which one arm is
  empty contributes nothing (that arm’s stratum curve is `NA`, with the
  same single warning) instead of stopping the call.
  [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)
  asks for `Maxtau` when an arm has no events.
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  warns when an arm has no events (the estimate is not finite). (B4)

- [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  returns the Kaplan-Meier value at every requested time point, whatever
  `at_points` holds. With `modify_tau = TRUE`, with a user grid that
  omits event times, and through
  `plotKM.band_subgroups(modify_tau = TRUE)`, the curves, differences,
  variances, resampled processes and bands now equal the Kaplan-Meier
  estimate at those points; before they could be far off (errors above
  0.4 and negative survival were seen). `at_points` must now be numeric,
  finite and non-decreasing (unsorted points used to give values outside
  \[0, 1\]). With the default grid, and in
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
  nothing changes. (B1)

- Quantiles
  ([`km_quantile()`](https://larry-leon.github.io/weightedsurv/reference/km_quantile.md),
  [`kmq_calculations()`](https://larry-leon.github.io/weightedsurv/reference/kmq_calculations.md),
  `df_counting()$quantile_results`) now equal
  `quantile(survfit(...), probs = 1 - qprob)`: when the curve equals the
  level on an interval the midpoint is returned (before, the left end;
  for example a median of 5.0 is now 5.5), a plateau at the level that
  runs to the end of follow-up is handled as
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) does, and
  the confidence limits are read the same way when the limit curve is
  not monotone. `qprob` keeps its meaning, the time at which survival
  drops to `qprob`, and is documented. `df_counting()$quantile_check`
  now asks [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
  for the same quantity (`probs = 1 - qprob`, at `conf_level`); a real
  discrepancy is a warning even when `verbose = FALSE`. (B3)

- [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md):
  the simultaneous-band critical value is the `1 - alpha` quantile (it
  was always the 0.95 quantile); results at `alpha = 0.05` are
  unchanged.
  [`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
  uses the pooled survival at `tzero` (it used the next grid point), so
  `sig2_dhat` and `cov_wlr_dhat` change slightly (by 0-8% on the test
  datasets; the correlation is unchanged); `tzero` before the first time
  point gives a difference of 0, and after the last time point the
  values there (before, an empty result). Magirr-Burman weights with
  `mb_tstar` at or beyond the last time point now use the survival at
  the last time point (before, all weights were 1, or infinite where
  survival was 0). (B5)

- The Kaplan-Meier plots
  ([`plot_weighted_km()`](https://larry-leon.github.io/weightedsurv/reference/plot_weighted_km.md),
  [`KM_plot_2sample_weighted_counting()`](https://larry-leon.github.io/weightedsurv/reference/KM_plot_2sample_weighted_counting.md),
  [`plot_km_curves_counting()`](https://larry-leon.github.io/weightedsurv/reference/plot_km_curves_counting.md))
  now draw each curve from time zero at S = 1, as a step function
  identical to
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)’s, and
  the confidence band as a step band that stays at or below 1 and closes
  onto the curve where survival reaches 0. Before, each curve started at
  its first event (the first drop and the line under early censoring
  marks were missing) and the band joined the event times with straight
  lines and could rise above 1. Everything else on the plot is
  unchanged. (B10)

- Functions that draw random numbers
  ([`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  with `draws` or `draws.band`,
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  and
  [`calculate_risk_event_counts()`](https://larry-leon.github.io/weightedsurv/reference/calculate_risk_event_counts.md)
  with `draws`,
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  with `draws`) no longer reset the random-number stream of the calling
  session; results for a given `seedstart` are unchanged. (B6)

- Weighted standard errors: when the weights are not all equal, the
  standard error of a weighted Kaplan-Meier curve is now the robust
  (infinitesimal-jackknife) standard error,
  `S(t)^2 * sum_i w_i^2 psi_i(t)^2` with `psi_i` subject i’s influence
  contribution, in
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  and
  [`calculate_risk_event_counts()`](https://larry-leon.github.io/weightedsurv/reference/calculate_risk_event_counts.md).
  It equals the standard error of
  `survival::survfit(weights = w, robust = TRUE)`. **Every weighted
  standard error changes**, and with it everything computed from one:
  pointwise limits, simultaneous-band limits and critical values, the
  standardised resampled curves, and the confidence limits of weighted
  quantiles. **Unweighted results do not change**, and neither do
  results with all weights equal (both use Greenwood’s formula; all
  weights equal to a constant give the unweighted standard error, where
  in 0.1.0 a constant weight c gave the standard error of the data with
  every row repeated c times). Survival estimates, differences, risk
  tables and the RMST point estimate do not change (the log-rank test
  and Cox estimate with subject weights change separately, see L1-L4
  below). The standard error does not depend on the scale of the weights
  (weights `w` and `10 * w` give the same answer), so weights are not
  frequency weights: a weight of 2 is not two copies of the row. By
  default weights are treated as fixed and known; for weights estimated
  by a logistic propensity model see `weight_variance = "estimated"`
  below. With extreme weights (low effective sample size) the standard
  error can be too small. (M, R1)

- Weighted resampling: when the weights are not all equal, resampled
  Kaplan-Meier curves are now built from the same per-subject influence
  contributions, `-S(t) * sum_i G_i w_i psi_i(t)` for the same
  standard-normal multipliers `G_i` as before (before, the multipliers
  perturbed the events only). This applies to
  [`resampling_survival()`](https://larry-leon.github.io/weightedsurv/reference/resampling_survival.md),
  the resampled curves of
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  (`surv0_star`, `surv1_star`, `dhat_star`, `Zdhat_star`), and the
  `draws > 0` variance of
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  and
  [`calculate_risk_event_counts()`](https://larry-leon.github.io/weightedsurv/reference/calculate_risk_event_counts.md).
  **Weighted simultaneous bands and weighted RMST bands change**
  ([`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md):
  the RMST variance, pointwise limits, critical value and band limits;
  not the RMST estimate), and the resampling variance now agrees with
  the robust standard error. **Unweighted results, and results with all
  weights equal, do not change**: they keep the events-only
  perturbation, with the same seeds and draws. (R2)

- Log-rank test with subject weights: when the subject weights are not
  all equal, the variance of the weighted log-rank score is now the sum
  of squared per-subject score contributions, `sum_i psi_i(0)^2` (the
  formula is in
  [`?df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)),
  in
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  (unstratified and per stratum) and in
  [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  when it is given a
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  result. For Fleming-Harrington (0, 0) weights the z-score squared is
  the robust score test of
  `survival::coxph(weights = w, robust = TRUE, ties = "breslow")`.
  **z-scores and p-values of subject-weighted log-rank tests change**
  (`z.score`, `sig2_lr`, `get_score`, `zlogrank_text`,
  `z.score_stratified`, `sig2_lr_stratified`; the cumulative
  `sig2.score` and `z.score` of
  [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md));
  the score itself (`lr`, `score`, `lr_stratified`) does not. They no
  longer depend on how the weights are scaled: `w`, `10 * w` and
  `w / sum(w)` give the same z-score (before, multiplying the weights by
  10 multiplied z by about 3.2, and normalised weights gave `NaN`).
  **Unweighted results, and results with all weights equal, do not
  change.**
  [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  on a bare table of counts (no subject-level data) is unchanged too. By
  default weights are treated as fixed and known; for weights estimated
  by a logistic propensity model see `weight_variance = "estimated"`
  below. With extreme weights (low effective sample size) the test can
  reject too often. (L1)

- [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  with subject weights: when the weights are not all equal, the standard
  error of the log hazard ratio is now
  `sqrt(sum_i psi_i(bhat)^2) / i_bhat`, from the same per-subject score
  contributions; for Fleming-Harrington (0, 0) weights it is the robust
  standard error of
  `survival::coxph(weights = w, robust = TRUE, ties = "breslow")`.
  **Standard errors, confidence intervals, z-scores and p-values of
  subject-weighted fits change** (`fit$sig_bhat_asy`, `hr_ci_asy`,
  `cox_text_asy`, `fit$sig2_score`, `fit$z.score`, `zlogrank_text`, and
  `sig2_score` and `sig2_beta_asy` of
  `cox_score_rhogamma(score_only = FALSE)`), and no longer depend on how
  the weights are scaled. The estimate (`bhat`), the scores (`u.beta`,
  `u.zero`), the information and the time weights do not change.
  **Unweighted fits, and fits with all weights equal, do not change.**
  (L2)

- `cox_rhogamma(draws > 0)` with subject weights: when the weights are
  not all equal, the resampled scores of
  [`cox_rhogamma_resample()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma_resample.md)
  are now `-sum_i G_i psi_i(bhat)` and `-sum_i G_i psi_i(0)`, built from
  the same per-subject contributions with the same multipliers `G_i` as
  before (before, the multipliers perturbed the events only). **The
  resampling standard error, the de-biased estimate and the resampling
  interval of subject-weighted fits change** (`fit_resamples`,
  `fit$sig_bhat_star`, `fit$bhat_debiased`, `fit$wald_debiased1`,
  `fit$wald_debiased2`, `hr_ci_star`, `sig2_score_star`,
  `z.score_debiased`, `cox_text_star`); the resampling standard error
  now agrees with `sig_bhat_asy`. **Unweighted fits, and fits with all
  weights equal, do not change**, with the same seeds and draws. (L3)

- [`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
  with subject weights: on a
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  result whose weights are not all equal, the variance of the log-rank
  statistic (as L1), the variance of the survival difference at `tzero`
  and their covariance are now computed from the same per-subject
  contributions; `sig2_dhat` equals `KM_diff(draws = 0)$sig2_dhat` at
  that time point. **`sig2.score`, `z.score`, `sig2_dhat`,
  `cov_wlr_dhat` and `cor_wlr_dhat` of subject-weighted calls change**;
  `score` and `dhat` do not. **Unweighted results, and results with all
  weights equal, do not change.** (L4)

- [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md):
  the default of `draws` is now 0 (it was 20). With weights that are not
  all equal and `draws` not supplied, the pointwise standard error, the
  pointwise limits and the standardisation of the simultaneous band now
  use the analytic (robust) variance instead of a variance from 20
  resampling draws, so **weighted band plots that do not pass `draws`
  change**. Plots without weights, or with all weights equal, are
  numerically unaffected. Pass `draws = 20` for the earlier behaviour.
  (B12)

- [`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
  has a new last argument `mean_rule`. The default, `"trapezoid"`, gives
  the same results as before. `"rectangle"` integrates the survival
  difference as a step function; it then equals the difference of
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)’s
  restricted means whenever the time grid holds every event time (as the
  default
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  grid does; not with `modify_tau = TRUE`). The result gains an element
  `mean_rule`. (B8)

- [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  no longer plots the resampled difference curves by default
  (`show_resamples = FALSE`); use `show_resamples = TRUE` to see them.
  Results are unchanged. (B9)

- The two confidence bands of the Kaplan-Meier plots
  ([`plot_weighted_km()`](https://larry-leon.github.io/weightedsurv/reference/plot_weighted_km.md),
  [`KM_plot_2sample_weighted_counting()`](https://larry-leon.github.io/weightedsurv/reference/KM_plot_2sample_weighted_counting.md),
  [`plot_km_curves_counting()`](https://larry-leon.github.io/weightedsurv/reference/plot_km_curves_counting.md))
  are now semi-transparent, in the same colours, so each stays visible
  where they overlap. Nothing else on the plot changes. (B11)

- New option for estimated weights.
  [`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
  (new export) fits a logistic propensity model, or takes a fitted one,
  and returns inverse-probability or stabilised weights together with
  the influence of the model’s coefficients and the derivative of each
  weight. Passed as `weight_model` with `weight_variance = "estimated"`
  to
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
  [`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
  or
  [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)
  (two new last arguments), it carries the estimation of the weights
  into the Kaplan-Meier standard errors, the variance of the survival
  difference (no longer the sum of the arm variances), the simultaneous
  band, the RMST variance and band, the log-rank z-score (stratified
  too), `logrank_results`, `cox_results$score`, `$se`, `$hr_ci` and
  `$pval`;
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
  [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  and
  [`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
  take the option from the
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  result they are given. **The default is `weight_variance = "known"`: a
  call that does not ask for `"estimated"` returns what it returned
  before.** No estimate changes under either option. `"estimated"` needs
  a logistic `glm` fitted on exactly the rows analysed (subgroup
  analyses that reuse full-sample weights must use `"known"`) and
  `draws = 0`; see the section “Estimated weights” of
  [`?df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md).
  Results gain an element `weight_variance`. (W1-W4)

- Zero weights: whether the weights are “all equal” is now decided on
  the positive weights. **Calls whose weight column contains a zero
  change**: a 0/1 weight column now gives the variances, z-scores and
  intervals of the positive-weight rows analysed without weights
  (before, the zeros made the weights count as unequal and the robust
  formulas ran; for example a log-rank z of -0.7569 against -0.7554 for
  the same rows without weights). With zeros among unequal positive
  weights the results are those of the positive-weight rows with the
  same weights, as before. Resampled quantities keep a multiplier per
  row and are not equal to the subset’s draw by draw. Calls without a
  zero weight do not change. (Rd)

- [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  with subject weights: `logrank_results` and `cox_results$score` are
  now robust tests of the weighted analysis. **Calls whose positive
  weights are not all equal change in these two outputs only**:
  `logrank_results` is the robust score test at 0 of the weighted Cox
  model (before, an unweighted
  [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) test
  that ignored the weights), and `cox_results$score` is
  `summary(cox_fit)$robscore` (before, the model-based score test, which
  reads the weights as counts). With subject weights `logrank_results`
  is the `rho = 0` test whatever the scheme.
  [`check_results()`](https://larry-leon.github.io/weightedsurv/reference/check_results.md)
  now agrees on weighted Fleming-Harrington (0, 0) results instead of
  warning. New elements `logrank_results$method` (`"survdiff"` or
  `"coxph_robust_score"`) and `cox_results$se` (the standard error of
  the log hazard ratio behind `hr_ci` and `pval`). Without weights, and
  with all positive weights equal, both outputs are unchanged. (Rb)

- [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md):
  the default of `draws.band` is now 1000 (it was 20), so the critical
  value of the simultaneous band no longer rests on 20 draws. **Calls
  that do not pass `draws.band` change in the band (critical value, band
  limits, resampled curves)**; pointwise limits and estimates do not.
  Calls that pass `draws.band` are unchanged; pass `draws.band = 20` for
  the earlier behaviour. (Rc)

- [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  returns a new element `ess`: by arm, the number of rows with positive
  weight, the effective sample size `(sum w)^2 / sum(w^2)` and their
  ratio. A low ratio means a few subjects carry most of the weight;
  standard errors, intervals and tests can then be too optimistic. No
  existing element changes and no warning depends on it. (Ra)

- [`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
  has the arguments `formula` and `data` (the logistic propensity model
  is fitted by [`glm()`](https://rdrr.io/r/stats/glm.html)), `model` (a
  fitted logistic `glm`, in place of `formula` and `data`) and
  `type = c("ipw", "stabilized")`. The result has class `"ps_weights"`
  with a [`print()`](https://rdrr.io/r/base/print.html) method
  ([`print.ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md))
  that shows, by arm, the number of subjects, the smallest and largest
  weight, the effective sample size and its ratio to the number of
  subjects.

### Bug fixes

- Stratified
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  no longer depends on `survival` being attached or recent enough to
  find [`strata()`](https://rdrr.io/pkg/survival/man/strata.html). (B7)

- Weighted Kaplan-Meier curves that end on an event now end at
  exactly 0. Before, a residue of about 1e-17, sometimes negative, made
  [`plot_weighted_km()`](https://larry-leon.github.io/weightedsurv/reference/plot_weighted_km.md)
  stop with “values outside \[0,1\]” and gave “NaNs produced” warnings.
  Confidence limits are `NA` where survival is 0 (as in
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)), and the
  plotted band no longer receives `NaN`. (B2)

- [`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
  and
  [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  with `scheme = "schemper"` now check for the censoring curve and say
  what is missing;
  [`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
  says that it needs a
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  fit made with `draws.band > 0`;
  `KM_diff(draws.band < 20, show_resamples = TRUE)` no longer fails;
  [`get_validated_weights()`](https://larry-leon.github.io/weightedsurv/reference/get_validated_weights.md)
  no longer prints to the console. (B5, B4, B6)

- Internal tidying with no change in results: a local function inside
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  no longer shadows the exported
  [`check_km_curve()`](https://larry-leon.github.io/weightedsurv/reference/check_km_curve.md),
  and unused code was removed. (V)

### Documentation

- The vignettes and articles were revised for 0.2.0: the Rotterdam
  example now uses postmenopausal women without chemotherapy and
  [`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
  with `weight_variance = "estimated"`; effective sample size tables
  were added; a new article gives the statistical details of estimated
  weights.
- Vignette chunks that use a suggested package run only when it is
  installed; website-only packages are listed in `Config/Needs/website`.
