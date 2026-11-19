# weightedsurv 0.2.0

## Changes in behaviour

* Input checks (`KM_diff()`, `df_counting()`, `plotKM.band_subgroups()`): missing values in the
  time, event, treatment, weight or strata column, non-finite or negative times, negative weights,
  a missing column (including `strata.name` and `weight.name`) and an empty treatment arm are now
  errors that name the column or arm (and the subgroup in `plotKM.band_subgroups()`). Before, some
  of these dropped into cryptic errors and some gave wrong results without a warning (`NA` survival
  for an `NA` weight, curves from negative weights or negative times). Rows are never dropped
  silently. The event column may be 0/1, logical, or 1/2 as `survival::Surv()` reads it, in all
  three functions (before, `KM_diff()` rejected 1/2). Zero weights are allowed: in `df_counting()`
  those rows no longer cause a false "Discrepancy in KM curve fit" warning, and the Cox model and
  log-rank test are fitted on the positive-weight rows (before, `coxph()` stopped on zero weights).
  On degenerate data (an arm without events, no events at all, one subject per arm)
  `df_counting()` now returns a result; quantities that cannot be computed (for example the hazard
  ratio when an arm has no events) are `NA`, with one warning that says why. In a stratified
  analysis a stratum in which one arm is empty contributes nothing instead of stopping the call.
  `plotKM.band_subgroups()` asks for `Maxtau` when an arm has no events. `cox_rhogamma()` warns
  when an arm has no events (the estimate is not finite). (B4)

* `KM_diff()` returns the Kaplan-Meier value at every requested time point, whatever `at_points`
  holds. With `modify_tau = TRUE`, with a user grid that omits event times, and through
  `plotKM.band_subgroups(modify_tau = TRUE)`, the curves, differences, variances, resampled
  processes and bands now equal the Kaplan-Meier estimate at those points; before they could be
  far off (errors above 0.4 and negative survival were seen). `at_points` must now be numeric,
  finite and non-decreasing (unsorted points used to give values outside [0, 1]). With the default
  grid, and in `df_counting()`, nothing changes. (B1)

* Quantiles (`km_quantile()`, `kmq_calculations()`, `df_counting()$quantile_results`) now equal
  `quantile(survfit(...), probs = 1 - qprob)`: when the curve equals the level on an interval the
  midpoint is returned (before, the left end; for example a median of 5.0 is now 5.5), a plateau at
  the level that runs to the end of follow-up is handled as `survfit()` does, and the confidence
  limits are read the same way when the limit curve is not monotone. `qprob` keeps its meaning, the
  time at which survival drops to `qprob`, and is documented. `df_counting()$quantile_check` now
  asks `survfit()` for the same quantity (`probs = 1 - qprob`, at `conf_level`); a real discrepancy
  is a warning even when `verbose = FALSE`. (B3)

* `KM_diff()`: the simultaneous-band critical value is the `1 - alpha` quantile (it was always the
  0.95 quantile); results at `alpha = 0.05` are unchanged. `wlr_dhat_estimates()` uses the pooled
  survival at `tzero` (it used the next grid point), so `sig2_dhat` and `cov_wlr_dhat` change
  slightly (by 0-8% on the test datasets; the correlation is unchanged); `tzero` before the first
  time point gives a difference of 0, and after the last time point the values there (before, an
  empty result). Magirr-Burman weights with `mb_tstar` at or beyond the last time point now use
  the survival at the last time point (before, all weights were 1, or infinite where survival was
  0). (B5)

* The Kaplan-Meier plots (`plot_weighted_km()`, `KM_plot_2sample_weighted_counting()`,
  `plot_km_curves_counting()`) now draw each curve from time zero at S = 1, as a step function
  identical to `survfit()`'s, and the confidence band as a step band that stays at or below 1 and
  closes onto the curve where survival reaches 0. Before, each curve started at its first event
  (the first drop and the line under early censoring marks were missing) and the band joined the
  event times with straight lines and could rise above 1. Everything else on the plot is
  unchanged. (B10)

* Functions that draw random numbers (`KM_diff()` with `draws` or `draws.band`, `df_counting()`
  and `calculate_risk_event_counts()` with `draws`, `cox_rhogamma()` with `draws`) no longer reset
  the random-number stream of the calling session; results for a given `seedstart` are unchanged.
  (B6)

## Bug fixes

* Stratified `df_counting()` no longer depends on `survival` being attached or recent enough to
  find `strata()`. (B7)

* Weighted Kaplan-Meier curves that end on an event now end at exactly 0. Before, a residue of
  about 1e-17, sometimes negative, made `plot_weighted_km()` stop with "values outside [0,1]" and
  gave "NaNs produced" warnings. Confidence limits are `NA` where survival is 0 (as in `survfit()`),
  and the plotted band no longer receives `NaN`. (B2)

* `wlr_dhat_estimates()` and `wlr_cumulative()` with `scheme = "schemper"` now check for the
  censoring curve and say what is missing; `cumulative_rmst_bands()` says that it needs a
  `KM_diff()` fit made with `draws.band > 0`; `KM_diff(draws.band < 20, show_resamples = TRUE)` no
  longer fails; `get_validated_weights()` no longer prints to the console. (B5, B4, B6)

* Internal tidying with no change in results: a local function inside `df_counting()` no longer
  shadows the exported `check_km_curve()`, and unused code was removed. (V)
