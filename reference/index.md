# Package index

## Data Preparation

Create counting-process datasets for survival analysis

- [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  : Weighted and Stratified Survival Analysis
- [`get_dfcounting()`](https://larry-leon.github.io/weightedsurv/reference/get_dfcounting.md)
  : Get df_counting results
- [`extract_group_data()`](https://larry-leon.github.io/weightedsurv/reference/extract_group_data.md)
  : Extract time, event, and weight data for a group
- [`get_event_risk_matrices()`](https://larry-leon.github.io/weightedsurv/reference/get_event_risk_matrices.md)
  : Event and Risk Matrices for Survival Analysis
- [`get_censoring_and_events()`](https://larry-leon.github.io/weightedsurv/reference/get_censoring_and_events.md)
  : Get censoring and event times and their indices
- [`calculate_risk_event_counts()`](https://larry-leon.github.io/weightedsurv/reference/calculate_risk_event_counts.md)
  : Calculate risk set and event counts at time points
- [`count_weighted()`](https://larry-leon.github.io/weightedsurv/reference/count_weighted.md)
  : Weighted counting process
- [`risk_weighted()`](https://larry-leon.github.io/weightedsurv/reference/risk_weighted.md)
  : Weighted risk set

## Weighted Cox Models

Fit weighted Cox proportional hazards models with resampling inference

- [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  : Weighted Cox Model with Rho-Gamma Weights
- [`cox_rhogamma_resample()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma_resample.md)
  : Resampling for Weighted Cox Model (rho, gamma)
- [`cox_score_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_score_rhogamma.md)
  : Cox score with rho-gamma weights
- [`find_cox_root()`](https://larry-leon.github.io/weightedsurv/reference/find_cox_root.md)
  : Root-finding for Cox score function
- [`ci_cox()`](https://larry-leon.github.io/weightedsurv/reference/ci_cox.md)
  : Confidence interval for Cox model estimate
- [`score_calculation()`](https://larry-leon.github.io/weightedsurv/reference/score_calculation.md)
  : Score calculation for weighted Cox model
- [`N_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/N_rhogamma.md)
  : Weighted event count for Cox model

## Weighted Log-Rank Tests

Weighted log-rank statistics and diagnostics

- [`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
  : Weighted Log-Rank and Difference Estimate at a Specified Time
- [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  : Weighted log-rank cumulative statistics
- [`check_results()`](https://larry-leon.github.io/weightedsurv/reference/check_results.md)
  : Check and Compare Statistical Test Results

## Kaplan-Meier Estimation

Weighted KM curves, survival differences, and confidence bands

- [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  : Kaplan-Meier Difference Between Groups
- [`KM_estimates()`](https://larry-leon.github.io/weightedsurv/reference/KM_estimates.md)
  : Kaplan-Meier Survival Estimates and Variance
- [`resampling_survival()`](https://larry-leon.github.io/weightedsurv/reference/resampling_survival.md)
  : Resampling Survival Curves for Confidence Bands
- [`km_quantile()`](https://larry-leon.github.io/weightedsurv/reference/km_quantile.md)
  : Kaplan-Meier quantile and confidence interval
- [`km_quantile_table()`](https://larry-leon.github.io/weightedsurv/reference/km_quantile_table.md)
  : Table of KM quantiles for two groups
- [`kmq_calculations()`](https://larry-leon.github.io/weightedsurv/reference/kmq_calculations.md)
  : Kaplan-Meier quantile calculation
- [`check_km_curve()`](https://larry-leon.github.io/weightedsurv/reference/check_km_curve.md)
  : Check KM curve for validity

## KM Plotting

Plot weighted Kaplan-Meier curves with risk tables and annotations

- [`plot_weighted_km()`](https://larry-leon.github.io/weightedsurv/reference/plot_weighted_km.md)
  : Plot weighted Kaplan-Meier curves
- [`plot_km()`](https://larry-leon.github.io/weightedsurv/reference/plot_km.md)
  : Plot Kaplan-Meier curves
- [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)
  : Plot Kaplan-Meier Survival Difference Curves with Subgroups and
  Confidence Bands
- [`plot_km_curves_counting()`](https://larry-leon.github.io/weightedsurv/reference/plot_km_curves_counting.md)
  : Plot KM curves for two groups with optional confidence intervals and
  censoring marks
- [`plot_km_confint_polygon()`](https://larry-leon.github.io/weightedsurv/reference/plot_km_confint_polygon.md)
  : Plot confidence interval polygon for KM curve
- [`KM_plot_2sample_weighted_counting()`](https://larry-leon.github.io/weightedsurv/reference/KM_plot_2sample_weighted_counting.md)
  : Plot Weighted Kaplan-Meier Curves for Two Samples (Counting Process
  Format)
- [`add_risk_table()`](https://larry-leon.github.io/weightedsurv/reference/add_risk_table.md)
  : Add risk table annotation to KM plot
- [`add_legends()`](https://larry-leon.github.io/weightedsurv/reference/add_legends.md)
  : Add legends to KM plot
- [`add_median_annotation()`](https://larry-leon.github.io/weightedsurv/reference/add_median_annotation.md)
  : Add median annotation to KM plot
- [`get_riskpoints()`](https://larry-leon.github.io/weightedsurv/reference/get_riskpoints.md)
  : Get risk set counts at specified risk points

## RMST

Restricted Mean Survival Time calculations with confidence bands

- [`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
  : Cumulative RMST bands for survival curves

## Weighting Schemes

Weight function computation, validation, and visualization

- [`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)
  : Compute Time-Dependent Weights for Survival Analysis
- [`get_weights()`](https://larry-leon.github.io/weightedsurv/reference/get_weights.md)
  : Get weights for a weighting scheme
- [`get_validated_weights()`](https://larry-leon.github.io/weightedsurv/reference/get_validated_weights.md)
  : Get validated weights for a data frame
- [`extract_and_calc_weights()`](https://larry-leon.github.io/weightedsurv/reference/extract_and_calc_weights.md)
  : Extract and calculate weights for multiple schemes
- [`validate_scheme_params()`](https://larry-leon.github.io/weightedsurv/reference/validate_scheme_params.md)
  : Validate weighting scheme parameters
- [`plot_weight_schemes()`](https://larry-leon.github.io/weightedsurv/reference/plot_weight_schemes.md)
  : Plot weight schemes for survival analysis

## Baseline Tables

Publication-ready baseline characteristics tables

- [`create_baseline_table()`](https://larry-leon.github.io/weightedsurv/reference/create_baseline_table.md)
  : Create Baseline Characteristics Table by Treatment Arm

## Utilities

Input validation and formatting helpers

- [`validate_input()`](https://larry-leon.github.io/weightedsurv/reference/validate_input.md)
  : Validate required columns in a data frame
- [`format_pval()`](https://larry-leon.github.io/weightedsurv/reference/format_pval.md)
  : Format p-value for display
- [`safe_run()`](https://larry-leon.github.io/weightedsurv/reference/safe_run.md)
  : Safe execution wrapper
