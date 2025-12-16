# weightedsurv

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://cran.r-project.org/web/checks/check_results_weightedsurv.html)
[![CRAN status](https://www.r-pkg.org/badges/version/weightedsurv)](https://CRAN.R-project.org/package=weightedsurv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview
  
**weightedsurv** provides tools for weighted and stratified survival analysis, including Cox proportional hazards models, weighted log-rank tests, Kaplan-Meier curves, and Restricted Mean Survival Time (RMST) calculations with flexible weighting schemes.

The package is particularly useful for clinical trial analyses where treatment effects may vary over time, such as in oncology trials with delayed treatment effects or crossing hazards.

## Installation

Install the released version from CRAN:

```r
install.packages("weightedsurv")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("larry-leon/weightedsurv")
```

## Key Features

- **Flexible weighting schemes** for log-rank tests and Cox models
- **Stratified analysis** support for adjusting confounders
- **Subject-specific weights** for complex sampling designs
- **Resampling-based inference** for improved small-sample properties
- **Simultaneous confidence bands** for survival curves and differences
- **RMST analysis** with confidence intervals
- **Publication-ready tables** for baseline characteristics

## Weighting Schemes

The package supports multiple weighting schemes for weighted log-rank tests:
 
| Scheme | Description | Use Case |
|--------|-------------|----------|
| Fleming-Harrington (fh) | w(t) = S(t)^ρ × (1-S(t))^γ | Emphasize early (ρ>0) or late (γ>0) differences |
| Magirr-Burman (MB) | w(t) = 1/max(S(t), S(t*)) | Modest downweighting after cutoff time |
| Schemper | w(t) = S(t)/G(t) | Adjust for informative censoring |
| Xu-O'Quigley (XO) | w(t) = S(t)/Y(t) | Adjust for risk set size |
| Custom time-based | Step function weights | Pre-specified time-based weighting |

## Usage

### Basic Survival Analysis

```r
library(weightedsurv)
library(survival)

# Prepare data (treatment must be coded 0/1)
data(veteran)
veteran$treat <- as.numeric(veteran$trt) - 1

# Comprehensive survival analysis
result <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  arms = c("Treatment", "Control")
)

# View results
print(result$cox_results$cox_text)
print(result$zlogrank_text)
print(result$quantile_results)
```

### Weighted Log-Rank Tests

```r
# Fleming-Harrington (0,1) weights - emphasize late differences
result_fh <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 1)
)

# Magirr-Burman weights with 6-month cutoff
result_mb <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  scheme = "MB",
  scheme_params = list(mb_tstar = 6)
)
```

### Stratified Analysis

```r
# Stratified by cell type
result_strat <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  strata.name = "celltype"
)

print(result_strat$z.score_stratified)
```

### Kaplan-Meier Curves

```r
# Plot KM curves with risk table
plot_weighted_km(result, 
                 arms = c("Treatment", "Control"),
                 show.med = TRUE,
                 risk.set = TRUE)
```
  
### Survival Difference with Confidence Bands

```r
# Calculate KM difference with simultaneous confidence bands
km_diff <- KM_diff(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  draws.band = 1000
)

# Plot difference curve
plot(km_diff$at_points, km_diff$dhat, type = "s",
     xlab = "Time", ylab = "Survival Difference (Treatment - Control)")
lines(km_diff$at_points, km_diff$sb_lower, lty = 2)
lines(km_diff$at_points, km_diff$sb_upper, lty = 2)
abline(h = 0, col = "grey")
```

### Weighted Cox Model

```r
# Weighted Cox model with resampling-based inference
cox_fit <- cox_rhogamma(
  dfcount = result,
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 0.5),
  draws = 1000,
  verbose = TRUE
)

# Results with resampling-based CI
print(cox_fit$cox_text_star)
print(cox_fit$hr_ci_star)
```

### RMST Analysis

```r
# Cumulative RMST with confidence bands
rmst_result <- cumulative_rmst_bands(

  df = veteran,
  fit = km_diff,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  draws_sb = 1000
)

print(rmst_result$rmst_text)
```

### Baseline Characteristics Table

```r
# Create publication-ready baseline table
baseline_table <- create_baseline_table(
  data = veteran,
  treat_var = "treat",
  vars_continuous = c("age", "karno"),
  vars_categorical = "celltype",
  vars_binary = "prior",
  var_labels = c(
    age = "Age (years)",
    karno = "Karnofsky Score",
    celltype = "Cell Type",
    prior = "Prior Therapy"
  )
)

print(baseline_table)
```

## Typical Workflow

1. **Prepare data** with treatment indicator (0/1), time-to-event, and event indicator
2. **Run `df_counting()`** for comprehensive survival analysis
3. **Visualize** with `plot_weighted_km()` for Kaplan-Meier curves
4. **Calculate differences** with `KM_diff()` for treatment effect visualization
5. **Fit weighted Cox model** with `cox_rhogamma()` for hazard ratio estimation
6. **Compute RMST** with `cumulative_rmst_bands()` for absolute treatment effects

## Documentation

For detailed examples and methodological background, see the package vignette:

```r
# View the main vignette
vignette("weightedsurv_examples", package = "weightedsurv")
```

The vignette demonstrates:

- Analysis of the GBSG randomized trial data
- Propensity score weighted analysis of Rotterdam observational data
- Comparison of randomized and observational study results
- Simultaneous confidence bands for survival differences
- Subgroup analysis with overlaid ITT confidence bands
- Weighted Cox models with multiple weighting schemes
- Cross-trial comparisons with propensity score adjustment
- Validation against the `survival` and `adjustedCurves` packages

You can also browse vignettes online at the [CRAN package page](https://CRAN.R-project.org/package=weightedsurv).

## Dependencies

**Imports:**
- survival
- stats
- graphics

**Suggests:**
- gt (for formatted tables
- ggplot2 (for weight scheme visualization)

## References

Fleming, T. R. and Harrington, D. P. (1991). *Counting Processes and Survival Analysis*. Wiley.

Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests. *Statistics in Medicine*, 38(20), 3782-3790.

Schemper, M., Wakounig, S., and Heinze, G. (2009). The estimation of average hazard ratios by weighted Cox regression. *Statistics in Medicine*, 28(19), 2473-2489.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This package is licensed under the MIT License. See the LICENSE file for details.
