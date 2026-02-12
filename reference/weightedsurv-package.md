# Weighted Survival Analysis Tools

Provides functions for weighted and stratified survival analysis,
including Cox proportional hazards models, weighted log-rank tests,
Kaplan-Meier curves, and RMST calculations with various weighting
schemes.

## Main Functions

- [`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md):
  Main function for weighted survival analysis

- [`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md):
  Kaplan-Meier difference calculations

- [`cox_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md):
  Weighted Cox model with flexible weights

- [`plot_weight_schemes`](https://larry-leon.github.io/weightedsurv/reference/plot_weight_schemes.md):
  Visualize weighting schemes

- [`cumulative_rmst_bands`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md):
  RMST analysis with confidence bands

## Weighting Schemes

The package supports multiple weighting schemes for log-rank tests:

- Fleming-Harrington (fh): Emphasizes early or late differences

- Magirr-Burman (MB): Modest downweighting after cutoff time

- Schemper: Adjusts for censoring distribution

- Xu-O'Quigley (XO): Adjusts for risk set size

- Custom time-based weights

## Typical Workflow

1.  Prepare data with treatment (0/1), time-to-event, and event
    indicator

2.  Run
    [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
    for comprehensive analysis

3.  Use
    [`plot_weighted_km()`](https://larry-leon.github.io/weightedsurv/reference/plot_weighted_km.md)
    to visualize survival curves

4.  Calculate survival differences with
    [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)

5.  Perform weighted Cox regression with
    [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)

6.  Compute RMST with
    [`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)

## Key Features

- Flexible weighting schemes for hypothesis testing

- Resampling-based inference for improved small-sample properties

- Simultaneous confidence bands for survival curves

- Stratified analysis support

- Weighted observations

- Comprehensive diagnostic checks

## References

Fleming, T. R. and Harrington, D. P. (1991). Counting Processes and
Survival Analysis. Wiley.

Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests.
Statistics in Medicine, 38(20), 3782-3790.

## See also

Useful links:

- <https://larry-leon.github.io/weightedsurv/>

- <https://github.com/larry-leon/weightedsurv>

- Report bugs at <https://github.com/larry-leon/weightedsurv/issues>

## Author

**Maintainer**: Larry Leon <larry.leon.05@post.harvard.edu>

## Examples

``` r
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

# Basic analysis
result <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat"
)

# Plot results
plot_weighted_km(result)


# Weighted log-rank emphasizing late differences
result_fh <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat",
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 1)
)
```
