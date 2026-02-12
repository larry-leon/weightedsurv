# Weighted Cox Model with Rho-Gamma Weights

Fits a weighted Cox proportional hazards model using flexible
time-dependent weights (e.g., Fleming-Harrington, Magirr-Burman).
Supports resampling-based inference for variance estimation and bias
correction.

## Usage

``` r
cox_rhogamma(
  dfcount,
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 0.5),
  draws = 0,
  alpha = 0.05,
  verbose = FALSE,
  lr.digits = 4
)
```

## Arguments

- dfcount:

  List; output from
  [`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  containing counting process data including time, delta, z, w_hat,
  survP_all, survG_all, etc.

- scheme:

  Character; weighting scheme. See
  [`wt.rg.S`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)
  for options. Default: "fh".

- scheme_params:

  List; parameters for the selected scheme. Default:
  `list(rho = 0, gamma = 0.5)`. Required parameters depend on scheme:

  - "fh": `rho`, `gamma`

  - "MB": `mb_tstar`

  - "custom_time": `t.tau`, `w0.tau`, `w1.tau`

- draws:

  Integer; number of resampling draws for variance estimation and bias
  correction. If 0, only asymptotic inference is performed. Default: 0.

- alpha:

  Numeric; significance level for confidence intervals. Default: 0.05.

- verbose:

  Logical; whether to print detailed output. Default: FALSE.

- lr.digits:

  Integer; number of decimal places for formatted output. Default: 4.

## Value

A list containing:

- fit:

  List with fitted model components:

  - `bhat`: Estimated log hazard ratio

  - `sig_bhat_asy`: Asymptotic standard error

  - `u.zero`: Score statistic at beta=0 (log-rank)

  - `z.score`: Standardized score statistic

  - `sig2_score`: Variance of score statistic

  - `wt_rg`: Vector of time-dependent weights

  - `bhat_debiased`: Bias-corrected estimate (if draws \> 0)

  - `sig_bhat_star`: Resampling-based standard error (if draws \> 0)

- hr_ci_asy:

  Data frame with asymptotic HR and CI

- hr_ci_star:

  Data frame with resampling-based HR and CI (if draws \> 0)

- cox_text_asy:

  Formatted string with HR and asymptotic CI

- cox_text_star:

  Formatted string with HR and resampling CI (if draws \> 0)

- z.score_debiased:

  Bias-corrected score statistic (if draws \> 0)

- zlogrank_text:

  Formatted log-rank test result

## Details

This function solves the weighted Cox partial likelihood score equation:
U(beta) = sum_i w_i K_i (dN_0/Y_0 - dN_1/Y_1) = 0

where K_i are time-dependent weights and Y_j, dN_j are risk sets and
event counts for group j.

When `draws > 0`, the function performs resampling to:

- Estimate finite-sample variance (more accurate than asymptotic)

- Compute bias correction for \\\hat{\beta}\\

- Provide improved confidence intervals for small samples

The score test at \\\beta=0\\ corresponds to the weighted log-rank test.

## Note

The treatment variable in `dfcount` must be coded as 0=control,
1=experimental.

## References

Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests.
Statistics in Medicine, 38(20), 3782-3790.

## See also

[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
for preprocessing
[`wt.rg.S`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)
for weighting schemes
[`cox_score_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_score_rhogamma.md)
for score function

Other survival_analysis:
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)

Other weighted_tests:
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)

## Examples

``` r
# First get counting process data
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

dfcount <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat"
)

# Fit weighted Cox model with FH(0,0.5) weights
fit <- cox_rhogamma(
  dfcount = dfcount,
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 0.5),
  draws = 1000,
  verbose = TRUE
)
#> z-statistic:  logrank (1-sided) p = 0.3167 
#> Hazard Ratio (HR): 0.912
#> 95% CI: [0.624, 1.332]
#> Hazard Ratio (HR): 0.914
#> 95% CI: [0.626, 1.335]

print(fit$cox_text_star)  # Resampling-based CI
#> [1] "HR = 0.912 (0.6261, 1.3351)"
print(fit$zlogrank_text)  # Weighted log-rank test
#> [1] "logrank (1-sided) p = 0.3167"

# Compare asymptotic and resampling CIs
print(fit$hr_ci_asy)
#>          beta  sig_bhat        hr     lower    upper
#> 1 -0.09211436 0.1932503 0.9120009 0.6244538 1.331957
print(fit$hr_ci_star)
#>          beta  sig_bhat        hr     lower    upper
#> 1 -0.08962629 0.1931806 0.9142728 0.6260949 1.335093
```
