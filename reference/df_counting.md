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
  censoring_allmarks = TRUE
)
```

## Arguments

- df:

  Data frame containing survival data.

- tte.name:

  Character; name of the time-to-event column in `df`.

- event.name:

  Character; name of the event indicator column in `df` (1=event,
  0=censored).

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

  Numeric in (0,1); quantile probability for KM quantile table. Default:
  0.5 (median).

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

## Value

A list with the following components:

- cox_results:

  List with Cox model results including hazard ratio, confidence
  interval, p-value, and formatted text

- logrank_results:

  List with log-rank test chi-square statistic, p-value, and formatted
  text

- z.score:

  Standardized weighted log-rank test statistic

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

  Weighted log-rank statistic and its variance

- riskpoints0, riskpoints1:

  Risk set counts at specified time points

- z.score_stratified:

  Stratified z-score (if stratified analysis)

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
handle tied event times. For stratified analyses, stratum-specific
estimates are computed and combined appropriately.

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
