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

- weight_variance:

  `"known"`, or `"estimated"` when `dfcount` was made with
  `weight_variance = "estimated"`

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

With subject weights that are not all equal (`dfcount$w_hat`), the
variance of the score at 0 and at \\\hat\beta\\ is the sum of squared
per-subject score contributions (section below), so `fit$sig2_score`,
`fit$z.score`, `fit$sig_bhat_asy`, `hr_ci_asy`, `zlogrank_text` and
`cox_text_asy` differ from versions before 0.2.0; with `draws > 0` the
resampled scores are built from the same contributions
([`cox_rhogamma_resample`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma_resample.md)),
so the resampling standard error estimates the same quantity as
`sig_bhat_asy`. The estimate `bhat`, the scores `u.beta` and `u.zero`
and the time weights `wt_rg` do not change. `fit$z.score` equals
`df_counting()$z.score` for the same scheme. With no weights or all
weights equal, every result is as before.

When `dfcount` was made with `weight_variance = "estimated"`
([`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
with a
[`ps_weights`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
object), the estimation of the weights is carried into `fit$sig2_score`,
`fit$z.score`, `z.score`, `fit$sig_bhat_asy`, `hr_ci_asy`,
`zlogrank_text` and `cox_text_asy`, and with `draws > 0` into the
resampled scores and everything computed from them (section "Estimated
weights"). `bhat`, the scores, the information and the time weights are
unchanged. The result's element `weight_variance` says which option was
used. A direct call to
[`cox_score_rhogamma`](https://larry-leon.github.io/weightedsurv/reference/cox_score_rhogamma.md)
is not affected.

## Note

The treatment variable in `dfcount` must be coded as 0=control,
1=experimental.

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
or
[`plotKM.band_subgroups`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md).
`cox_rhogamma`,
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
  `cox_rhogamma()` are.

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
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md),
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
