# Resampling for Weighted Cox Model (rho, gamma)

Performs resampling to estimate uncertainty for the weighted Cox model
(rho, gamma).

## Usage

``` r
cox_rhogamma_resample(
  fit_rhogamma,
  i_bhat,
  K_wt_rg,
  i_zero,
  K_zero,
  G1.draws = NULL,
  G0.draws = NULL,
  draws = 100,
  seedstart = 8316951
)
```

## Arguments

- fit_rhogamma:

  List with fitted Cox model results.

- i_bhat:

  Information at estimated beta.

- K_wt_rg:

  Weights at estimated beta.

- i_zero:

  Information at beta=0.

- K_zero:

  Weights at beta=0.

- G1.draws:

  Optional: pre-generated random draws for groups.

- G0.draws:

  Optional: pre-generated random draws for groups.

- draws:

  Number of resampling iterations (default: 100).

- seedstart:

  Random seed for reproducibility (default: 8316951).

## Value

List with resampling results (score, beta, standard error, etc.).

## Details

The multipliers are standard normal, one row per subject (`G0.draws` for
controls, `G1.draws` for treated, in the order of `fit_rhogamma$time`),
drawn under `set.seed(seedstart)`; the caller's random-number stream is
restored afterwards. With `w_hat` absent or all equal, the multipliers
perturb the events only, as before, and these results are unchanged.
When `w_hat` is not all equal, each resampled score is built from the
per-subject score contributions described below: `score_star` is
\\-\sum_i G_i \psi_i(\hat\beta)\\ and `score_star_null` is \\-\sum_i G_i
\psi_i(0)\\, so their variance over draws is the variance of the score
used for `sig_bhat_asy` and the log-rank z-score (`K_wt_rg` and `K_zero`
are then not used). `bhat_center_star`, `Z_bstar` and `sig_bhat_star`
are formed from `score_star` as before. When `fit_rhogamma` comes from a
[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
result made with `weight_variance = "estimated"` (it then holds
`ew_parts`), the contributions carry the estimation of the weights,
\\\psi_i + IF_i' B\\ (section "Estimated weights" of
[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)).

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
