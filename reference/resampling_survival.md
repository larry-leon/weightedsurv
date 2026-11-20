# Resampling Survival Curves for Confidence Bands

Performs resampling to generate survival curves for a group, used for
constructing confidence bands.

## Usage

``` r
resampling_survival(U, W, D, at.points, draws.band, surv, G_draws)
```

## Arguments

- U:

  Vector of observed times.

- W:

  Vector of weights.

- D:

  Vector of event indicators (0/1).

- at.points:

  Vector of time points for evaluation.

- draws.band:

  Number of resampling draws.

- surv:

  Vector of survival estimates.

- G_draws:

  Matrix of random draws for resampling.

## Value

Matrix of resampled survival curves.

## Details

Each column is one resampled deviation \\\hat S^\*(t) - \hat S(t)\\ at
`at.points`, for one column of multipliers `G_draws` (one row per
subject). With no weights or all weights equal, the multipliers perturb
the events only: \\-\hat S(t) \sum\_{s \le t} \sum_i G_i w_i dN_i(s) /
Y_w(s)\\; these results are unchanged. With weights that are not all
equal, the deviation is built from the per-subject influence
contributions \\\psi_i(t)\\ of the robust variance (see
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
section "Variance of the Kaplan-Meier estimate"): \\-\hat S(t) \sum_i
G_i w_i \psi_i(t)\\, so its variance over draws is the robust variance
\\\hat S(t)^2 \sum_i w_i^2 \psi_i(t)^2\\, which equals that of
`survival::survfit(weights = w, robust = TRUE)`. The weights are treated
as fixed and known. Both forms need `at.points` to hold every event time
of the group
([`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
arranges this).
