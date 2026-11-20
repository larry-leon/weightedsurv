# Kaplan-Meier Survival Estimates and Variance

Computes Kaplan-Meier survival estimates and their variances given risk
and event counts. Survival is returned in \[0, 1\]; when everything
still at risk fails at a time point, survival is exactly 0 from there on
and so is its variance.

## Usage

``` r
KM_estimates(ybar, nbar, sig2w_multiplier = NULL)
```

## Arguments

- ybar:

  Vector of risk set sizes at each time point.

- nbar:

  Vector of cumulative event counts at each time point.

- sig2w_multiplier:

  Optional vector for variance calculation. If NULL, calculated
  internally.

## Value

List with survival estimates and variances.

## Details

The result is the Kaplan-Meier estimate only when the time points behind
`ybar` and `nbar` include every event time of the group (see
[`KM_diff`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
which arranges this for any `at_points`).
