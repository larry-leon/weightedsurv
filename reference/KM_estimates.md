# Kaplan-Meier Survival Estimates and Variance

Computes Kaplan-Meier survival estimates and their variances given risk
and event counts.

## Usage

``` r
KM_estimates(ybar, nbar, sig2w_multiplier = NULL)
```

## Arguments

- ybar:

  Vector of risk set sizes at each time point.

- nbar:

  Vector of event counts at each time point.

- sig2w_multiplier:

  Optional vector for variance calculation. If NULL, calculated
  internally.

## Value

List with survival estimates and variances.
