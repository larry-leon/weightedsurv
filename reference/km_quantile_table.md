# Table of KM quantiles for two groups

Returns a data frame of quantiles and confidence intervals for two
groups.

## Usage

``` r
km_quantile_table(
  time_points,
  surv0,
  se0,
  surv1,
  se1,
  arms = c("treat", "control"),
  qprob = 0.5,
  type = c("midpoint", "min"),
  conf_level = 0.95
)
```

## Arguments

- time_points:

  Vector of time points.

- surv0:

  Survival probabilities for group 0.

- se0:

  Standard errors for group 0.

- surv1:

  Survival probabilities for group 1.

- se1:

  Standard errors for group 1.

- arms:

  Group labels.

- qprob:

  Quantile probability.

- type:

  Calculation type.

- conf_level:

  Confidence level.

## Value

Data frame of quantiles and CIs for each group.
