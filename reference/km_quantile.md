# Kaplan-Meier quantile and confidence interval

Calculates the quantile and confidence interval for a Kaplan-Meier
curve.

## Usage

``` r
km_quantile(
  time_points,
  survival_probs,
  se_probs = NULL,
  qprob = 0.5,
  type = c("midpoint", "min"),
  conf_level = 0.95
)
```

## Arguments

- time_points:

  Vector of time points.

- survival_probs:

  Vector of survival probabilities.

- se_probs:

  Standard errors of survival probabilities.

- qprob:

  Quantile probability (default 0.5).

- type:

  Calculation type (midpoint or min).

- conf_level:

  Confidence level (default 0.95).

## Value

List with quantile and confidence interval.
