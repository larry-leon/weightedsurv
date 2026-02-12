# Kaplan-Meier quantile calculation

Calculates the quantile time for a Kaplan-Meier curve.

## Usage

``` r
kmq_calculations(time_points, survival_probs, qprob = 0.5, type = "midpoint")
```

## Arguments

- time_points:

  Vector of time points.

- survival_probs:

  Vector of survival probabilities.

- qprob:

  Quantile probability (default 0.5).

- type:

  Calculation type (midpoint or min).

## Value

Estimated quantile time.
