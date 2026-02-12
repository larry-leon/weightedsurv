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
