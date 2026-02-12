# Weighted log-rank cumulative statistics

Calculates cumulative weighted log-rank statistics for survival data.

## Usage

``` r
wlr_cumulative(
  df_weights,
  scheme,
  scheme_params = list(rho = 0, gamma = 0),
  return_cumulative = FALSE
)
```

## Arguments

- df_weights:

  Data frame with weights and risk/event counts.

- scheme:

  Weighting scheme.

- scheme_params:

  List of scheme parameters.

- return_cumulative:

  Logical; whether to return cumulative statistics.

## Value

List with score and variance.
