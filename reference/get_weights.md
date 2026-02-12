# Get weights for a weighting scheme

Calculates weights for a specified scheme at given time points.

## Usage

``` r
get_weights(scheme, scheme_params, S.pool, tpoints)
```

## Arguments

- scheme:

  Character string specifying the weighting scheme.

- scheme_params:

  List of parameters for the scheme.

- S.pool:

  Numeric vector of pooled survival probabilities.

- tpoints:

  Numeric vector of time points.

## Value

Numeric vector of weights.
