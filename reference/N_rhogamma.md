# Weighted event count for Cox model

Calculates the weighted event count for the Cox model with rho-gamma
weights.

## Usage

``` r
N_rhogamma(x, error, delta, weight = 1)
```

## Arguments

- x:

  Numeric vector of time points.

- error:

  Numeric vector of error terms.

- delta:

  Numeric vector of event indicators.

- weight:

  Numeric vector of weights (default: 1).

## Value

Numeric value of weighted event count.
