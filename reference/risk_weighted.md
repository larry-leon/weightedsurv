# Weighted risk set

Computes the weighted number at risk at a specified time.

## Usage

``` r
risk_weighted(x, y, w = rep(1, length(y)))
```

## Arguments

- x:

  Time point.

- y:

  Vector of event/censoring times.

- w:

  Weights (default 1).

## Value

Weighted number at risk at time x.
