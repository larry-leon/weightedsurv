# Weighted counting process

Computes the weighted count of events up to a specified time.

## Usage

``` r
count_weighted(x, y, w = rep(1, length(y)))
```

## Arguments

- x:

  Time point.

- y:

  Vector of event/censoring times.

- w:

  Weights (default 1).

## Value

Weighted count of events up to time x.
