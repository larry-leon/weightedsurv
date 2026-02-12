# Root-finding for Cox score function

Finds the root of the Cox model score equation for weighted log-rank
statistics.

## Usage

``` r
find_cox_root(time, delta, z, w_hat, wt_rg)
```

## Arguments

- time:

  Numeric vector of event times.

- delta:

  Numeric vector of event indicators.

- z:

  Numeric vector of group indicators.

- w_hat:

  Numeric vector of weights.

- wt_rg:

  Numeric vector of rho-gamma weights.

## Value

List with root and additional information, or NA if not found.
