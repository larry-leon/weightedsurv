# Score calculation for weighted Cox model

Calculates the score and variance for the weighted Cox model.

## Usage

``` r
score_calculation(ybar1, ybar0, dN1, dN0, wt_rg)
```

## Arguments

- ybar1:

  Numeric vector of event counts for group 1.

- ybar0:

  Numeric vector of event counts for group 0.

- dN1:

  Numeric vector of event increments for group 1.

- dN0:

  Numeric vector of event increments for group 0.

- wt_rg:

  Numeric vector of rho-gamma weights.

## Value

List with score, variance, information, and weights.
