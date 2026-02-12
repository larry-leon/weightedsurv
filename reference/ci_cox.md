# Confidence interval for Cox model estimate

Calculates the confidence interval for a Cox model log hazard ratio
estimate.

## Usage

``` r
ci_cox(bhat, sig_bhat, alpha = 0.05, verbose = FALSE)
```

## Arguments

- bhat:

  Estimated log hazard ratio.

- sig_bhat:

  Standard error of the estimate.

- alpha:

  Significance level (default: 0.05).

- verbose:

  Logical; if TRUE, prints interval.

## Value

Data frame with log hazard ratio, standard error, hazard ratio, lower
and upper confidence limits.
