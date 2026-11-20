# Plot confidence interval polygon for KM curve

Plots a shaded polygon representing the confidence interval for a
Kaplan-Meier survival curve. Limits use the log transform; the upper
limit is capped at 1 and, where survival is 0, the band closes onto the
curve.

## Usage

``` r
plot_km_confint_polygon(x, surv, se, conf_level, col)
```

## Arguments

- x:

  Numeric vector of time points.

- surv:

  Numeric vector of survival probabilities.

- se:

  Numeric vector of standard errors of survival probabilities.

- conf_level:

  Numeric; confidence level for interval (default 0.95).

- col:

  Color for the polygon.

## Value

Invisibly returns NULL. Used for plotting side effects.
