# Plot KM curves for two groups with optional confidence intervals and censoring marks

Plots Kaplan-Meier survival curves for two groups, with options for
confidence intervals and censoring marks.

## Usage

``` r
plot_km_curves_counting(
  at_points,
  S0.KM,
  idx0,
  idv0,
  S1.KM,
  idx1,
  idv1,
  col.0,
  col.1,
  ltys,
  lwds,
  Xlab,
  Ylab,
  ylim,
  xlim,
  show.ticks = FALSE,
  cens0 = NULL,
  risk.points,
  risk.points.label,
  cens1 = NULL,
  se0.KM = NULL,
  se1.KM = NULL,
  conf.int = FALSE,
  conf_level = 0.95,
  censor.cex = 1,
  time.zero = 0,
  tpoints.add = c(0),
  ...
)
```

## Arguments

- at_points:

  Numeric vector of time points for plotting.

- S0.KM:

  Numeric vector of survival probabilities for group 0.

- idx0:

  Indices for censoring in group 0.

- idv0:

  Indices for events in group 0.

- S1.KM:

  Numeric vector of survival probabilities for group 1.

- idx1:

  Indices for censoring in group 1.

- idv1:

  Indices for events in group 1.

- col.0:

  Color for group 0.

- col.1:

  Color for group 1.

- ltys:

  Line types for groups.

- lwds:

  Line widths for groups.

- Xlab:

  X-axis label.

- Ylab:

  Y-axis label.

- ylim:

  Y-axis limits.

- xlim:

  X-axis limits.

- show.ticks:

  Logical; show censoring marks (default FALSE).

- cens0:

  Numeric vector of censoring times for group 0.

- risk.points:

  Numeric vector of risk time points.

- risk.points.label:

  Character vector of labels for risk time points.

- cens1:

  Numeric vector of censoring times for group 1.

- se0.KM:

  Numeric vector of standard errors for group 0.

- se1.KM:

  Numeric vector of standard errors for group 1.

- conf.int:

  Logical; show confidence intervals (default FALSE).

- conf_level:

  Confidence level (default 0.95).

- censor.cex:

  Numeric; censoring mark size (default 1.0).

- time.zero:

  Numeric; time zero value for risk table alignment (default 0).

- tpoints.add:

  Numeric vector; additional time points to include (default c(0)).

- ...:

  Additional arguments to plot.

## Value

Invisibly returns NULL. Used for plotting side effects.
