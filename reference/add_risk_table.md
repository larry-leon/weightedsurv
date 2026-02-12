# Add risk table annotation to KM plot

Adds risk set counts for two groups to a Kaplan-Meier plot.

## Usage

``` r
add_risk_table(
  risk.points,
  rpoints0,
  rpoints1,
  col.0,
  col.1,
  risk.cex,
  ymin,
  risk_offset,
  risk_delta,
  y.risk0,
  y.risk1
)
```

## Arguments

- risk.points:

  Numeric vector of risk time points.

- rpoints0:

  Numeric vector of risk set counts for group 0.

- rpoints1:

  Numeric vector of risk set counts for group 1.

- col.0:

  Color for group 0.

- col.1:

  Color for group 1.

- risk.cex:

  Numeric; text size for risk table.

- ymin:

  Numeric; minimum y value.

- risk_offset:

  Numeric; offset for risk table.

- risk_delta:

  Numeric; delta for risk table.

- y.risk0:

  Numeric; y position for group 0 risk table.

- y.risk1:

  Numeric; y position for group 1 risk table.

## Value

Invisibly returns NULL. Used for plotting side effects.
