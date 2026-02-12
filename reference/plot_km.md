# Plot Kaplan-Meier curves

Plots Kaplan-Meier survival curves for groups in the data.

## Usage

``` r
plot_km(df, tte.name, event.name, treat.name, weights = NULL, ...)
```

## Arguments

- df:

  Data frame containing survival data.

- tte.name:

  Name of time-to-event column.

- event.name:

  Name of event indicator column.

- treat.name:

  Name of treatment/group column.

- weights:

  Optional; name of weights column.

- ...:

  Additional arguments passed to plot().

## Value

Kaplan-Meier fit object (invisible).
