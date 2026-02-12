# Cumulative RMST bands for survival curves

Calculates cumulative Restricted Mean Survival Time (RMST) and
confidence bands for survival curves using resampling. Optionally plots
the cumulative RMST curve, pointwise confidence intervals, and
simultaneous confidence bands.

## Usage

``` r
cumulative_rmst_bands(
  df,
  fit,
  tte.name,
  event.name,
  treat.name,
  weight.name = NULL,
  draws_sb = 1000,
  xlab = "months",
  ylim_pad = 0.5,
  rmst_max_legend = "left",
  rmst_max_cex = 0.7,
  plot = TRUE
)
```

## Arguments

- df:

  Data frame containing survival data.

- fit:

  Survival fit object (output from KM_diff).

- tte.name:

  Name of time-to-event variable in `df`.

- event.name:

  Name of event indicator variable in `df`.

- treat.name:

  Name of treatment group variable in `df`.

- weight.name:

  Optional name of weights variable in `df`.

- draws_sb:

  Number of resampling draws for simultaneous bands (default: 1000).

- xlab:

  Label for x-axis (default: "months").

- ylim_pad:

  Padding for y-axis limits (default: 0.5).

- rmst_max_legend:

  Position for RMST legend (default: "left").

- rmst_max_cex:

  Text size for RMST legend (default: 0.7).

- plot:

  Logical; if TRUE, plot the results. Default is TRUE.

## Value

A list with elements:

- at_points:

  Time points used for RMST calculation

- rmst_time:

  Cumulative RMST estimates

- sig2_rmst_time:

  Variance of RMST estimates

- rmst_time_lower:

  Pointwise lower confidence interval

- rmst_time_upper:

  Pointwise upper confidence interval

- rmst_maxtau_ci:

  RMST and CI at maximum time

- rmst_text:

  Text summary for legend

- c_alpha_band:

  Critical value for simultaneous band

- rmst_time_sb_lower:

  Simultaneous band lower bound

- rmst_time_sb_upper:

  Simultaneous band upper bound
