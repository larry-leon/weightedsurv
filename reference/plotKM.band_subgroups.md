# Plot Kaplan-Meier Survival Difference Curves with Subgroups and Confidence Bands

Plots the difference in Kaplan-Meier survival curves between two groups,
optionally including simultaneous confidence bands and subgroup curves.
Also displays risk tables for the overall population and specified
subgroups.

## Usage

``` r
plotKM.band_subgroups(
  df,
  tte.name,
  event.name,
  treat.name,
  weight.name = NULL,
  sg_labels = NULL,
  ltype = "s",
  lty = 1,
  draws = 20,
  lwd = 2,
  sg_colors = NULL,
  color = "lightgrey",
  ymax.pad = 0.01,
  ymin.pad = -0.01,
  taus = c(-Inf, Inf),
  yseq_length = 5,
  cex_Yaxis = 0.8,
  risk_cex = 0.8,
  by.risk = 6,
  risk.add = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  ymin.del = 0.035,
  y.risk1 = NULL,
  y.risk2 = NULL,
  ymin2 = NULL,
  risk_offset = NULL,
  risk.pad = 0.01,
  risk_delta = 0.0275,
  tau_add = NULL,
  time.zero.pad = 0,
  time.zero.label = 0,
  xlabel = NULL,
  ylabel = NULL,
  Maxtau = NULL,
  seedstart = 8316951,
  ylim = NULL,
  draws.band = 20,
  qtau = 0.025,
  show_resamples = FALSE,
  modify_tau = FALSE
)
```

## Arguments

- df:

  Data frame containing survival data.

- tte.name:

  Name of the time-to-event column.

- event.name:

  Name of the event indicator column (0/1).

- treat.name:

  Name of the treatment group column (0/1).

- weight.name:

  Optional name of the weights column.

- sg_labels:

  Character vector of subgroup definitions (as logical expressions).

- ltype:

  Line type for curves (default: "s").

- lty:

  Line style for curves (default: 1).

- draws:

  Number of draws for resampling (default: 20).

- lwd:

  Line width for curves (default: 2).

- sg_colors:

  Colors for subgroup curves.

- color:

  Color for confidence band polygon (default: "lightgrey").

- ymax.pad:

  Padding for y-axis limits.

- ymin.pad:

  Padding for y-axis limits.

- taus:

  Vector for time truncation (default: c(-Inf, Inf)).

- yseq_length:

  Number of y-axis ticks (default: 5).

- cex_Yaxis:

  Text size for axis.

- risk_cex:

  Text size for risk table.

- by.risk:

  Interval for risk table time points (default: 6).

- risk.add:

  Additional time points for risk table.

- xmax:

  Additional graphical and calculation parameters.

- ymin:

  Additional graphical and calculation parameters.

- ymax:

  Additional graphical and calculation parameters.

- ymin.del:

  Additional graphical and calculation parameters.

- y.risk1:

  Additional graphical and calculation parameters.

- y.risk2:

  Additional graphical and calculation parameters.

- ymin2:

  Additional graphical and calculation parameters.

- risk_offset:

  Additional graphical and calculation parameters.

- risk.pad:

  Additional graphical and calculation parameters.

- risk_delta:

  Additional graphical and calculation parameters.

- tau_add:

  Additional graphical and calculation parameters.

- time.zero.pad:

  Additional graphical and calculation parameters.

- time.zero.label:

  Additional graphical and calculation parameters.

- xlabel:

  Additional graphical and calculation parameters.

- ylabel:

  Additional graphical and calculation parameters.

- Maxtau:

  Additional graphical and calculation parameters.

- seedstart:

  Additional graphical and calculation parameters.

- ylim:

  Additional graphical and calculation parameters.

- draws.band:

  Number of draws for simultaneous confidence bands (default: 20).

- qtau:

  Quantile for time range in simultaneous bands (default: 0.025).

- show_resamples:

  Logical; whether to plot resampled curves (default: FALSE).

- modify_tau:

  Logical; restrict time range for bands.

## Value

(Invisible) list containing KM_diff results, time points, subgroup
curves, risk tables, and confidence intervals.
