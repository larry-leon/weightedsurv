# Plot weight schemes for survival analysis

Plots the weights for different schemes over time using ggplot2.

## Usage

``` r
plot_weight_schemes(
  dfcount,
  tte.name = "time_months",
  event.name = "status",
  treat.name = "treat",
  arms = c("treat", "control"),
  weights_spec_list = list(`MB(12)` = list(scheme = "MB", mb_tstar = 12), `MB(6)` =
    list(scheme = "MB", mb_tstar = 6), `FH(0,1)` = list(scheme = "fh", rho = 0, gamma =
    1), `FH(0.5,0.5)` = list(scheme = "fh", rho = 0.5, gamma = 0.5), custom_time =
    list(scheme = "custom_time", t.tau = 25, w0.tau = 1, w1.tau = 1.5), FHexp2 =
    list(scheme = "fh_exp2")),
  custom_colors = c(`FH(0,1)` = "grey", FHexp2 = "black", `MB(12)` = "#1b9e77", `MB(6)` =
    "#d95f02", `FH(0.5,0.5)` = "#7570b3", custom_time = "green"),
  custom_sizes = c(`FH(0,1)` = 2, FHexp2 = 1, `MB(12)` = 1, `MB(6)` = 1, `FH(0.5,0.5)` =
    1, custom_time = 1),
  transform_fh = FALSE,
  rescheme_fhexp2 = TRUE
)
```

## Arguments

- dfcount:

  Data frame containing counting process results, including time points
  and survival probabilities.

- tte.name:

  Name of the time-to-event variable (default: 'time_months').

- event.name:

  Name of the event indicator variable (default: 'status').

- treat.name:

  Name of the treatment group variable (default: 'treat').

- arms:

  Character vector of group names (default: c('treat', 'control')).

- weights_spec_list:

  List of weighting scheme specifications.

- custom_colors:

  Named character vector of colors for each scheme.

- custom_sizes:

  Named numeric vector of line sizes for each scheme.

- transform_fh:

  Logical; whether to transform FH weights (default: FALSE).

- rescheme_fhexp2:

  Logical; whether to rescheme FHexp2 and custom_time (default: TRUE).

## Value

A ggplot object showing the weight schemes over time.

## Details

This function visualizes the weights used in various survival analysis
schemes (e.g., FH, MB, custom) using ggplot2. Facets and colors are
customizable.
