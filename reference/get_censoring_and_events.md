# Get censoring and event times and their indices

Extracts censoring and event times and their indices for a group at
specified time points.

## Usage

``` r
get_censoring_and_events(time, delta, z, group, censoring_allmarks, at_points)
```

## Arguments

- time:

  Numeric vector of times.

- delta:

  Numeric vector of event indicators.

- z:

  Numeric vector of group indicators.

- group:

  Value of group to extract.

- censoring_allmarks:

  Logical; if FALSE, remove events from censored.

- at_points:

  Numeric vector of time points.

## Value

List with cens (censored times), ev (event times), idx_cens, idx_ev,
idx_ev_full.
