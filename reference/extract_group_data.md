# Extract time, event, and weight data for a group

Extracts time, event, and weight vectors for a specified group.

## Usage

``` r
extract_group_data(time, delta, wgt, z, group = 1)
```

## Arguments

- time:

  Numeric vector of times.

- delta:

  Numeric vector of event indicators (1=event, 0=censored).

- wgt:

  Numeric vector of weights.

- z:

  Numeric vector of group indicators.

- group:

  Value of group to extract (default 1).

## Value

List with U (times), D (events), W (weights).
