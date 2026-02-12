# Calculate risk set and event counts at time points

Calculates risk set and event counts for a group at specified time
points, with variance estimation.

## Usage

``` r
calculate_risk_event_counts(U, D, W, at_points, draws = 0, seedstart = 816951)
```

## Arguments

- U:

  Numeric vector of times for group.

- D:

  Numeric vector of event indicators for group.

- W:

  Numeric vector of weights for group.

- at_points:

  Numeric vector of time points.

- draws:

  Number of draws for variance estimation (default 0).

- seedstart:

  Random seed for draws (default 816951).

## Value

List with ybar (risk set counts), nbar (event counts), sig2w_multiplier
(variance term).
