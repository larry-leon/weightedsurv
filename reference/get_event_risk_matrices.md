# Event and Risk Matrices for Survival Analysis

Constructs matrices indicating event and risk status for each subject at
specified time points.

## Usage

``` r
get_event_risk_matrices(U, at.points)
```

## Arguments

- U:

  Vector of observed times (e.g., time-to-event).

- at.points:

  Vector of time points at which to evaluate events and risk.

## Value

A list with event and risk matrices.
