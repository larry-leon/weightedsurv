# Extract and calculate weights for multiple schemes

Extracts and calculates weights for multiple schemes and returns a
combined data frame.

## Usage

``` r
extract_and_calc_weights(atpoints, S.pool, weights_spec_list)
```

## Arguments

- atpoints:

  Numeric vector of time points.

- S.pool:

  Numeric vector of pooled survival probabilities.

- weights_spec_list:

  List of weighting scheme specifications.

## Value

Data frame with weights for each scheme.
