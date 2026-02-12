# Get validated weights for a data frame

Validates and returns weights for a data frame according to specified
schemes.

## Usage

``` r
get_validated_weights(
  df_weights,
  scheme = "fh",
  scheme_params = list(rho = 0, gamma = 0),
  details = FALSE
)
```

## Arguments

- df_weights:

  Data frame containing weights and related data.

- scheme:

  Character string specifying the weighting scheme.

- scheme_params:

  List of parameters for the scheme.

- details:

  Logical; if TRUE, returns detailed output.

## Value

Numeric vector or list of validated weights.
