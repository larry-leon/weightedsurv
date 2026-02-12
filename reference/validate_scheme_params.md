# Validate weighting scheme parameters

Checks and validates the parameters for a given weighting scheme.

## Usage

``` r
validate_scheme_params(scheme, scheme_params, S.pool)
```

## Arguments

- scheme:

  Character string specifying the weighting scheme.

- scheme_params:

  List of parameters for the scheme.

- S.pool:

  Numeric vector of pooled survival probabilities.

## Value

Logical indicating if parameters are valid, or stops with error.
