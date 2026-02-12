# Check KM curve for validity

Checks that a Kaplan-Meier curve is valid (values in \[0,1\],
non-increasing).

## Usage

``` r
check_km_curve(S.KM, group_name = "Group", stop_on_error = TRUE)
```

## Arguments

- S.KM:

  Numeric vector of survival probabilities.

- group_name:

  Character; name of the group.

- stop_on_error:

  Logical; whether to stop on error (default TRUE).

## Value

None. Stops or warns if invalid.
