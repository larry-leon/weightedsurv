# Safe execution wrapper

Executes an R expression safely, returning NULL and printing an error
message if an error occurs.

## Usage

``` r
safe_run(expr)
```

## Arguments

- expr:

  An R expression to evaluate.

## Value

The result of expr, or NULL if an error occurs.
