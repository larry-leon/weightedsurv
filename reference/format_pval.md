# Format p-value for display

Formats a p-value for display, showing "\<eps" for small values.

## Usage

``` r
format_pval(pval, eps = 0.001, digits = 3)
```

## Arguments

- pval:

  Numeric p-value.

- eps:

  Threshold for small p-values.

- digits:

  Number of digits to display.

## Value

Formatted p-value as character.
