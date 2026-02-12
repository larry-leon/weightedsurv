# Check and Compare Statistical Test Results

Calculates and compares three equivalent test statistics from weighted
survival analysis: the squared standardized weighted log-rank statistic,
the log-rank chi-squared statistic, and the squared Cox model z-score.
These should be approximately equal under correct implementation.

## Usage

``` r
check_results(dfcount, verbose = TRUE)
```

## Arguments

- dfcount:

  A list or data frame from
  [`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  containing:

  lr

  :   Weighted log-rank statistic

  sig2_lr

  :   Variance of the weighted log-rank statistic

  z.score

  :   Standardized z-score from Cox model

  logrank_results

  :   List containing `chisq`, the log-rank chi-squared statistic

- verbose:

  Logical; if `TRUE`, prints the comparison table to the console.
  Default: `TRUE`.

## Value

A data frame with one row and three columns, returned invisibly:

- zlr_sq:

  Squared standardized weighted log-rank: \\(lr / \sqrt{sig2\\lr})^2\\

- logrank_chisq:

  Chi-squared statistic from log-rank test

- zCox_sq:

  Squared z-score from Cox model: \\z.score^2\\

## Details

This function serves as a diagnostic tool to verify computational
consistency. The three statistics should be numerically equivalent
(within rounding error): \$\$(lr / \sqrt{sig2\\lr})^2 \approx
logrank\\chisq \approx z.score^2\$\$

Discrepancies between these values may indicate:

- Numerical instability in variance estimation

- Incorrect weighting scheme application

- Data processing errors

## Note

This function is primarily used for package development and validation.
End users typically don't need to call it directly.

## See also

[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
for generating the input object

## Examples

``` r
# \donttest{
# After running df_counting
library(survival)
data(veteran)
veteran$treat <- as.numeric(veteran$trt) - 1

result <- df_counting(
  df = veteran,
  tte.name = "time",
  event.name = "status",
  treat.name = "treat"
)

# Check consistency of test statistics
check_results(result)
#> 
#> Test Statistic Comparison:
#> (These should be approximately equal)
#> 
#>       zlr_sq logrank_chisq     zCox_sq
#>  0.008227343   0.008227343 0.008227343
#> 
#> Relative difference: 1.423e-11%

# Store results without printing
stats_comparison <- check_results(result, verbose = FALSE)
print(stats_comparison)
#>        zlr_sq logrank_chisq     zCox_sq
#> 1 0.008227343   0.008227343 0.008227343
# }

# Simple example with constructed data
dfcount_example <- list(
  lr = 2.5,
  sig2_lr = 1.0,
  z.score = 2.5,
  logrank_results = list(chisq = 6.25)
)
check_results(dfcount_example)
#> 
#> Test Statistic Comparison:
#> (These should be approximately equal)
#> 
#>  zlr_sq logrank_chisq zCox_sq
#>    6.25          6.25    6.25
#> 
#> Relative difference: 0%
```
