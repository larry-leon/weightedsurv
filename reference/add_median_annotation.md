# Add median annotation to KM plot

Adds median survival annotation to a Kaplan-Meier plot.

## Usage

``` r
add_median_annotation(
  medians_df,
  med.digits,
  med.cex,
  med.font,
  xmed.fraction,
  ymed.offset
)
```

## Arguments

- medians_df:

  Data frame with quantile results. Should contain columns `quantile`,
  `lower`, `upper`, and `group`.

- med.digits:

  Integer; number of digits to display for median and confidence
  interval.

- med.cex:

  Numeric; text size for median annotation.

- med.font:

  Integer; font for median annotation.

- xmed.fraction:

  Numeric; fraction of the x-axis for annotation placement (e.g., 0.8
  for 80\\ to the right).

- ymed.offset:

  Numeric; offset from the top of the plot for annotation placement.

## Value

Invisibly returns NULL. Used for plotting side effects.
