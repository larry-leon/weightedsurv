# Add legends to KM plot

Adds legends for Cox model, log-rank test, and arms to a Kaplan-Meier
plot.

## Usage

``` r
add_legends(
  dfcount,
  show.cox,
  cox.cex,
  put.legend.cox,
  show.logrank,
  logrank.cex,
  put.legend.lr,
  show_arm_legend,
  arms,
  col.1,
  col.0,
  ltys,
  lwds,
  arm.cex,
  put.legend.arms
)
```

## Arguments

- dfcount:

  List with results, typically output from a survival analysis function.
  Should contain elements such as `cox_results`, `z.score`, and
  `zlogrank_text`.

- show.cox:

  Logical; show Cox legend.

- cox.cex:

  Numeric; Cox legend size.

- put.legend.cox:

  Character; Cox legend position (e.g., "topright").

- show.logrank:

  Logical; show logrank legend.

- logrank.cex:

  Numeric; logrank legend size.

- put.legend.lr:

  Character; logrank legend position (e.g., "topleft").

- show_arm_legend:

  Logical; show arm legend.

- arms:

  Character vector of arm labels.

- col.1:

  Color for group 1.

- col.0:

  Color for group 0.

- ltys:

  Line types.

- lwds:

  Line widths.

- arm.cex:

  Numeric; arm legend size.

- put.legend.arms:

  Character; arm legend position (e.g., "left").

## Value

Invisibly returns NULL. Used for plotting side effects.
