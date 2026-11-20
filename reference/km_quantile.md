# Kaplan-Meier quantile and confidence interval

Calculates the quantile and confidence interval for a Kaplan-Meier
curve. The quantile is the time at which survival drops to `qprob` (see
[`kmq_calculations`](https://larry-leon.github.io/weightedsurv/reference/kmq_calculations.md));
the limits are the same quantity read from the log-transform pointwise
limits of the curve (upper limit capped at 1, NA where S = 0), as
[`survival::quantile.survfit`](https://rdrr.io/pkg/survival/man/quantile.survfit.html)
does with `conf.type = "log"`.

## Usage

``` r
km_quantile(
  time_points,
  survival_probs,
  se_probs = NULL,
  qprob = 0.5,
  type = c("midpoint", "min"),
  conf_level = 0.95
)
```

## Arguments

- time_points:

  Vector of time points.

- survival_probs:

  Vector of survival probabilities.

- se_probs:

  Standard errors of survival probabilities.

- qprob:

  Survival level (default 0.5): the quantile is the time at which S
  drops to `qprob`.

- type:

  Calculation type (midpoint or min).

- conf_level:

  Confidence level (default 0.95).

## Value

List with quantile and confidence interval.
