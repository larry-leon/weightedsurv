# Kaplan-Meier quantile calculation

Calculates the quantile time for a Kaplan-Meier curve: the time at which
survival drops to `qprob` (`qprob = 0.5` is the median; `qprob = 0.25`
is the time at which 25\\

## Usage

``` r
kmq_calculations(time_points, survival_probs, qprob = 0.5, type = "midpoint")
```

## Arguments

- time_points:

  Vector of time points (non-decreasing).

- survival_probs:

  Vector of survival probabilities at `time_points`.

- qprob:

  Survival level (default 0.5): the quantile is the time at which S
  drops to `qprob`.

- type:

  Calculation type (midpoint or min).

## Value

Estimated quantile time.

## Details

With `type = "midpoint"` the rule of
[`survival::quantile.survfit`](https://rdrr.io/pkg/survival/man/quantile.survfit.html)
is used: the first time at which S(t) \<= qprob (comparisons within
`sqrt(.Machine$double.eps)`); when S(t) equals `qprob` on an interval,
the midpoint of that interval; when that plateau runs to the last time
point, the midpoint between its start and the last time point. NA when
the curve never reaches `qprob`. With `type = "min"`, the first time at
which S(t) \<= qprob. The curve is taken to start at S = 1 at time 0 and
to end at the last element of `time_points`; values may be NA (they are
skipped).
