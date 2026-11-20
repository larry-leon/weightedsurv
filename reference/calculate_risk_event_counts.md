# Calculate risk set and event counts at time points

Calculates risk set and event counts for a group at specified time
points, with variance estimation. Passed to
[`KM_estimates`](https://larry-leon.github.io/weightedsurv/reference/KM_estimates.md),
these give the Kaplan-Meier estimate only when `at_points` holds every
event time of the group.

## Usage

``` r
calculate_risk_event_counts(U, D, W, at_points, draws = 0, seedstart = 816951)
```

## Arguments

- U:

  Numeric vector of times for group.

- D:

  Numeric vector of event indicators for group.

- W:

  Numeric vector of weights for group.

- at_points:

  Numeric vector of time points.

- draws:

  Number of draws for variance estimation (default 0).

- seedstart:

  Random seed for draws (default 816951).

## Value

List with ybar (risk set counts), nbar (event counts), sig2w_multiplier
(variance term: its cumulative sum times \\\hat S(t)^2\\ is the variance
of the Kaplan-Meier estimate, see
[`KM_estimates`](https://larry-leon.github.io/weightedsurv/reference/KM_estimates.md)).

## Details

`sig2w_multiplier` is, with no weights or all weights equal, the
per-time-point Greenwood term \\dN(s) / \[Y(s)\\Y(s) - dN(s)\\\]\\ (0
where \\Y(s) = 0\\ or \\Y(s) - dN(s) \le 10^{-10} Y(s)\\). With weights
that are not all equal it is the successive differences, along
`at_points`, of \\\sum_i w_i^2 \psi_i(t)^2\\ (the robust variance below
divided by \\\hat S(t)^2\\), or, when `draws > 0`, of the variance over
draws of \\\sum_i G_i w_i \psi_i(t)\\; individual differences may be
negative, only their cumulative sum is meaningful.

## Variance of the Kaplan-Meier estimate

With no weights, or with all weights equal, the variance of a
Kaplan-Meier estimate is Greenwood's formula, and these results are
unchanged. When the weights are not all equal it is the robust
(infinitesimal-jackknife) variance \$\$\widehat{Var}\[\hat S(t)\] = \hat
S(t)^2 \sum_i w_i^2 \psi_i(t)^2,\$\$ \$\$\psi_i(t) = \frac{D_i 1(U_i \le
t)}{Y_w(U_i) - dN_w(U_i)} - \sum\_{s \le \min(t, U_i)}
\frac{dN_w(s)}{Y_w(s)\\Y_w(s) - dN_w(s)\\},\$\$ where \\U_i\\, \\D_i\\
and \\w_i\\ are subject \\i\\'s time, event indicator and weight,
\\Y_w(s)\\ is the weighted number at risk, \\dN_w(s)\\ the weighted
number of events at \\s\\, the sum runs over event times, and a term is
0 where \\Y_w(s) - dN_w(s) \le 10^{-10} Y_w(s)\\. It equals the standard
error of `survival::survfit(weights = w, robust = TRUE)` (on the
time-fixed data), reduces to Greenwood's formula at equal weights, and
does not depend on the scale of the weights. With `draws > 0` and
unequal weights the variance is instead the variance over resampling
draws of \\\sum_i G_i w_i \psi_i(t)\\ (standard-normal multipliers
\\G_i\\), which estimates the same quantity; resampled curves for
weighted data are built from the same per-subject contributions (see
[`resampling_survival`](https://larry-leon.github.io/weightedsurv/reference/resampling_survival.md)).
By default the weights are treated as fixed and known. For weights
estimated by a logistic propensity model,
`weight_variance = "estimated"` carries their estimation into the
standard error (section "Estimated weights" of
[`df_counting`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)).
When the weights are extreme (low effective sample size) the standard
error can be too small under either option.
