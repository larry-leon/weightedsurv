# Estimated propensity-score weights: statistical details and properties

## The two options

Propensity-score weights are estimated from the data. The argument
`weight_variance` of
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)
and
[`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
says how standard errors, intervals, bands and tests treat them.
[`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
and
[`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
take the option from the
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
result they are given.

- `weight_variance = "known"` (the default). The weights are treated as
  fixed constants. The variance of each estimator is the sum of squared
  per-subject contributions computed with the weights held fixed.
- `weight_variance = "estimated"`. The weights come from a logistic
  regression for the treatment, supplied as `weight_model` (the object
  returned by
  [`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)).
  The sampling variability of the fitted coefficients is carried into
  the variance by a linearization: each subject’s contribution receives
  an added term, the subject’s influence on the fitted coefficients
  times the derivative of the estimator in those coefficients.

Three facts hold for every estimator in the package.

1.  No point estimate differs between the two options: survival curves,
    their difference, the RMST, the log-rank score, the hazard ratio,
    counts and risk tables are identical.
2.  Under `"estimated"` the correction is analytic. The propensity model
    is fitted once.
3.  Resampling (`draws.band`, `draws_sb`, and `draws` in
    [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md))
    applies standard normal multipliers to the corrected per-subject
    contributions. **The propensity model is not refitted within a
    draw**, and the weights, the influence values and the derivatives
    are the same in every draw.

This article gives the derivation and the properties. Each statement
about large-sample behaviour is labelled as *proved here*, *cited*, or a
*heuristic argument*. Formulas labelled *background* are not computed by
the package. The displayed steps of a derivation (the two lemmas, the
expansions written with $`\approx`$ or $`o_p`$, and the projection
argument of section 5) are statements about the computed quantities;
they are background in the same sense.

## Setup and notation

The notation is that of the `weightedkm_methodology` article.

**Data.** There are $`n`$ subjects, the rows analysed. Subject $`i`$ has
treatment indicator $`A_i \in \{0, 1\}`$, baseline covariates $`X_i`$,
observed time $`T_i = \min(\tilde{T}_i, C_i)`$ and event indicator
$`\Delta_i = I(\tilde{T}_i \leq C_i)`$. The at-risk and counting
processes are $`Y_i(t) = I(T_i \geq t)`$ and
$`N_i(t) = I(T_i \leq t, \Delta_i = 1)`$.

**Weights.** With $`\hat{e}_i`$ the fitted probability of treatment
(section 3),
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
returns one of two weight types:

``` math
\texttt{"ipw"}:\ \ w_i = \frac{A_i}{\hat{e}_i} + \frac{1 - A_i}{1 - \hat{e}_i}, \qquad
\texttt{"stabilized"}:\ \ w_i = \hat{\pi}_{A_i} \left\{\frac{A_i}{\hat{e}_i} + \frac{1 - A_i}{1 - \hat{e}_i}\right\}, \qquad \hat{\pi}_1 = \frac{1}{n} \sum_{i=1}^{n} A_i, \quad \hat{\pi}_0 = 1 - \hat{\pi}_1 .
```

The weights are never truncated or rescaled.

**Weighted counts and curves.** For arm $`a`$,

``` math
\bar{Y}_a^w(s) = \sum_{i: A_i = a} w_i\, Y_i(s), \qquad d\bar{N}_a^w(s) = \sum_{i: A_i = a} w_i\, \Delta_i\, I(T_i = s), \qquad
\hat{S}^a(t) = \prod_{s \leq t} \left\{1 - \frac{d\bar{N}_a^w(s)}{\bar{Y}_a^w(s)}\right\},
```

the product running over the event times of the arm. The survival
difference is $`\hat{\Delta}(t) = \hat{S}^1(t) - \hat{S}^0(t)`$.

**Estimand.** Both weight types weight each arm to the covariate
distribution of all $`n`$ subjects. Under the identifying assumptions of
the `weightedkm_methodology` article, $`\hat{S}^a(t)`$ estimates the
counterfactual survival function $`S^a(t) = P\{\tilde{T}(a) > t\}`$ in
the population the rows were sampled from, and $`\hat{\Delta}(t)`$
estimates $`\Delta(t) = S^1(t) - S^0(t)`$. The two weight types give the
same curves, because $`\hat{S}^a(t)`$ does not change when the weights
of an arm are multiplied by a constant. The weighted log-rank statistic
and the weighted Cox estimate (section 4.4) are computed in the weighted
sample of both arms together; the log hazard ratio estimates the root
$`\beta_w^*`$ of the limit of the weighted score. The two weight types
weight the arms differently in the pooled risk set, so they give
different log-rank and Cox results. No other estimand (treated
population, overlap population) is implemented by
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md).

## The propensity model

**Model.**
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
fits, or takes, a logistic regression for the treatment,
`stats::glm(formula, family = binomial(), data = data)` with the default
control. With $`x_i`$ the row of the model matrix of subject $`i`$
($`p`$ columns, aliased columns dropped) and coefficients $`\gamma`$,

``` math
e_i(\gamma) = \frac{\exp(x_i^\top \gamma)}{1 + \exp(x_i^\top \gamma)}, \qquad \hat{e}_i = e_i(\hat{\gamma}).
```

**Score and information.** The score of subject $`i`$ and the
information are

``` math
s_i(\gamma) = x_i\, \{A_i - e_i(\gamma)\}, \qquad \mathcal{I}_n(\gamma) = \sum_{j=1}^{n} x_j x_j^\top\, e_j(\gamma)\, \{1 - e_j(\gamma)\}.
```

For the logistic model the observed and the expected information
coincide, and $`\mathcal{I}_n`$ does not involve $`A`$. The formulas
$`s_i`$ and $`\mathcal{I}_n`$ are computed inside the expression for
$`IF_i`$ below.

**Influence of the coefficients.** The estimate $`\hat{\gamma}`$ solves
$`\sum_i s_i(\gamma) = 0`$. With $`\gamma^*`$ its limit,

``` math
\hat{\gamma} - \gamma^* = \sum_{i=1}^{n} IF_i + o_p(n^{-1/2}), \qquad IF_i = \mathcal{I}_n(\hat{\gamma})^{-1}\, x_i\, (A_i - \hat{e}_i).
```

This is the linear representation of an M-estimator (*cited*: Lunceford
and Davidian (2004), section 3.1, apply it to the logistic propensity
model through the theory of M-estimation). It does not require the model
to be correct: when the model is wrong, $`\gamma^*`$ is the limit of the
fit and the representation still holds, because $`\mathcal{I}_n`$ does
not involve $`A`$ (*heuristic argument*; not proved here).
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
returns the $`n \times p`$ matrix `IF` with rows $`IF_i^\top`$.

**Derivative of the weights.** For the `"ipw"` weights, differentiating
$`1 / e_i(\gamma)`$ and $`1 / \{1 - e_i(\gamma)\}`$ with
$`\partial e_i / \partial \gamma = e_i (1 - e_i)\, x_i`$ gives (*proved
here*, by differentiation)

``` math
dw_i = \frac{\partial w_i}{\partial \gamma}\bigg|_{\hat{\gamma}} = -\frac{1 - \hat{e}_i}{\hat{e}_i}\, x_i \ \ (A_i = 1), \qquad dw_i = \frac{\hat{e}_i}{1 - \hat{e}_i}\, x_i \ \ (A_i = 0).
```

[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
returns the $`n \times p`$ matrix `dw` with rows $`dw_i^\top`$.

**Stabilized weights.** With `type = "stabilized"`, $`w_i`$ and $`dw_i`$
are those of `"ipw"` multiplied by $`\hat{\pi}_{A_i}`$, and `IF` is the
same. The proportion $`\hat{\pi}_a`$ is treated as a fixed number: **its
estimation is not included in the correction.** For the Kaplan–Meier
curves, their difference and the RMST this loses nothing: those
estimators do not depend on $`\hat{\pi}_a`$ (section 4.1, *proved
here*). For the log-rank score and the Cox estimate the variability of
$`\hat{\pi}_a`$ is not carried.

## Linearization, estimator by estimator

### Weighted Kaplan–Meier curve of one arm

**Contribution under known weights.** For subject $`i`$ of arm $`a`$,

``` math
\psi_i(t) = \frac{\Delta_i\, I(T_i \leq t)}{\bar{Y}_a^w(T_i) - d\bar{N}_a^w(T_i)} - \sum_{s \leq \min(t,\, T_i)} \frac{d\bar{N}_a^w(s)}{\bar{Y}_a^w(s)\bigl\{\bar{Y}_a^w(s) - d\bar{N}_a^w(s)\bigr\}},
```

the sum running over the event times of the arm; a term is 0 when
$`\bar{Y}_a^w(s) = 0`$ or
$`\bar{Y}_a^w(s) - d\bar{N}_a^w(s) \leq 10^{-10}\, \bar{Y}_a^w(s)`$.
Equivalently,

``` math
\psi_i(t) = \sum_{s \leq t} \frac{dN_i(s) - Y_i(s)\, d\bar{N}_a^w(s) / \bar{Y}_a^w(s)}{\bar{Y}_a^w(s) - d\bar{N}_a^w(s)},
```

a sum of the subject’s residual increments (*proved here*: the two
expressions are equal term by term).

**Lemma 1** (*proved here*). For a subject $`i`$ of arm $`a`$, and a
time $`t`$ at which no excluded term occurs,

``` math
\frac{\partial \hat{S}^a(t)}{\partial w_i} = -\hat{S}^a(t)\, \psi_i(t).
```

*Proof.*
$`\log \hat{S}^a(t) = \sum_{s \leq t} \log\{1 - d\bar{N}_a^w(s) / \bar{Y}_a^w(s)\}`$,
and
$`\partial \{d\bar{N}_a^w(s) / \bar{Y}_a^w(s)\} / \partial w_i = \Delta_i I(T_i = s) / \bar{Y}_a^w(s) - Y_i(s)\, d\bar{N}_a^w(s) / \bar{Y}_a^w(s)^2`$.
Multiplying by
$`-\bar{Y}_a^w(s) / \{\bar{Y}_a^w(s) - d\bar{N}_a^w(s)\}`$ and summing
over $`s \leq t`$ gives $`-\psi_i(t)`$. $`\square`$

**Two consequences** (*proved here*). $`\hat{S}^a(t)`$ is unchanged when
all weights of the arm are multiplied by a constant, so by Euler’s
relation for homogeneous functions

``` math
\sum_{i: A_i = a} w_i\, \psi_i(t) = 0 .
```

The contributions are centred, and the curve does not depend on
$`\hat{\pi}_a`$: the derivative of $`\hat{S}^a(t)`$ in a common factor
of the arm’s weights is
$`-\hat{S}^a(t) \sum_{i: A_i = a} w_i \psi_i(t) / \hat{\pi}_a = 0`$.

**Representation under known weights** (*cited*; *background*). With the
weights evaluated at the true propensity score, the weighted estimator
has the martingale representation of Deng and Wang (2025) (see also
Fleming and Harrington (1991) for the unweighted case),

``` math
\hat{S}^a(t) - S^a(t) = -S^a(t) \sum_{i: A_i = a} w_i\, \psi_i^{0}(t) + o_p(n^{-1/2}), \qquad \psi_i^{0}(t) = \int_0^t \frac{dN_i(s) - Y_i(s)\, d\Lambda^a(s)}{n\, y_a(s)},
```

with $`\Lambda^a`$ the cumulative hazard of $`S^a`$ and $`y_a(s)`$ the
limit of $`\bar{Y}_a^w(s) / n`$. $`\psi_i(t)`$ is $`\psi_i^{0}(t)`$ with
$`\Lambda^a`$ and $`n\, y_a`$ replaced by their estimates; by Lemma 1,
$`-\hat{S}^a(t)\, w_i\, \psi_i(t)`$ is also the infinitesimal-jackknife
value of subject $`i`$, the derivative of $`\hat{S}^a(t)`$ in a case
weight multiplying $`w_i`$. The variance
$`\hat{S}^a(t)^2 \sum_{i: A_i = a} w_i^2 \psi_i(t)^2`$ is the robust
variance of `survival::survfit(weights = w, robust = TRUE)`.

**Added term under estimated weights.** Write $`\hat{S}^a(t; \gamma)`$
for the curve computed with the weights $`w_i(\gamma)`$. By the chain
rule and Lemma 1 (*proved here*),

``` math
\frac{\partial \hat{S}^a(t; \gamma)}{\partial \gamma}\bigg|_{\hat{\gamma}} = \sum_{k: A_k = a} \frac{\partial \hat{S}^a(t)}{\partial w_k}\, dw_k = -\hat{S}^a(t)\, D_a(t), \qquad D_a(t) = \sum_{k: A_k = a} \psi_k(t)\, dw_k .
```

A first-order expansion of $`\hat{S}^a(t; \hat{\gamma})`$ around
$`\gamma^*`$, with the representation of $`\hat{\gamma} - \gamma^*`$ of
section 3, gives (*heuristic argument*: the remainder is not bounded
here)

``` math
\hat{S}^a(t; \hat{\gamma}) - S^a(t) \approx -S^a(t) \sum_{i=1}^{n} \Phi_i^a(t), \qquad \Phi_i^a(t) = I(A_i = a)\, w_i\, \psi_i(t) + IF_i^\top D_a(t) \quad (i = 1, \ldots, n).
```

$`\Phi_i^a(t)`$ is the corrected contribution. Every subject
contributes, also those in the other arm, through the model. The same
construction applies to any group $`g`$ of subjects (an arm within a
stratum; all subjects, for the pooled curve), with $`g`$ in place of the
arm.

### Survival difference

``` math
\hat{\Delta}(t) - \Delta(t) \approx -\sum_{i=1}^{n} \xi_i(t), \qquad \xi_i(t) = \hat{S}^1(t)\, \Phi_i^1(t) - \hat{S}^0(t)\, \Phi_i^0(t).
```

Under `"known"`, $`\xi_i(t)`$ involves subject $`i`$’s own arm only; the
two curves are sums over disjoint sets of independent subjects, and the
variance of the difference is the sum of the two arm variances. Under
`"estimated"` each subject enters both $`\Phi_i^1`$ and $`\Phi_i^0`$
through $`IF_i`$, because one fitted model supplies the weights of both
arms. The two curves are then correlated, and the variance of the
difference is **not** the sum of the arm variances. The cross term is
(*proved here*, by expanding the product and using
$`I(A_i = 1)\, I(A_i = 0) = 0`$)

``` math
\sum_{i=1}^{n} \Phi_i^1(t)\, \Phi_i^0(t) = E_1(t)^\top D_0(t) + E_0(t)^\top D_1(t) + D_1(t)^\top M\, D_0(t), \qquad E_a(t) = \sum_{i: A_i = a} w_i\, \psi_i(t)\, IF_i, \quad M = \sum_{i=1}^{n} IF_i\, IF_i^\top .
```

### RMST and cumulative RMST difference

The cumulative RMST difference is the integral of the survival
difference over the time grid $`t_0 < t_1 < \cdots < t_m`$ of the
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
fit,

``` math
\hat{R}(t_j) = \sum_{k \leq j} h_k\, (t_k - t_{k-1}), \qquad h_k = \tfrac{1}{2}\{\hat{\Delta}(t_{k-1}) + \hat{\Delta}(t_k)\} \ \text{(trapezoid, the default)} \quad \text{or} \quad h_k = \hat{\Delta}(t_{k-1}) \ \text{(rectangle)},
```

with $`\hat{R}(t_0) = 0`$. The map from the difference curve to
$`\hat{R}`$ is linear, so the linearization of section 4.2 carries over:
$`\hat{R}(t_j) - R(t_j) \approx -\sum_i \mathcal{L}_j(\xi_i)`$, with
$`\mathcal{L}_j`$ the same rule applied to $`\xi_i`$ (*heuristic
argument*, as in section 4.1). The package does not form
$`\mathcal{L}_j(\xi_i)`$ for each subject (*background*); it applies the
rule to the resampled difference curves (section 6).

### Weighted log-rank and Cox score

**Notation.** $`\kappa(s)`$ is the time weight of the scheme (written
$`w(t)`$ in the `weightedcox_methodology` article, where the subject
weight is $`a_i`$ and the arm $`Z_i`$). For a log hazard ratio
$`\beta`$,

``` math
S_0(s, \beta) = \sum_j w_j\, Y_j(s)\, e^{\beta A_j}, \qquad \bar{A}(s, \beta) = \frac{\sum_j w_j\, Y_j(s)\, A_j\, e^{\beta A_j}}{S_0(s, \beta)}, \qquad d\hat{\Lambda}(s, \beta) = \frac{\sum_j w_j\, \Delta_j\, I(T_j = s)}{S_0(s, \beta)} .
```

**Contribution under known weights.**

``` math
\psi_i(\beta) = w_i \left[ \Delta_i\, \kappa(T_i) \bigl\{A_i - \bar{A}(T_i, \beta)\bigr\} - e^{\beta A_i} \sum_{s \le T_i} \kappa(s) \bigl\{A_i - \bar{A}(s, \beta)\bigr\}\, d\hat{\Lambda}(s, \beta) \right],
```

the sum running over the distinct event times. The score, in the sign
convention of the package, and the information are

``` math
U_w(\beta) = -\sum_{i=1}^{n} \psi_i(\beta) = -\sum_{i=1}^{n} w_i\, \Delta_i\, \kappa(T_i) \bigl\{A_i - \bar{A}(T_i, \beta)\bigr\}, \qquad
I_w(\beta) = \frac{\partial U_w(\beta)}{\partial \beta} = \sum_{s} \kappa(s)\, \bar{A}(s, \beta) \bigl\{1 - \bar{A}(s, \beta)\bigr\}\, S_0(s, \beta)\, d\hat{\Lambda}(s, \beta) .
```

The second expression for $`U_w`$ holds because the compensator terms
sum to zero:
$`\sum_i w_i Y_i(s) e^{\beta A_i} \{A_i - \bar{A}(s, \beta)\} = 0`$ for
every $`s`$ (*proved here*). $`U_w(0)`$ is the weighted log-rank
statistic and $`\hat{\beta}_w`$ solves $`U_w(\beta) = 0`$.

**Lemma 2** (*proved here*). With the time weights $`\kappa(s)`$ held
fixed,

``` math
\frac{\partial U_w(\beta)}{\partial w_j} = -\frac{\psi_j(\beta)}{w_j} .
```

*Proof.* Differentiate the second expression for $`U_w`$. The term $`j`$
gives $`-\Delta_j \kappa(T_j) \{A_j - \bar{A}(T_j, \beta)\}`$. Through
$`\bar{A}`$, with
$`\partial \bar{A}(s, \beta) / \partial w_j = Y_j(s)\, e^{\beta A_j} \{A_j - \bar{A}(s, \beta)\} / S_0(s, \beta)`$,
the other terms give
$`\sum_s \kappa(s)\, d\hat{\Lambda}(s, \beta)\, Y_j(s)\, e^{\beta A_j} \{A_j - \bar{A}(s, \beta)\}`$.
The sum of the two is $`-\psi_j(\beta) / w_j`$. $`\square`$

**Corrected contribution.** By the chain rule and Lemma 2,
$`\partial U_w(\beta) / \partial \gamma = -B(\beta)`$ with

``` math
B(\beta) = \sum_{j=1}^{n} \frac{\psi_j(\beta)}{w_j}\, dw_j, \qquad \psi_i^{ew}(\beta) = \psi_i(\beta) + IF_i^\top B(\beta),
```

and the first-order expansion in $`\gamma`$ gives
$`U_w(\beta; \hat{\gamma}) \approx -\sum_i \psi_i^{ew, 0}(\beta)`$, a
sum of independent mean-zero terms that $`\psi_i^{ew}(\beta)`$ estimates
(*heuristic argument*).

**Test, standard error, interval.**

``` math
z = \frac{U_w(0)}{\bigl\{\sum_{i} \psi_i^{ew}(0)^2\bigr\}^{1/2}}, \qquad
\widehat{\text{se}}(\hat{\beta}_w) = \frac{\bigl\{\sum_{i} \psi_i^{ew}(\hat{\beta}_w)^2\bigr\}^{1/2}}{I_w(\hat{\beta}_w)}, \qquad
\exp\bigl\{\hat{\beta}_w \pm z_{1 - \alpha/2}\, \widehat{\text{se}}(\hat{\beta}_w)\bigr\} .
```

The standard error follows from expanding
$`0 = U_w(\hat{\beta}_w; \hat{\gamma})`$ in $`\beta`$:
$`\hat{\beta}_w - \beta_w^* \approx -U_w(\beta_w^*; \hat{\gamma}) / I_w(\beta_w^*)`$
(*heuristic argument*). Under `"known"` the same formulas hold with
$`\psi_i`$ in place of $`\psi_i^{ew}`$; for FH(0, 0) they are then the
robust score test and robust standard error of
`survival::coxph(weights = , robust = TRUE, ties = "breslow")` (Lin and
Wei 1989).

**Time-weight schemes.** The scheme enters only through $`\kappa(s)`$,
evaluated at the observed times. As coded in
[`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md),
with $`\hat{S}(s-)`$ the pooled weighted Kaplan–Meier curve just before
$`s`$, $`\hat{G}`$ the same for the censoring distribution and
$`\bar{Y}^w`$ the pooled weighted risk set:

| `scheme` | $`\kappa(s)`$ | depends on the subject weights |
|----|----|----|
| `"fh"` | $`\hat{S}(s-)^{\rho}\, \{1 - \hat{S}(s-)\}^{\gamma}`$ | yes, unless $`\rho = \gamma = 0`$ |
| `"MB"` | $`1 / \max\{\hat{S}(s-), \hat{S}(t^*)\}`$ | yes |
| `"schemper"` | $`\hat{S}(s-) / \hat{G}(s-)`$ | yes |
| `"XO"` | $`\hat{S}(s-) / \bar{Y}^w(s)`$ | yes |
| `"custom_time"` | $`w_0`$ for $`s \leq t_\tau`$, $`w_1`$ after | no |
| `"fh_exp1"`, `"fh_exp2"` | exponential variants of the FH weight | yes |

`"estimated"` runs for every scheme, with the formulas above. **The
dependence of $`\kappa(s)`$ on the estimated subject weights is
ignored**: Lemma 2 holds $`\kappa`$ fixed. For FH(0, 0) and
`"custom_time"` nothing is ignored. For the other schemes the omitted
term is not derived here and its size is not bounded here.

**Strata.** In a stratified
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
call, $`\psi_i(\beta)`$ is computed within the subject’s stratum,
$`B(\beta)`$ is the sum over all subjects of all strata, and the
variance is the single sum $`\sum_{i=1}^n \psi_i^{ew}(0)^2`$ over all
subjects. It is not a sum of per-stratum variances: the strata share the
coefficients of the model.

**Outputs of
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
from the `survival` package.** `logrank_results`, `cox_results$score`,
`$se`, `$hr_ci` and `$pval` carry the same correction, applied to a
vector $`r`$ of weighted score or dfbeta residuals of the
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) fit:

``` math
r_i^{ew} = r_i + IF_i^\top \sum_{j=1}^{n} \frac{r_j}{w_j}\, dw_j .
```

**Cumulative statistic and covariance with the survival difference.**
With $`c_i(t)`$ the contribution of subject $`i`$ to the score up to
$`t`$ ($`\psi_i(0)`$ when $`T_i \leq t`$, the compensator part
accumulated up to $`t`$ otherwise),
`wlr_cumulative(return_cumulative = TRUE)` returns

``` math
V^{ew}(t) = \sum_{i=1}^{n} \bigl\{c_i(t) + IF_i^\top B(t)\bigr\}^2, \qquad B(t) = \sum_{j=1}^{n} \frac{c_j(t)}{w_j}\, dw_j ,
```

and
[`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
returns, at the grid point read for `tzero`, the variance
$`\sum_i \xi_i(t)^2`$ of the survival difference and the covariance
$`\sum_i \psi_i^{ew}(0)\, \xi_i(t)`$ of the log-rank score with it.

## Variance

**Sum of squared corrected contributions.**

``` math
\widehat{\text{Var}}\bigl[\hat{S}^a(t)\bigr] = \hat{S}^a(t)^2 \sum_{i=1}^{n} \Phi_i^a(t)^2, \qquad
\widehat{\text{Var}}\bigl[\hat{\Delta}(t)\bigr] = \sum_{i=1}^{n} \xi_i(t)^2, \qquad
\widehat{\text{Var}}\bigl[U_w(\beta)\bigr] = \sum_{i=1}^{n} \psi_i^{ew}(\beta)^2 .
```

The variance of the RMST difference is taken over resampling draws
(section 6). That each sum estimates the variance of its estimator
follows from the linearizations of section 4, the contributions being
independent across subjects (*heuristic argument*).

**Relation to the `"known"` variance** (*proved here*, by expanding the
square). For contributions $`c_i`$ under known weights and a derivative
vector $`b`$ (for the score, $`c_i = \psi_i(\beta)`$ and
$`b = B(\beta)`$; for the survival difference,
$`c_i = \hat{S}^1 I(A_i = 1) w_i \psi_i - \hat{S}^0 I(A_i = 0) w_i \psi_i`$
and $`b = \hat{S}^1 D_1 - \hat{S}^0 D_0`$),

``` math
\sum_{i=1}^{n} \bigl(c_i + IF_i^\top b\bigr)^2 = \underbrace{\sum_{i=1}^{n} c_i^2}_{\text{"known"}} + 2\, b^\top C + b^\top M\, b, \qquad C = \sum_{i=1}^{n} c_i\, IF_i, \quad M = \sum_{i=1}^{n} IF_i\, IF_i^\top .
```

The last term is non-negative. The middle term can have either sign, so
in a given sample the `"estimated"` variance can be smaller or larger
than the `"known"` one.

**Why `"known"` tends to be conservative when the propensity model is
correct** (*heuristic argument*; the result for
inverse-probability-weighted means is *cited*). Let $`m_i(\gamma)`$ be
the population version of the known-weight contribution, computed with
the weights $`w_i(\gamma)`$, and let $`P_\gamma`$ be the distribution of
the data when the treatment follows the model with coefficients
$`\gamma`$. When the model is correct, weighting by the true propensity
score removes the dependence of the weighted sample on the treatment
mechanism, so $`E_\gamma\{m_i(\gamma)\} = 0`$ for every $`\gamma`$.
Differentiating this identity in $`\gamma`$ under the integral sign
gives

``` math
E\Bigl\{\frac{\partial m_i}{\partial \gamma}\Bigr\} = -E\bigl\{m_i\, s_i\bigr\},
```

with $`s_i`$ the score of the propensity model. The derivative $`b`$
estimates $`n\, E\{\partial m_i / \partial \gamma\}`$ and $`IF_i`$
estimates $`\{n\, E(s_i s_i^\top)\}^{-1} s_i`$, so

``` math
c_i + IF_i^\top b \ \approx\ m_i - E\bigl\{m_i\, s_i^\top\bigr\}\, \bigl\{E(s_i s_i^\top)\bigr\}^{-1} s_i ,
```

the residual of the projection of the known-weight contribution on the
space spanned by the score of the propensity model. A residual has
variance no larger than the variable projected:

``` math
\text{Var}\bigl(m_i - \Pi[m_i \mid s_i]\bigr) = \text{Var}(m_i) - E\bigl\{m_i\, s_i^\top\bigr\}\, \bigl\{E(s_i s_i^\top)\bigr\}^{-1} E\bigl\{s_i\, m_i\bigr\} \ \leq\ \text{Var}(m_i) .
```

The estimator computed with estimated weights therefore has a
large-sample variance no larger than the `"known"` formula estimates,
and the `"known"` intervals are too wide by the projected part. The two
agree when the known-weight contribution is uncorrelated with the score
of the propensity model. For inverse-probability-weighted means with a
logistic propensity model this is the result of Lunceford and Davidian
(2004) (section 3.1): the variance with the estimated propensity score
is the variance with the known one minus a quadratic form of this kind.
Hirano et al. (2003) show that weighting by an estimated propensity
score, in place of the true one, leads to an efficient estimate of the
average treatment effect. The steps not proved here are the interchange
of differentiation and expectation for the survival contributions, and
the convergence of $`b / n`$ and of the estimated contributions.

The argument uses the correctness of the model at one point only, the
identity $`E_\gamma\{m_i(\gamma)\} = 0`$. When the model is wrong, the
linearization of section 4 still describes the variability of the
estimator around its limit, but the ordering of the two variances is not
implied, and the limit is not the causal quantity.

The sample counterpart of the projected variance,
$`\sum_i c_i^2 - C^\top M^{-1} C`$, is *background*: the package does
not compute it. Section 8 shows it beside the package value.

## Multiplier resampling with corrected contributions

**Resampled survival difference.** Let $`G_1, \ldots, G_n`$ be
independent standard normal multipliers, independent of the data
(generated by
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
as one matrix for the controls and one for the treated, one row per
subject). One resampled deviation of the curve of arm $`a`$, and of the
difference, is

``` math
\hat{S}^{a*}(t) - \hat{S}^a(t) = -\hat{S}^a(t) \sum_{i=1}^{n} G_i\, \Phi_i^a(t) = -\hat{S}^a(t) \Bigl\{\sum_{i: A_i = a} G_i\, w_i\, \psi_i(t) + D_a(t)^\top \sum_{i=1}^{n} G_i\, IF_i\Bigr\}, \qquad
\hat{\Delta}^{*}(t) = -\sum_{i=1}^{n} G_i\, \xi_i(t).
```

The same multiplier $`G_i`$ multiplies the subject’s own-arm term and
its term through the model, in both arms. Conditionally on the data,
$`\hat{\Delta}^{*}`$ is a mean-zero Gaussian process with covariance
$`\sum_i \xi_i(s)\, \xi_i(t)`$; its variance at $`t`$ is exactly
$`\widehat{\text{Var}}[\hat{\Delta}(t)]`$ of section 5 (*proved here*:
the $`G_i`$ are independent with variance 1).

**Band for the survival difference.** With $`B`$ draws (`draws.band`)
and
$`\widehat{\text{se}}(t) = \widehat{\text{Var}}[\hat{\Delta}(t)]^{1/2}`$,

``` math
c_{1 - \alpha} = \text{the } (1 - \alpha) \text{ quantile over } b = 1, \ldots, B \text{ of } \max_{t} \frac{|\hat{\Delta}^{*(b)}(t)|}{\widehat{\text{se}}(t)}, \qquad \hat{\Delta}(t) \pm c_{1 - \alpha}\, \widehat{\text{se}}(t),
```

the maximum over the time grid of the fit (points where the ratio is
undefined are skipped).

**RMST.** Each resampled difference curve is integrated with the rule of
section 4.3,
$`\hat{R}^{*(b)}(t_j) = \mathcal{L}_j(\hat{\Delta}^{*(b)})`$.
[`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
uses two sets of draws. The pointwise variance $`\hat{\sigma}_R^2(t_j)`$
is the sample variance over the draws of the
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
fit it is given. A second, independent set (`draws_sb` draws, its own
seed) is standardized by $`\hat{\sigma}_R(t_j)`$ and gives the critical
value:

``` math
\hat{R}(t_j) \pm 1.96\, \hat{\sigma}_R(t_j) \ \text{(pointwise)}, \qquad \hat{R}(t_j) \pm c_R\, \hat{\sigma}_R(t_j), \quad c_R = \text{the } 0.95 \text{ quantile of } \max_j \frac{|\hat{R}^{*(b)}(t_j)|}{\hat{\sigma}_R(t_j)} .
```

Both sets are built from the corrected contributions when the fit and
the call use `"estimated"`; a fit made with the other option is an
error.

**Resampled score in
[`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md).**
With `draws > 0`,

``` math
U_w^{\dagger}(\hat{\beta}_w) = -\sum_{i=1}^{n} G_i\, \psi_i^{ew}(\hat{\beta}_w), \qquad U_w^{\dagger}(0) = -\sum_{i=1}^{n} G_i\, \psi_i^{ew}(0), \qquad \hat{\beta}^{\dagger} = \frac{U_w^{\dagger}(\hat{\beta}_w)}{I_w(\hat{\beta}_w)} .
```

The resampling standard error is the standard deviation of
$`\hat{\beta}^{\dagger}`$ over the draws, the de-biased estimate is
$`\hat{\beta}_w`$ minus their mean, and the resampling interval is the
de-biased estimate plus and minus the normal quantile times the
resampling standard error.

**What is and is not resampled.** $`\hat{\gamma}`$, $`IF_i`$, $`dw_i`$,
$`D_a(t)`$, $`B(\beta)`$, the curves and the information are computed
once from the data. Only the multipliers change between draws. The
estimation of the weights enters each draw through the term
$`IF_i^\top D_a(t)`$ or $`IF_i^\top B(\beta)`$ of the contribution, not
through a new fit.

**Justification.** For the unweighted estimators, replacing the
martingale increments by normal multipliers times observed increments is
justified by the results of Lin et al. (1993) and Dobler et al. (2017)
(*cited*). For the corrected contributions the argument is the
conditional multiplier central limit theorem applied to the linear
representation $`-\sum_i \xi_i(t)`$, with the estimated contributions in
place of the true ones (*heuristic argument*; not proved here). It
needs:

1.  the linear representations of section 4 to hold uniformly in $`t`$
    over the time range of the grid, which requires the weighted risk
    set of each arm to stay away from zero at the end of the range;
2.  the estimated contributions to converge to the true ones, in the
    sense that $`\sum_i \sup_t \{\xi_i(t) - \xi_i^{0}(t)\}^2 \to 0`$ in
    probability;
3.  contributions with finite second moments, which needs weights with
    finite second moments (overlap, section 7);
4.  $`\hat{\gamma}`$ to be asymptotically linear (section 3).

Under these conditions the conditional distribution of
$`\hat{\Delta}^{*}`$ given the data and the distribution of
$`\hat{\Delta} - \Delta`$ have the same Gaussian limit, and the same
holds for continuous functionals: the standardized supremum and the
integral.

**Pointwise resampling variance (`draws > 0`) under `"estimated"`.** In
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
and
[`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md),
`draws > 0` replaces the analytic pointwise variance of a curve by the
variance over draws of $`\sum_{i \in g} G_i\, w_i\, \psi_i(t)`$. In the
source that route draws multipliers for the members of one group at a
time and uses that group’s contributions only; $`\Phi_i^a(t)`$ needs the
multipliers and the influence values of all $`n`$ subjects together.
With `weight_variance = "estimated"` the argument check stops with an
error unless `draws` is 0. The development record of this version lists
`draws = 0` as a validation rule and gives no further reason. The
analytic variance $`\hat{S}^a(t)^2 \sum_i \Phi_i^a(t)^2`$ is the
quantity such a resampling variance would estimate.

## Conditions and properties

**Correct propensity model.** The estimand of section 2 requires it. The
variance formulas of section 5 describe the variability of the estimator
around its limit whether or not the model is correct; the ordering of
`"known"` and `"estimated"` in section 5 requires a correct model. The
model must be a logistic `glm` fitted on exactly the rows analysed.

**Overlap.** The fitted probabilities must stay away from 0 and 1.
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
stops when a fitted probability is numerically 0 or 1; it does not
otherwise restrict them.

**Sample size relative to the concentration of the weights.** The
linearizations are first-order. When a few subjects carry most of the
weight of an arm, the sums of section 5 are dominated by a few terms and
the normal approximation is poor. The diagnostic is the effective sample
size by arm,

``` math
\text{ESS}_a = \frac{\bigl(\sum_{i: A_i = a} w_i\bigr)^2}{\sum_{i: A_i = a} w_i^2},
```

returned with its ratio to the arm size as `df_counting()$ess`. A ratio
near 1 means evenly spread weights. The package sets no threshold on it,
and no warning depends on it.

### Simulation results by time

**Design.** Each cell is a Monte Carlo study of 2,000 replicates of one
data-generating model. A covariate $`X \sim N(0, 1)`$ drives both
treatment and outcome. The event time of a subject in arm $`a`$ is
$`\Lambda_a^{-1}(H)`$ with $`H = -\log \Phi(0.6 X + 0.8 Z)`$,
$`Z \sim N(0, 1)`$ independent of $`X`$; $`H`$ is standard exponential,
so each arm’s marginal survival is $`S_a(t) = \exp\{-\Lambda_a(t)\}`$,
with control hazard 0.04 per month and treated hazard
$`0.04 \exp(-0.4)`$ (proportional hazards, hazard ratio 0.670) or 0.04
(null). Censoring is uniform on (0, 72) months. Treatment is Bernoulli
with probability $`e(X)`$, $`\text{logit}\, e = \alpha + \beta X`$, and
the weights are built in every replicate by `ps_weights(treat ~ X)` from
the fitted model; the analyses are run under `weight_variance = "known"`
and `"estimated"` on the same data and weights. The regimes, named by
the treated arm’s ESS ratio (mean of ESS / $`n`$ over replicates,
control / treated in the tables):

- **Near 0.9**: $`\text{logit}\, e = -0.44 + 0.51 X`$, $`N = 400`$,
  stabilized weights (treated ESS ratio 0.90). The restricted Rotterdam
  example of section 8 has weights spread about this evenly. The
  coefficients were set so that the population ESS ratio of the treated
  arm is about 0.91 and the treated proportion matches that of the next
  regime.
- **About 0.69**: $`\text{logit}\, e = -0.5 + X`$, $`N = 400`$,
  stabilized weights (treated ESS ratio 0.69).
- **Strong confounding**: $`\text{logit}\, e = -0.5 + 2X`$, $`N = 400`$,
  inverse-probability weights (treated ESS ratio 0.37–0.39).
- **14% treated**: $`\text{logit}\, e = -2.3 + 1.3 X`$ (14.3% treated),
  $`N = 1{,}450`$, stabilized weights (treated ESS ratio 0.35–0.36); the
  control arm’s weights are nearly even and the treated arm’s are
  concentrated.
- **Reference, true weights**: the designs of the second and third
  regimes with the true propensities as inverse-probability weights
  (treated ESS ratio 0.69 and 0.39), analysed under `"known"` only (with
  true weights there is no estimation to carry).

Every regime has a proportional-hazards cell and a null cell, except the
reference, which has the proportional-hazards cell of each design. The
truth is in closed form: $`S_1(t) - S_0(t)`$; the RMST difference
$`\int_0^t \{S_1(u) - S_0(u)\}\, du`$; the log hazard ratio $`-0.4`$ or
0 (both arms are exponential marginally). A replicate counts at month
$`t`$ when each arm has a subject with observed time at least $`t`$. The
RMST difference at month $`t`$ is read on the grid of the
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
fit (the observed times) at the last grid time $`t_j \le t`$ of the
replicate, so its bias and coverage are against $`R(t_j)`$, the true
RMST difference at that time; the mean of $`R(t) - R(t_j)`$ over
replicates is at most 0.068 months (month 48). A Kaplan–Meier curve is
constant between grid times, so the survival difference at $`t_j`$ is
that at $`t`$; the bootstrap RMST of Table 5 integrates each curve to
month 36 itself and is scored against $`R(36)`$. Coverage is that of the
nominal 0.95 pointwise interval; an interval that could not be computed
counts as not covering. The Monte Carlo standard error (MCSE) of a
coverage near 0.95 is 0.005 at 2,000 replicates and 0.010 at 500. Where
the reported months fall among observed event times (percentile, from
the formula, marginal over $`X`$; the same in every regime because the
outcome model is the same): month 1 at percentile 6 (control) and 5
(treated, PH); month 6 at percentile 30 (control) and 26 (treated, PH);
month 12 at percentile 52 (control) and 45 (treated, PH); month 24 at
percentile 79 (control) and 72 (treated, PH); month 36 at percentile 92
(control) and 88 (treated, PH); month 48 at percentile 97 (control) and
96 (treated, PH). The 10th and 90th percentiles are months 1.8 and 33.7
(control) and 2.2 and 38.6 (treated, PH); months 1 and 48 lie outside
them.

Provenance: run for weightedsurv 0.2.0 (the near-0.9 cells at package
commit f55001f; the other cells with estimated weights at c195953 and
the true-weight cells at d9b7558); the scripts and records are kept in
the package repository under `dev/audit/`. This article reads only the
stored summary `sim_results/ps_weight_results_by_time.csv`.

Columns of Tables 1 to 3: true value; bias (mean estimate minus truth;
the estimate is the same under both options); empirical SD of the
estimate; mean standard error under each option; coverage with its MCSE
in brackets; mean interval length. ESS c / t: mean ESS ratio, control /
treated.

| Regime | Cell | Scenario | ESS c / t | Month | True | Bias | Emp SD | SE known | SE estimated | Cov known | Cov estimated | Length known | Length estimated |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Near 0.9 | Es9-400 | PH | 0.96 / 0.90 | 1 | 0.0128 | 0.0007 | 0.0182 | 0.0183 | 0.0176 | 0.945 (0.005) | 0.939 (0.005) | 0.072 | 0.069 |
|  |  |  |  | 6 | 0.0648 | 0.0001 | 0.0358 | 0.0399 | 0.0361 | 0.975 (0.004) | 0.953 (0.005) | 0.156 | 0.141 |
|  |  |  |  | 12 | 0.1061 | -0.0003 | 0.0452 | 0.0500 | 0.0443 | 0.962 (0.004) | 0.943 (0.005) | 0.196 | 0.174 |
|  |  |  |  | 24 | 0.1426 | 0.0002 | 0.0536 | 0.0578 | 0.0518 | 0.962 (0.004) | 0.934 (0.006) | 0.227 | 0.203 |
|  |  |  |  | 36 | 0.1440 | 0.0009 | 0.0563 | 0.0601 | 0.0554 | 0.965 (0.004) | 0.947 (0.005) | 0.236 | 0.217 |
|  |  |  |  | 48 | 0.1295 | 0.0008 | 0.0599 | 0.0617 | 0.0585 | 0.951 (0.005) | 0.935 (0.005) | 0.242 | 0.229 |
|  | Es9-400-null | null | 0.95 / 0.90 | 1 | 0.0000 | -0.0009 | 0.0195 | 0.0198 | 0.0190 | 0.952 (0.005) | 0.945 (0.005) | 0.078 | 0.074 |
|  |  |  |  | 6 | 0.0000 | -0.0001 | 0.0386 | 0.0434 | 0.0387 | 0.977 (0.003) | 0.952 (0.005) | 0.170 | 0.152 |
|  |  |  |  | 12 | 0.0000 | 0.0003 | 0.0477 | 0.0535 | 0.0469 | 0.975 (0.004) | 0.949 (0.005) | 0.210 | 0.184 |
|  |  |  |  | 24 | 0.0000 | 0.0003 | 0.0543 | 0.0586 | 0.0526 | 0.968 (0.004) | 0.951 (0.005) | 0.230 | 0.206 |
|  |  |  |  | 36 | 0.0000 | -0.0011 | 0.0542 | 0.0571 | 0.0532 | 0.959 (0.004) | 0.943 (0.005) | 0.224 | 0.208 |
|  |  |  |  | 48 | 0.0000 | -0.0021 | 0.0539 | 0.0543 | 0.0520 | 0.944 (0.005) | 0.929 (0.006) | 0.213 | 0.204 |
| About 0.69 | Es-400 | PH | 0.83 / 0.69 | 1 | 0.0128 | -0.0001 | 0.0215 | 0.0211 | 0.0202 | 0.928 (0.006) | 0.923 (0.006) | 0.083 | 0.079 |
|  |  |  |  | 6 | 0.0648 | -0.0001 | 0.0389 | 0.0439 | 0.0387 | 0.970 (0.004) | 0.948 (0.005) | 0.172 | 0.152 |
|  |  |  |  | 12 | 0.1061 | 0.0004 | 0.0464 | 0.0543 | 0.0467 | 0.978 (0.003) | 0.954 (0.005) | 0.213 | 0.183 |
|  |  |  |  | 24 | 0.1426 | -0.0010 | 0.0569 | 0.0636 | 0.0559 | 0.969 (0.004) | 0.945 (0.005) | 0.249 | 0.219 |
|  |  |  |  | 36 | 0.1440 | -0.0033 | 0.0640 | 0.0678 | 0.0620 | 0.953 (0.005) | 0.934 (0.006) | 0.266 | 0.243 |
|  |  |  |  | 48 | 0.1295 | -0.0058 | 0.0692 | 0.0703 | 0.0664 | 0.937 (0.005) | 0.928 (0.006) | 0.275 | 0.260 |
|  | Es-400-null | null |  | 1 | 0.0000 | -0.0000 | 0.0225 | 0.0223 | 0.0212 | 0.940 (0.005) | 0.931 (0.006) | 0.087 | 0.083 |
|  |  |  |  | 6 | 0.0000 | -0.0022 | 0.0425 | 0.0476 | 0.0415 | 0.973 (0.004) | 0.946 (0.005) | 0.186 | 0.163 |
|  |  |  |  | 12 | 0.0000 | -0.0017 | 0.0514 | 0.0590 | 0.0506 | 0.975 (0.003) | 0.947 (0.005) | 0.231 | 0.198 |
|  |  |  |  | 24 | 0.0000 | -0.0015 | 0.0606 | 0.0662 | 0.0589 | 0.970 (0.004) | 0.946 (0.005) | 0.259 | 0.231 |
|  |  |  |  | 36 | 0.0000 | -0.0038 | 0.0654 | 0.0655 | 0.0608 | 0.946 (0.005) | 0.929 (0.006) | 0.257 | 0.238 |
|  |  |  |  | 48 | 0.0000 | -0.0042 | 0.0654 | 0.0619 | 0.0593 | 0.912 (0.006) | 0.900 (0.007) | 0.243 | 0.232 |
| Strong confounding, 0.33-0.39 | Fe-400 | PH | 0.50 / 0.39 | 1 | 0.0128 | -0.0015 | 0.0461 | 0.0282 | 0.0268 | 0.741 (0.010) | 0.737 (0.010) | 0.111 | 0.105 |
|  |  |  |  | 6 | 0.0648 | -0.0079 | 0.0648 | 0.0615 | 0.0542 | 0.932 (0.006) | 0.901 (0.007) | 0.241 | 0.213 |
|  |  |  |  | 12 | 0.1061 | -0.0091 | 0.0687 | 0.0733 | 0.0622 | 0.971 (0.004) | 0.939 (0.005) | 0.287 | 0.244 |
|  |  |  |  | 24 | 0.1426 | -0.0129 | 0.0807 | 0.0832 | 0.0722 | 0.955 (0.005) | 0.927 (0.006) | 0.326 | 0.283 |
|  |  |  |  | 36 | 0.1440 | -0.0157 | 0.0922 | 0.0869 | 0.0789 | 0.913 (0.006) | 0.893 (0.007) | 0.341 | 0.309 |
|  |  |  |  | 48 | 0.1295 | -0.0202 | 0.0986 | 0.0863 | 0.0811 | 0.868 (0.008) | 0.854 (0.008) | 0.338 | 0.318 |
|  | Fe-400-null | null | 0.51 / 0.37 | 1 | 0.0000 | -0.0039 | 0.0448 | 0.0283 | 0.0267 | 0.780 (0.009) | 0.765 (0.009) | 0.111 | 0.105 |
|  |  |  |  | 6 | 0.0000 | -0.0070 | 0.0685 | 0.0654 | 0.0567 | 0.932 (0.006) | 0.899 (0.007) | 0.256 | 0.222 |
|  |  |  |  | 12 | 0.0000 | -0.0089 | 0.0811 | 0.0808 | 0.0683 | 0.943 (0.005) | 0.901 (0.007) | 0.317 | 0.268 |
|  |  |  |  | 24 | 0.0000 | -0.0097 | 0.0949 | 0.0898 | 0.0792 | 0.932 (0.006) | 0.898 (0.007) | 0.352 | 0.310 |
|  |  |  |  | 36 | 0.0000 | -0.0128 | 0.1023 | 0.0850 | 0.0784 | 0.865 (0.008) | 0.842 (0.008) | 0.333 | 0.307 |
|  |  |  |  | 48 | 0.0000 | -0.0081 | 0.0974 | 0.0743 | 0.0704 | 0.822 (0.009) | 0.808 (0.009) | 0.291 | 0.276 |
| 14% treated, 0.35 | V-1450-ph | PH | 0.95 / 0.35 | 1 | 0.0128 | -0.0001 | 0.0111 | 0.0111 | 0.0104 | 0.968 (0.004) | 0.950 (0.005) | 0.044 | 0.041 |
|  |  |  |  | 6 | 0.0648 | -0.0005 | 0.0287 | 0.0312 | 0.0276 | 0.971 (0.004) | 0.935 (0.006) | 0.122 | 0.108 |
|  |  |  |  | 12 | 0.1061 | -0.0017 | 0.0439 | 0.0466 | 0.0420 | 0.962 (0.004) | 0.935 (0.006) | 0.183 | 0.165 |
|  |  |  |  | 24 | 0.1426 | -0.0054 | 0.0652 | 0.0656 | 0.0614 | 0.946 (0.005) | 0.926 (0.006) | 0.257 | 0.241 |
|  |  |  |  | 36 | 0.1440 | -0.0079 | 0.0779 | 0.0754 | 0.0724 | 0.920 (0.006) | 0.911 (0.006) | 0.296 | 0.284 |
|  |  |  |  | 48 | 0.1295 | -0.0100 | 0.0902 | 0.0786 | 0.0764 | 0.866 (0.008) | 0.860 (0.008) | 0.308 | 0.299 |
|  | V-1450-null | null | 0.95 / 0.36 | 1 | 0.0000 | -0.0002 | 0.0126 | 0.0132 | 0.0121 | 0.963 (0.004) | 0.943 (0.005) | 0.052 | 0.047 |
|  |  |  |  | 6 | 0.0000 | -0.0015 | 0.0358 | 0.0385 | 0.0343 | 0.966 (0.004) | 0.930 (0.006) | 0.151 | 0.134 |
|  |  |  |  | 12 | 0.0000 | -0.0037 | 0.0547 | 0.0558 | 0.0509 | 0.949 (0.005) | 0.920 (0.006) | 0.219 | 0.200 |
|  |  |  |  | 24 | 0.0000 | -0.0054 | 0.0743 | 0.0712 | 0.0676 | 0.927 (0.006) | 0.915 (0.006) | 0.279 | 0.265 |
|  |  |  |  | 36 | 0.0000 | -0.0080 | 0.0809 | 0.0724 | 0.0702 | 0.876 (0.007) | 0.871 (0.008) | 0.284 | 0.275 |
|  |  |  |  | 48 | 0.0000 | -0.0060 | 0.0785 | 0.0661 | 0.0647 | 0.822 (0.009) | 0.817 (0.009) | 0.259 | 0.254 |
| Reference: true weights | D-400 | PH | 0.84 / 0.69 | 1 | 0.0128 | 0.0000 | 0.0217 | 0.0207 |  | 0.917 (0.006) |  | 0.081 |  |
|  |  |  |  | 6 | 0.0648 | -0.0006 | 0.0434 | 0.0433 |  | 0.943 (0.005) |  | 0.170 |  |
|  |  |  |  | 12 | 0.1061 | -0.0010 | 0.0531 | 0.0538 |  | 0.951 (0.005) |  | 0.211 |  |
|  |  |  |  | 24 | 0.1426 | -0.0013 | 0.0623 | 0.0632 |  | 0.951 (0.005) |  | 0.248 |  |
|  |  |  |  | 36 | 0.1440 | -0.0018 | 0.0678 | 0.0674 |  | 0.938 (0.005) |  | 0.264 |  |
|  |  |  |  | 48 | 0.1295 | -0.0017 | 0.0729 | 0.0703 |  | 0.932 (0.006) |  | 0.275 |  |
|  | F-400 |  | 0.51 / 0.39 | 1 | 0.0128 | -0.0034 | 0.0410 | 0.0269 |  | 0.740 (0.010) |  | 0.106 |  |
|  |  |  |  | 6 | 0.0648 | -0.0057 | 0.0678 | 0.0605 |  | 0.915 (0.006) |  | 0.237 |  |
|  |  |  |  | 12 | 0.1061 | -0.0072 | 0.0777 | 0.0722 |  | 0.936 (0.005) |  | 0.283 |  |
|  |  |  |  | 24 | 0.1426 | -0.0130 | 0.0923 | 0.0822 |  | 0.916 (0.006) |  | 0.322 |  |
|  |  |  |  | 36 | 0.1440 | -0.0159 | 0.0997 | 0.0863 |  | 0.887 (0.007) |  | 0.338 |  |
|  |  |  |  | 48 | 0.1295 | -0.0189 | 0.1046 | 0.0861 |  | 0.858 (0.008) |  | 0.338 |  |

Table 1. Survival difference $`S_1(t) - S_0(t)`$. {.table}

Evaluable replicates are fewer than 2,000 only at months 36, 42 and 48,
in 6 cells (fewest: 1,919, Fe-400-null, month 48).

![Figure 1. Coverage of the survival difference at the 11 stored months.
Solid: "known"; dashed: "estimated" (absent with true weights). Grey
band: 0.95 plus or minus 2 MCSE at 2,000
replicates.](estimated_weights_methodology_files/figure-html/simres-fig1-1.png)

Figure 1. Coverage of the survival difference at the 11 stored months.
Solid: “known”; dashed: “estimated” (absent with true weights). Grey
band: 0.95 plus or minus 2 MCSE at 2,000 replicates.

| Regime | Cell | Scenario | ESS c / t | Month | True | Bias | Emp SD | SE known | SE estimated | Cov known | Cov estimated | Length known | Length estimated |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Near 0.9 | Es9-400 | PH | 0.96 / 0.90 | 1 | 0.006 | 0.000 | 0.010 | 0.010 | 0.010 | 0.948 (0.005) | 0.942 (0.005) | 0.04 | 0.04 |
|  |  |  |  | 6 | 0.208 | 0.003 | 0.130 | 0.143 | 0.130 | 0.973 (0.004) | 0.954 (0.005) | 0.56 | 0.51 |
|  |  |  |  | 12 | 0.730 | 0.002 | 0.328 | 0.372 | 0.323 | 0.975 (0.003) | 0.953 (0.005) | 1.46 | 1.27 |
|  |  |  |  | 24 | 2.271 | -0.003 | 0.782 | 0.899 | 0.758 | 0.974 (0.004) | 0.941 (0.005) | 3.52 | 2.97 |
|  |  |  |  | 36 | 4.014 | 0.006 | 1.236 | 1.434 | 1.207 | 0.976 (0.003) | 0.944 (0.005) | 5.62 | 4.73 |
|  |  |  |  | 48 | 5.664 | 0.019 | 1.685 | 1.948 | 1.656 | 0.975 (0.003) | 0.945 (0.005) | 7.64 | 6.49 |
|  | Es9-400-null | null | 0.95 / 0.90 | 1 | 0.000 | -0.000 | 0.011 | 0.011 | 0.011 | 0.949 (0.005) | 0.943 (0.005) | 0.04 | 0.04 |
|  |  |  |  | 6 | 0.000 | -0.003 | 0.140 | 0.156 | 0.139 | 0.970 (0.004) | 0.948 (0.005) | 0.61 | 0.54 |
|  |  |  |  | 12 | 0.000 | -0.003 | 0.346 | 0.401 | 0.342 | 0.975 (0.004) | 0.950 (0.005) | 1.57 | 1.34 |
|  |  |  |  | 24 | 0.000 | 0.011 | 0.798 | 0.944 | 0.786 | 0.984 (0.003) | 0.952 (0.005) | 3.70 | 3.08 |
|  |  |  |  | 36 | 0.000 | 0.003 | 1.246 | 1.456 | 1.219 | 0.978 (0.003) | 0.945 (0.005) | 5.71 | 4.78 |
|  |  |  |  | 48 | 0.000 | -0.022 | 1.660 | 1.902 | 1.614 | 0.976 (0.003) | 0.941 (0.005) | 7.45 | 6.33 |
| About 0.69 | Es-400 | PH | 0.83 / 0.69 | 1 | 0.006 | 0.000 | 0.013 | 0.012 | 0.011 | 0.893 (0.007) | 0.884 (0.007) | 0.05 | 0.04 |
|  |  |  |  | 6 | 0.208 | 0.003 | 0.150 | 0.162 | 0.144 | 0.961 (0.004) | 0.939 (0.005) | 0.64 | 0.57 |
|  |  |  |  | 12 | 0.730 | 0.005 | 0.356 | 0.414 | 0.349 | 0.978 (0.003) | 0.954 (0.005) | 1.62 | 1.37 |
|  |  |  |  | 24 | 2.271 | 0.005 | 0.817 | 0.993 | 0.811 | 0.980 (0.003) | 0.954 (0.005) | 3.89 | 3.18 |
|  |  |  |  | 36 | 4.014 | -0.012 | 1.330 | 1.601 | 1.314 | 0.980 (0.003) | 0.947 (0.005) | 6.27 | 5.15 |
|  |  |  |  | 48 | 5.664 | -0.064 | 1.866 | 2.193 | 1.833 | 0.976 (0.003) | 0.942 (0.005) | 8.60 | 7.19 |
|  | Es-400-null | null |  | 1 | 0.000 | 0.000 | 0.013 | 0.012 | 0.012 | 0.919 (0.006) | 0.911 (0.006) | 0.05 | 0.05 |
|  |  |  |  | 6 | 0.000 | -0.008 | 0.158 | 0.173 | 0.152 | 0.959 (0.004) | 0.933 (0.006) | 0.68 | 0.59 |
|  |  |  |  | 12 | 0.000 | -0.023 | 0.381 | 0.445 | 0.371 | 0.977 (0.003) | 0.937 (0.005) | 1.74 | 1.45 |
|  |  |  |  | 24 | 0.000 | -0.044 | 0.875 | 1.060 | 0.864 | 0.979 (0.003) | 0.948 (0.005) | 4.16 | 3.39 |
|  |  |  |  | 36 | 0.000 | -0.076 | 1.419 | 1.657 | 1.368 | 0.976 (0.003) | 0.945 (0.005) | 6.50 | 5.36 |
|  |  |  |  | 48 | 0.000 | -0.110 | 1.947 | 2.184 | 1.841 | 0.969 (0.004) | 0.939 (0.005) | 8.56 | 7.22 |
| Strong confounding, 0.33-0.39 | Fe-400 | PH | 0.50 / 0.39 | 1 | 0.006 | -0.001 | 0.029 | 0.015 | 0.014 | 0.686 (0.010) | 0.678 (0.010) | 0.06 | 0.05 |
|  |  |  |  | 6 | 0.208 | -0.021 | 0.286 | 0.233 | 0.208 | 0.892 (0.007) | 0.863 (0.008) | 0.91 | 0.82 |
|  |  |  |  | 12 | 0.730 | -0.072 | 0.619 | 0.587 | 0.500 | 0.945 (0.005) | 0.910 (0.006) | 2.30 | 1.96 |
|  |  |  |  | 24 | 2.271 | -0.202 | 1.291 | 1.375 | 1.125 | 0.961 (0.004) | 0.919 (0.006) | 5.39 | 4.41 |
|  |  |  |  | 36 | 4.014 | -0.366 | 2.041 | 2.184 | 1.796 | 0.957 (0.005) | 0.915 (0.006) | 8.56 | 7.04 |
|  |  |  |  | 48 | 5.664 | -0.575 | 2.840 | 2.945 | 2.466 | 0.938 (0.005) | 0.908 (0.006) | 11.54 | 9.67 |
|  | Fe-400-null | null | 0.51 / 0.37 | 1 | 0.000 | -0.002 | 0.031 | 0.015 | 0.014 | 0.736 (0.010) | 0.716 (0.010) | 0.06 | 0.06 |
|  |  |  |  | 6 | 0.000 | -0.030 | 0.287 | 0.240 | 0.210 | 0.905 (0.007) | 0.864 (0.008) | 0.94 | 0.83 |
|  |  |  |  | 12 | 0.000 | -0.079 | 0.660 | 0.623 | 0.520 | 0.930 (0.006) | 0.892 (0.007) | 2.44 | 2.04 |
|  |  |  |  | 24 | 0.000 | -0.195 | 1.498 | 1.487 | 1.217 | 0.936 (0.005) | 0.889 (0.007) | 5.83 | 4.77 |
|  |  |  |  | 36 | 0.000 | -0.320 | 2.407 | 2.310 | 1.917 | 0.926 (0.006) | 0.883 (0.007) | 9.06 | 7.51 |
|  |  |  |  | 48 | 0.000 | -0.394 | 3.249 | 2.996 | 2.531 | 0.910 (0.007) | 0.873 (0.008) | 11.74 | 9.92 |
| 14% treated, 0.35 | V-1450-ph | PH | 0.95 / 0.35 | 1 | 0.006 | -0.000 | 0.006 | 0.006 | 0.006 | 0.961 (0.004) | 0.943 (0.005) | 0.02 | 0.02 |
|  |  |  |  | 6 | 0.208 | -0.000 | 0.091 | 0.103 | 0.089 | 0.973 (0.004) | 0.942 (0.005) | 0.40 | 0.35 |
|  |  |  |  | 12 | 0.730 | -0.006 | 0.275 | 0.306 | 0.264 | 0.970 (0.004) | 0.941 (0.005) | 1.20 | 1.03 |
|  |  |  |  | 24 | 2.271 | -0.047 | 0.837 | 0.888 | 0.784 | 0.961 (0.004) | 0.930 (0.006) | 3.48 | 3.07 |
|  |  |  |  | 36 | 4.014 | -0.129 | 1.543 | 1.586 | 1.435 | 0.946 (0.005) | 0.926 (0.006) | 6.22 | 5.63 |
|  |  |  |  | 48 | 5.664 | -0.233 | 2.320 | 2.306 | 2.123 | 0.931 (0.006) | 0.913 (0.006) | 9.04 | 8.32 |
|  | V-1450-null | null | 0.95 / 0.36 | 1 | 0.000 | -0.000 | 0.007 | 0.007 | 0.007 | 0.968 (0.004) | 0.954 (0.005) | 0.03 | 0.03 |
|  |  |  |  | 6 | 0.000 | -0.003 | 0.111 | 0.126 | 0.108 | 0.975 (0.004) | 0.939 (0.005) | 0.49 | 0.43 |
|  |  |  |  | 12 | 0.000 | -0.019 | 0.341 | 0.375 | 0.325 | 0.964 (0.004) | 0.938 (0.005) | 1.47 | 1.27 |
|  |  |  |  | 24 | 0.000 | -0.073 | 1.004 | 1.031 | 0.923 | 0.949 (0.005) | 0.922 (0.006) | 4.04 | 3.62 |
|  |  |  |  | 36 | 0.000 | -0.158 | 1.753 | 1.718 | 1.574 | 0.931 (0.006) | 0.910 (0.006) | 6.74 | 6.17 |
|  |  |  |  | 48 | 0.000 | -0.173 | 2.457 | 2.340 | 2.175 | 0.914 (0.006) | 0.899 (0.007) | 9.17 | 8.53 |
| Reference: true weights | D-400 | PH | 0.84 / 0.69 | 1 | 0.006 | 0.000 | 0.013 | 0.011 |  | 0.896 (0.007) |  | 0.04 |  |
|  |  |  |  | 6 | 0.208 | -0.002 | 0.162 | 0.159 |  | 0.940 (0.005) |  | 0.62 |  |
|  |  |  |  | 12 | 0.730 | -0.007 | 0.406 | 0.407 |  | 0.947 (0.005) |  | 1.60 |  |
|  |  |  |  | 24 | 2.271 | -0.020 | 0.969 | 0.982 |  | 0.949 (0.005) |  | 3.85 |  |
|  |  |  |  | 36 | 4.014 | -0.038 | 1.558 | 1.586 |  | 0.952 (0.005) |  | 6.22 |  |
|  |  |  |  | 48 | 5.664 | -0.062 | 2.151 | 2.181 |  | 0.948 (0.005) |  | 8.55 |  |
|  | F-400 |  | 0.51 / 0.39 | 1 | 0.006 | -0.001 | 0.023 | 0.014 |  | 0.695 (0.010) |  | 0.06 |  |
|  |  |  |  | 6 | 0.208 | -0.022 | 0.279 | 0.226 |  | 0.876 (0.007) |  | 0.89 |  |
|  |  |  |  | 12 | 0.730 | -0.061 | 0.657 | 0.573 |  | 0.910 (0.006) |  | 2.25 |  |
|  |  |  |  | 24 | 2.271 | -0.187 | 1.507 | 1.348 |  | 0.919 (0.006) |  | 5.28 |  |
|  |  |  |  | 36 | 4.014 | -0.363 | 2.419 | 2.142 |  | 0.905 (0.007) |  | 8.40 |  |
|  |  |  |  | 48 | 5.664 | -0.552 | 3.321 | 2.903 |  | 0.898 (0.007) |  | 11.38 |  |

Table 2. RMST difference $`\int_0^t \{S_1(u) - S_0(u)\}\, du`$ (months).
{.table style="width:100%;"}

Table 2: True is $`R(t)`$ at the stated month; bias and coverage are
against $`R(t_j)`$, $`t_j`$ the last grid time $`\le t`$ of each
replicate (design paragraph above).

| Regime | Cell | Scenario | ESS c / t | Statistic | True | Bias | Emp SD | SE known | SE estimated | Cov known | Cov estimated | Length known | Length estimated |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Near 0.9 | Es9-400 | PH | 0.96 / 0.90 | FH(0,0) | -0.4 | -0.0030 | 0.1203 | 0.1387 | 0.1188 | 0.976 (0.003) | 0.953 (0.005) | 0.544 | 0.466 |
|  |  |  |  | FH(0,1) | -0.4 | -0.0056 | 0.1598 | 0.1659 | 0.1561 | 0.954 (0.005) | 0.945 (0.005) | 0.650 | 0.612 |
|  | Es9-400-null | null | 0.95 / 0.90 | FH(0,0) | 0.0 | 0.0037 | 0.1167 | 0.1332 | 0.1132 | 0.977 (0.003) | 0.943 (0.005) | 0.522 | 0.444 |
|  |  |  |  | FH(0,1) | 0.0 | 0.0089 | 0.1569 | 0.1599 | 0.1505 | 0.951 (0.005) | 0.939 (0.005) | 0.627 | 0.590 |
| About 0.69 | Es-400 | PH | 0.83 / 0.69 | FH(0,0) | -0.4 | 0.0059 | 0.1356 | 0.1565 | 0.1323 | 0.976 (0.003) | 0.944 (0.005) | 0.613 | 0.519 |
|  |  |  |  | FH(0,1) | -0.4 | 0.0120 | 0.1816 | 0.1868 | 0.1751 | 0.952 (0.005) | 0.936 (0.005) | 0.732 | 0.686 |
|  | Es-400-null | null |  | FH(0,0) | 0.0 | 0.0115 | 0.1399 | 0.1537 | 0.1301 | 0.968 (0.004) | 0.943 (0.005) | 0.603 | 0.510 |
|  |  |  |  | FH(0,1) | 0.0 | 0.0170 | 0.1921 | 0.1843 | 0.1734 | 0.940 (0.005) | 0.925 (0.006) | 0.723 | 0.680 |
| Strong confounding, 0.33-0.39 | Fe-400 | PH | 0.50 / 0.39 | FH(0,0) | -0.4 | 0.0415 | 0.2086 | 0.2099 | 0.1768 | 0.933 (0.006) | 0.901 (0.007) | 0.823 | 0.693 |
|  |  |  |  | FH(0,1) | -0.4 | 0.0482 | 0.2655 | 0.2341 | 0.2173 | 0.898 (0.007) | 0.880 (0.007) | 0.918 | 0.852 |
|  | Fe-400-null | null | 0.51 / 0.37 | FH(0,0) | 0.0 | 0.0351 | 0.2383 | 0.2138 | 0.1811 | 0.914 (0.006) | 0.871 (0.008) | 0.838 | 0.710 |
|  |  |  |  | FH(0,1) | 0.0 | 0.0415 | 0.2992 | 0.2365 | 0.2203 | 0.864 (0.008) | 0.847 (0.008) | 0.927 | 0.864 |
| 14% treated, 0.35 | V-1450-ph | PH | 0.95 / 0.35 | FH(0,0) | -0.4 | 0.0142 | 0.1780 | 0.1712 | 0.1588 | 0.920 (0.006) | 0.902 (0.007) | 0.671 | 0.623 |
|  | V-1450-null | null | 0.95 / 0.36 | FH(0,0) | 0.0 | 0.0181 | 0.1793 | 0.1658 | 0.1544 | 0.906 (0.007) | 0.889 (0.007) | 0.650 | 0.605 |
| Reference: true weights | D-400 | PH | 0.84 / 0.69 | FH(0,0) | -0.4 | 0.0020 | 0.1549 | 0.1560 |  | 0.952 (0.005) |  | 0.611 |  |
|  |  |  |  | FH(0,1) | -0.4 | 0.0004 | 0.1911 | 0.1865 |  | 0.945 (0.005) |  | 0.731 |  |
|  | F-400 |  | 0.51 / 0.39 | FH(0,0) | -0.4 | 0.0380 | 0.2414 | 0.2080 |  | 0.899 (0.007) |  | 0.815 |  |
|  |  |  |  | FH(0,1) | -0.4 | 0.0441 | 0.2823 | 0.2340 |  | 0.875 (0.007) |  | 0.917 |  |

Table 3. Log hazard ratio, weighted Cox estimate under FH(0,0) and
FH(0,1) (FH(0,1) not run in the 14%-treated cells). {.table}

| Regime | Cell | Scenario | Band | Cov known | Cov estimated | Width known | Width estimated |
|:---|:---|:---|:---|---:|---:|---:|---:|
| Near 0.9 | Es9-400 | PH | difference | 0.978 (0.003) | 0.955 (0.005) | 0.313 | 0.286 |
|  |  |  | RMST | 0.987 (0.003) | 0.971 (0.004) | 3.234 | 2.787 |
|  | Es9-400-null | null | difference | 0.983 (0.003) | 0.956 (0.005) | 0.324 | 0.294 |
|  |  |  | RMST | 0.987 (0.003) | 0.967 (0.004) | 3.223 | 2.750 |
| About 0.69 | Es-400 | PH | difference | 0.976 (0.003) | 0.949 (0.005) | 0.342 | 0.306 |
|  |  |  | RMST | 0.982 (0.003) | 0.958 (0.005) | 3.597 | 3.020 |
|  | Es-400-null | null | difference | 0.962 (0.004) | 0.939 (0.005) | 0.358 | 0.320 |
|  |  |  | RMST | 0.980 (0.003) | 0.958 (0.005) | 3.622 | 3.036 |
| Strong confounding, 0.33-0.39 | Fe-400 | PH | difference | 0.897 (0.007) | 0.864 (0.008) | 0.446 | 0.394 |
|  |  |  | RMST | 0.887 (0.007) | 0.864 (0.008) | 4.981 | 4.181 |
|  | Fe-400-null | null | difference | 0.865 (0.008) | 0.818 (0.009) | 0.475 | 0.420 |
|  |  |  | RMST | 0.910 (0.006) | 0.872 (0.007) | 5.038 | 4.226 |
| 14% treated, 0.35 | V-1450-ph | PH | difference | 0.932 (0.006) | 0.897 (0.007) | 0.313 | 0.292 |
|  |  |  | RMST | 0.973 (0.004) | 0.949 (0.005) | 3.101 | 2.774 |
|  | V-1450-null | null | difference | 0.875 (0.007) | 0.840 (0.008) | 0.350 | 0.329 |
|  |  |  | RMST | 0.962 (0.004) | 0.935 (0.006) | 3.486 | 3.150 |
| Reference: true weights | D-400 | PH | difference | 0.954 (0.005) |  | 0.339 |  |
|  |  |  | RMST | 0.966 (0.004) |  | 3.540 |  |
|  | F-400 |  | difference | 0.870 (0.008) |  | 0.441 |  |
|  |  |  | RMST | 0.873 (0.007) |  | 4.883 |  |

Table 4. Simultaneous 0.95 bands on \[3, 42\] months: coverage (the true
curve inside the band at every grid point of \[3, 42\] where the band
has positive width) and mean width over those points. {.table}

| Regime | Cell | Estimand | Cov known | Cov estimated | Cov boot normal | Cov boot percentile | Length known | Length estimated | Length boot normal | Length boot percentile |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| About 0.69 | E-400 | difference, month 18 | 0.972 | 0.934 | 0.936 (0.011) | 0.936 | 0.235 | 0.202 | 0.204 | 0.200 |
|  |  | RMST, month 36 | 0.974 | 0.940 | 0.942 (0.010) | 0.936 | 6.267 | 5.146 | 5.258 | 5.142 |
|  |  | log HR, FH(0,0) | 0.980 | 0.954 | 0.960 (0.009) | 0.956 | 0.613 | 0.517 | 0.531 | 0.519 |
| 14% treated, 0.35 | V-1450-ph | difference, month 18 | 0.948 | 0.932 | 0.932 (0.011) | 0.932 | 0.227 | 0.209 | 0.212 | 0.207 |
|  |  | RMST, month 36 | 0.938 | 0.922 | 0.920 (0.012) | 0.922 | 6.229 | 5.653 | 5.760 | 5.598 |
|  |  | log HR, FH(0,0) | 0.914 | 0.890 | 0.892 (0.014) | 0.894 | 0.669 | 0.623 | 0.635 | 0.616 |
|  | V-1450-null | difference, month 18 | 0.930 | 0.912 | 0.906 (0.013) | 0.912 | 0.256 | 0.240 | 0.242 | 0.236 |
|  |  | RMST, month 36 | 0.926 | 0.900 | 0.898 (0.014) | 0.906 | 6.677 | 6.122 | 6.203 | 6.034 |
|  |  | log HR, FH(0,0) | 0.898 | 0.884 | 0.882 (0.014) | 0.876 | 0.646 | 0.602 | 0.613 | 0.596 |

Table 5. Refitting bootstrap on the first 500 replicates of three cells,
beside “known” and “estimated” on the same replicates. {.table}

In Table 5 the bootstrap resamples subjects, refits the propensity model
and rebuilds the weights in every resample (200 resamples per
replicate); “normal” is the interval from the bootstrap standard
deviation, “percentile” the 2.5% and 97.5% quantiles. Cell E-400 has the
data of Es-400 with inverse-probability weights from the same fitted
model; its `"estimated"` column was computed by the evaluation’s own
script, which equals the package option on the replicates both ran.

**Reading.**

*Near 0.9 (treated ESS ratio 0.90).* The survival difference is covered
at 0.929–0.956 under `"estimated"` and 0.944–0.977 under `"known"`
across all 11 months (Figure 1). The `"estimated"` standard error is
0.96–1.01 times the empirical SD, the `"known"` one 1.01–1.12 times. The
bias is at most 0.0023 in absolute value at every month. RMST
difference: 0.939–0.954 (`"estimated"`) and 0.948–0.984 (`"known"`). Log
hazard ratio under FH(0,0): 0.943–0.953 and 0.976–0.977. Bands:
0.955–0.956 and 0.978–0.983 for the difference, 0.967–0.971 and 0.987
for RMST. Interval length under `"estimated"` is 0.88–0.96 of `"known"`
for the difference.

*About 0.69 (treated ESS ratio 0.69).* At months 3 to 30 the difference
is covered at 0.930–0.954 (`"estimated"`) and 0.955–0.978 (`"known"`).
At month 1: 0.923–0.931 and 0.928–0.940; at month 48: 0.900–0.928 and
0.912–0.937. The bias at month 48 is -0.0058 (PH) and -0.0042 (null),
against at most 0.0030 in absolute value at months 3 to 30.

*Strong confounding (treated ESS ratio 0.37–0.39).* Both options fall
short at month 1 (0.737–0.765 `"estimated"`, 0.741–0.780 `"known"`) and
at month 48 (0.808–0.854 and 0.822–0.868); at months 6 to 24 they cover
0.898–0.939 and 0.932–0.971. `"estimated"` is below `"known"` at every
month, by 0.004–0.042. The bias is negative and grows with time: -0.0015
at month 1 to -0.0202 at month 48 in Fe-400 (9.2 MCSE from zero). For
the RMST difference, against $`R(t_j)`$, it accumulates: -0.37 months at
month 36 and -0.57 at month 48 (8.0 and 9.1 MCSE from zero; true values
4.01 and 5.66).

*14% treated (treated ESS ratio 0.35–0.36).* Coverage of the difference
falls with time under both options: 0.943–0.950 (`"estimated"`) and
0.963–0.968 (`"known"`) at month 1, 0.915–0.926 and 0.927–0.946 at month
24, 0.817–0.860 and 0.822–0.866 at month 48. The two options are
0.018–0.020 apart at month 1 and 0.005–0.006 at month 48. The bias
reaches -0.0100 at month 48 (PH).

*True weights against estimated weights in the same design.* With strong
confounding, `"estimated"` on fitted weights (Fe-400) and `"known"` on
true weights (F-400) cover the difference at 0.737 and 0.740 at month 1,
0.927 and 0.916 at month 24, 0.854 and 0.858 at month 48. In the design
of the About-0.69 regime (Es-400 against D-400): 0.923 and 0.917 at
month 1, 0.945 and 0.951 at month 24, 0.928 and 0.932 at month 48.
`"known"` on the fitted weights of the same data covers 0.969 at month
24.

*Log hazard ratio and bands (Tables 3 and 4).* Under FH(0,0),
`"estimated"` covers 0.943–0.953 in the near-0.9 and 0.69 regimes and
0.871–0.902 with concentrated weights; `"known"` 0.968–0.977 and
0.906–0.933. The band for the difference: `"estimated"` 0.939–0.956 and
0.818–0.897, `"known"` 0.962–0.983 and 0.865–0.932.

*Bootstrap (Table 5).* On the same 500 replicates, the coverage under
`"estimated"` and that of the refitting bootstrap (normal interval)
differ by at most 0.006, and the `"estimated"` interval is 0.97–0.99
times the bootstrap’s length. `"known"` covers 0.016–0.036 more than the
bootstrap.

Which option is nearer the nominal level differs by regime and by month
(Figure 1).

## Numerical illustration

The restricted Rotterdam population of the `weightedsurv_examples`
vignette: node-positive, postmenopausal women without chemotherapy, with
stabilized weights from the model
`treat ~ age + size + grade3 + nodes + pgr + er`. The `"estimated"`
standard error is rebuilt from the formulas of sections 3 to 5 and shown
beside the package value, for the survival difference at one time point
and for the log hazard ratio under FH(0, 0). Only exported functions and
the objects they return are used.

``` r

library(survival)
library(weightedsurv)

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

df_rott <- within(rotterdam, {
  rfstime <- ifelse(recur == 1, rtime, dtime)
  tte <- rfstime / 30.4375
  ignore <- (recur == 0 & death == 1 & rtime < dtime)
  event <- ifelse(recur == 1 | ignore, recur, death)
  grade3 <- ifelse(grade == "3", 1, 0)
  treat <- hormon
})
df_rott <- subset(df_rott, nodes > 0 & meno == 1 & chemo == 0)

pw <- ps_weights(treat ~ age + size + grade3 + nodes + pgr + er,
                 data = df_rott, type = "stabilized")
df_rott$sw.weights <- pw$weights

dc_est <- df_counting(df_rott, "tte", "event", "treat", weight.name = "sw.weights",
                      weight_variance = "estimated", weight_model = pw)
dc_known <- df_counting(df_rott, "tte", "event", "treat", weight.name = "sw.weights")
```

The population has 829 women, 291 treated. Effective sample size by arm
(`df_counting()$ess`):

``` r

knitr::kable(dc_est$ess, digits = c(0, 0, 1, 3))
```

| group   |   n |   ess | ratio |
|:--------|----:|------:|------:|
| treat   | 291 | 270.7 | 0.930 |
| control | 538 | 529.3 | 0.984 |

The pieces.
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
returns the time-fixed data sorted by time (`time`, `delta`, `z`,
`w_hat`) and, under `"estimated"`, `ew_parts`: the matrices `IF` and
`dw` of the
[`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md)
object with their rows in that order.

``` r

U <- dc_est$time; D <- dc_est$delta; A <- dc_est$z; w <- dc_est$w_hat
IF <- dc_est$ew_parts$IF; dw <- dc_est$ew_parts$dw

# the same matrices as pw$IF and pw$dw, rows reordered
ord <- order(aeqSurv(Surv(df_rott$tte, df_rott$event))[, "time"])
c(IF = max(abs(IF - pw$IF[ord, ])), dw = max(abs(dw - pw$dw[ord, ])))
```

    ## IF dw 
    ##  0  0

### Survival difference at one time point

``` r

t_star <- 60    # months

# psi_i(t_star) of section 4.1 for the subjects of arm a (0 for the others), and S^a(t_star)
km_contrib <- function(a) {
  i <- which(A == a)
  s <- sort(unique(U[i][D[i] == 1 & U[i] <= t_star]))            # event times of the arm up to t_star
  Y <- sapply(s, function(x) sum(w[i] * (U[i] >= x)))            # weighted risk set
  dN <- sapply(s, function(x) sum(w[i] * D[i] * (U[i] == x)))    # weighted events
  psi <- numeric(length(U))
  psi[i] <- sapply(i, function(k) {
    own <- if (D[k] == 1 && U[k] <= t_star) 1 / (Y[s == U[k]] - dN[s == U[k]]) else 0
    own - sum((dN / (Y * (Y - dN)))[s <= U[k]])
  })
  list(S = prod(1 - dN / Y), psi = psi)
}
k1 <- km_contrib(1); k0 <- km_contrib(0)

# corrected contributions Phi_i^a = 1(A_i = a) w_i psi_i + IF_i' D_a, with D_a = sum_k psi_k dw_k
Phi <- function(k) w * k$psi + as.numeric(IF %*% colSums(k$psi * dw))
xi <- k1$S * Phi(k1) - k0$S * Phi(k0)
se_diff_hand <- sqrt(sum(xi^2))

fit_est <- KM_diff(df_rott, "tte", "event", "treat", weight.name = "sw.weights", at_points = t_star,
                   weight_variance = "estimated", weight_model = pw)
fit_known <- KM_diff(df_rott, "tte", "event", "treat", weight.name = "sw.weights", at_points = t_star)
se_diff_pkg <- sqrt(fit_est$sig2_dhat)
```

The contributions are centred (section 4.1):
$`\sum_{i: A_i = a} w_i \psi_i(t)`$ is -7.3e-18 in the treated arm and
-1.9e-17 in the control arm.

### Log hazard ratio, FH(0, 0)

``` r

cox_est <- cox_rhogamma(dc_est, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))
cox_known <- cox_rhogamma(dc_known, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))
bhat <- cox_est$fit$bhat

# psi_i(beta) and the information I_w(beta) of section 4.4, time weight kappa = 1
score_contrib <- function(beta) {
  s <- sort(unique(U[D == 1])); eb <- exp(beta * A)
  S0 <- sapply(s, function(x) sum(w * eb * (U >= x)))
  Abar <- sapply(s, function(x) sum(w * eb * A * (U >= x))) / S0
  dL <- sapply(s, function(x) sum(w * D * (U == x))) / S0
  psi <- sapply(seq_along(U), function(k) {
    own <- if (D[k] == 1) A[k] - Abar[s == U[k]] else 0
    w[k] * (own - eb[k] * sum(((A[k] - Abar) * dL)[s <= U[k]]))
  })
  list(psi = psi, info = sum(Abar * (1 - Abar) * S0 * dL))
}
sc_bhat <- score_contrib(bhat); sc_zero <- score_contrib(0)

# corrected contributions psi_i^ew = psi_i + IF_i' B, with B = sum_j (psi_j / w_j) dw_j
psi_ew <- function(psi) psi + as.numeric(IF %*% colSums((psi / w) * dw))
se_loghr_hand <- sqrt(sum(psi_ew(sc_bhat$psi)^2)) / sc_bhat$info
z_hand <- -sum(sc_zero$psi) / sqrt(sum(psi_ew(sc_zero$psi)^2))
```

The known-weight contributions computed here from the formula equal
those the package uses:
[`cox_rhogamma_resample()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma_resample.md)
returns $`-\sum_i G_i\, \psi_i`$ for supplied multipliers, so unit
multipliers return the contributions one by one.

``` r

n <- length(U)
unit <- cox_rhogamma_resample(cox_known$fit, i_bhat = sc_bhat$info, K_wt_rg = NULL, i_zero = NULL, K_zero = NULL,
                              G0.draws = diag(n)[A == 0, ], G1.draws = diag(n)[A == 1, ], draws = n)
c(psi_at_bhat = max(abs(-unit$score_star - sc_bhat$psi)), psi_at_zero = max(abs(-unit$score_star_null - sc_zero$psi)))
```

    ##  psi_at_bhat  psi_at_zero 
    ## 2.220446e-16 4.440892e-16

### By hand against the package

``` r

cmp <- data.frame(
  quantity = c(paste0("SE of the survival difference at ", t_star, " months"),
               "SE of the log hazard ratio, FH(0,0)",
               "Log-rank z, FH(0,0)"),
  by_hand = c(se_diff_hand, se_loghr_hand, z_hand),
  package = c(se_diff_pkg, cox_est$fit$sig_bhat_asy, cox_est$z.score),
  known = c(sqrt(fit_known$sig2_dhat), cox_known$fit$sig_bhat_asy, cox_known$z.score)
)
cmp$difference <- cmp$by_hand - cmp$package
knitr::kable(
  data.frame(cmp$quantity, fmt(cmp$by_hand, 10), fmt(cmp$package, 10),
             formatC(cmp$difference, format = "e", digits = 1), fmt(cmp$known, 10)),
  col.names = c("Quantity", "By hand, estimated", "Package, estimated", "Difference", "Package, known"),
  align = "lrrrr")
```

| Quantity | By hand, estimated | Package, estimated | Difference | Package, known |
|:---|---:|---:|---:|---:|
| SE of the survival difference at 60 months | 0.0351664399 | 0.0351664399 | 1.4e-17 | 0.0373479646 |
| SE of the log hazard ratio, FH(0,0) | 0.0840255194 | 0.0840255194 | 0.0e+00 | 0.0912795746 |
| Log-rank z, FH(0,0) | 4.8686437524 | 4.8686437524 | 8.9e-16 | 4.4490619446 |

``` r

stopifnot(all(abs(cmp$difference) < 1e-8))
```

The estimates are the same under the two options: the survival
difference at 60 months is 0.1503 (by hand 0.1503) and the hazard ratio
is 0.6648. In this example the standard error of the survival difference
under `"estimated"` is 5.8% smaller than under `"known"`, and that of
the log hazard ratio 7.9% smaller.

### The three terms of the variance

The decomposition of section 5, and the sample counterpart of the
projected variance (*background*: not computed by the package).

``` r

M <- crossprod(IF)
terms <- function(c0, b) {
  C <- colSums(c0 * IF)
  c(known = sum(c0^2), cross = 2 * sum(b * C), quadratic = as.numeric(b %*% M %*% b),
    projected = sum(c0^2) - as.numeric(C %*% solve(M, C)))
}
dec <- rbind(
  terms(w * (k1$S * k1$psi - k0$S * k0$psi), k1$S * colSums(k1$psi * dw) - k0$S * colSums(k0$psi * dw)),
  terms(sc_bhat$psi, colSums((sc_bhat$psi / w) * dw)))
dec <- cbind(dec, estimated = dec[, "known"] + dec[, "cross"] + dec[, "quadratic"])
dec_tab <- data.frame(
  c(paste0("Survival difference at ", t_star, " months"), "Score at the estimate, FH(0,0)"),
  formatC(dec[, "known"], format = "g", digits = 6), formatC(dec[, "cross"], format = "g", digits = 6),
  formatC(dec[, "quadratic"], format = "g", digits = 6), formatC(dec[, "estimated"], format = "g", digits = 6),
  formatC(dec[, "projected"], format = "g", digits = 6))
knitr::kable(dec_tab, align = "lrrrrr",
             col.names = c("Variance of", "Known: $\\sum c_i^2$", "$2\\, b^\\top C$", "$b^\\top M\\, b$",
                           "Estimated (sum)", "Projected: $\\sum c_i^2 - C^\\top M^{-1} C$"))
```

| Variance of | Known: $`\sum c_i^2`$ | $`2\, b^\top C`$ | $`b^\top M\, b`$ | Estimated (sum) | Projected: $`\sum c_i^2 - C^\top M^{-1} C`$ |
|:---|---:|---:|---:|---:|---:|
| Survival difference at 60 months | 0.00139487 | -0.000325549 | 0.000167357 | 0.00123668 | 0.00123547 |
| Score at the estimate, FH(0,0) | 118.593 | -36.8974 | 18.797 | 100.493 | 100.374 |

In both rows the cross term is negative and larger in absolute value
than the quadratic term, so the `"estimated"` variance is below the
`"known"` one, and it is close to the projected value (ratio 1.001 for
the survival difference and 1.001 for the score).

## Limits and summary

**Covered by `weight_variance = "estimated"`.**

- Model class: a logistic `glm` for a binary treatment (family binomial,
  logit link, no prior weights, no offset), fitted on exactly the rows
  analysed.
- Weight types: `"ipw"` and `"stabilized"`.
- Estimators: the weighted Kaplan–Meier curve of an arm, of an arm
  within a stratum and of all subjects; the survival difference with its
  pointwise limits and simultaneous band; the cumulative RMST difference
  with its limits and band; the weighted log-rank test, stratified or
  not, cumulative or not; the weighted Cox estimate of
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  with its asymptotic and resampling intervals; the `survival`-based
  outputs of
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md).
- Time-weight schemes: all of them run.

**Not covered.**

- Weights from any other model, and weights truncated or modified after
  [`ps_weights()`](https://larry-leon.github.io/weightedsurv/reference/ps_weights.md):
  use `"known"`.
- A subgroup analysed with weights from a model fitted on the full
  sample: the model’s rows are not the rows analysed; use `"known"`.
- The estimation of the arm proportion in stabilized weights (no effect
  on the Kaplan–Meier quantities; ignored for the log-rank and Cox
  score).
- The dependence of the time weights on the estimated subject weights
  (every scheme other than FH(0, 0) and `"custom_time"`).
- The resampling-based pointwise variance, `draws > 0` in
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
  [`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md)
  and
  [`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md):
  an error under `"estimated"`. `draws.band`, `draws_sb` and `draws` in
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)
  are available.
- Stratified results:
  [`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
  [`wlr_cumulative()`](https://larry-leon.github.io/weightedsurv/reference/wlr_cumulative.md)
  and
  [`wlr_dhat_estimates()`](https://larry-leon.github.io/weightedsurv/reference/wlr_dhat_estimates.md)
  compute unstratified statistics under both options; the stratified
  correction is in
  [`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
  only.
- Concentrated weights: the correction does not repair them.

**Summary.**

1.  `"estimated"` adds to each subject’s contribution the term
    $`IF_i^\top \times`$ (derivative of the estimator in the
    coefficients of the propensity model). The derivative is exact
    (Lemmas 1 and 2); the resulting variance rests on a first-order
    expansion.
2.  The variance is the sum of squared corrected contributions over all
    $`n`$ subjects. The two arms, and the strata, are no longer
    independent.
3.  With a correct model the large-sample variance is no larger than
    under `"known"`, by a projection argument.
4.  Resampling multiplies the corrected contributions; nothing is
    refitted.
5.  In the evaluation, `"estimated"` reproduces a refitting bootstrap,
    is close to nominal with evenly spread weights, and falls short, as
    `"known"` does, with concentrated weights.

Related articles: `weightedkm_methodology` (weighted Kaplan–Meier curves
and bands), `weightedcox_methodology` (weighted Cox estimation),
`weightedrmst_methodology` (RMST), and the `weightedsurv_examples`
vignette (the worked analyses). The IPTW-weighted Kaplan–Meier estimator
and log-rank test go back to Xie and Liu (2005).

## References

Deng, Yuhao, and Rui Wang. 2025. “Adjusted Nelson–Aalen Estimators by
Inverse Treatment Probability Weighting with an Estimated Propensity
Score.” *Statistics in Medicine* 44 (10–12): e70085.
<https://doi.org/10.1002/sim.70085>.

Dobler, D., J. Beyersmann, and M. Pauly. 2017. “Non-Strange Weird
Resampling for Complex Survival Data.” *Biometrika* 104 (3): 699–711.

Fleming, T. R., and D. P. Harrington. 1991. *Counting Processes and
Survival Analysis*. Wiley Series in Probability and Statistics. Wiley.

Hirano, Keisuke, Guido W. Imbens, and Geert Ridder. 2003. “Efficient
Estimation of Average Treatment Effects Using the Estimated Propensity
Score.” *Econometrica* 71 (4): 1161–89.
<https://doi.org/10.1111/1468-0262.00442>.

Lin, D. Y., and L. J. Wei. 1989. “The Robust Inference for the Cox
Proportional Hazards Model.” *Journal of the American Statistical
Association* 84: 1074–78.

Lin, Dan-Yu, L. J. Wei, and Zhiliang Ying. 1993. “Checking the Cox Model
with Cumulative Sums of Martingale-Based Residuals.” *Biometrika* 80
(3): 557–72.

Lunceford, Jared K., and Marie Davidian. 2004. “Stratification and
Weighting via the Propensity Score in Estimation of Causal Treatment
Effects: A Comparative Study.” *Statistics in Medicine* 23 (19):
2937–60. <https://doi.org/10.1002/sim.1903>.

Xie, Jun, and Chaofeng Liu. 2005. “Adjusted Kaplan–Meier Estimator and
Log-Rank Test with Inverse Probability of Treatment Weighting for
Survival Data.” *Statistics in Medicine* 24 (20): 3089–110.
<https://doi.org/10.1002/sim.2174>.
