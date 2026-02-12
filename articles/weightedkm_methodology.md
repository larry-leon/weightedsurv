# Weighted Kaplan-Meier Estimation Methodology

## Introduction

In clinical trials and observational studies, the Kaplan–Meier (KM)
estimator is the most widely used nonparametric method for estimating
survival functions from censored time-to-event data. When treatment
groups are formed by randomization, the standard KM estimator provides
consistent estimation of treatment-specific survival distributions.
However, in observational studies—where treatment assignment may depend
on baseline covariates—direct KM comparisons can be misleading due to
confounding.

Inverse probability of treatment weighting (IPTW) provides a principled
approach to address confounding by creating a pseudo-population in which
baseline covariates are balanced across treatment groups. The weighted
Kaplan–Meier estimator, constructed by applying propensity-score-derived
weights to the at-risk and event counting processes, yields consistent
estimation of population-level (counterfactual) survival distributions
under standard causal assumptions (Cole and Hernán 2004; Austin 2014).

This article reviews the methodology underlying weighted KM estimation
as implemented in the `weightedsurv` package, connecting the classical
counting-process framework to the causal inference formulation. We draw
extensively on the theoretical foundations established by Deng and Wang
(2025) for the adjusted Nelson–Aalen estimator, as well as the
martingale resampling approach of Dobler, Beyersmann, and Pauly (2017)
used in `weightedsurv` for simultaneous confidence bands.

## Notation and Setup

### Observed Data

We observe $`n`$ independent subjects from two treatment groups
($`A = 0`$ for control, $`A = 1`$ for active treatment). For each
subject $`i`$, the observed data are $`(A_i, X_i, T_i, \Delta_i)`$,
where:

- $`A_i \in \{0, 1\}`$: treatment indicator,
- $`X_i`$: vector of baseline covariates,
- $`T_i = \min(\tilde{T}_i, C_i)`$: observed follow-up time (the minimum
  of the event time $`\tilde{T}_i`$ and censoring time $`C_i`$),
- $`\Delta_i = I(\tilde{T}_i \leq C_i)`$: event indicator.

### Counting Process Formulation

In the counting process framework (Fleming and Harrington 1991), the
individual-level processes are:

- At-risk process: $`Y_i(t) = I(T_i \geq t)`$,
- Event counting process: $`N_i(t) = I(T_i \leq t, \Delta_i = 1)`$.

For two groups $`j = 0, 1`$, the aggregate processes are:

``` math
\bar{Y}_j(t) = \sum_{i: A_i = j} Y_i(t), \qquad \bar{N}_j(t) = \sum_{i: A_i = j} N_i(t).
```

The Nelson–Aalen estimator of the cumulative hazard for group $`j`$ is

``` math
\hat{\Lambda}_j(t) = \int_0^t \frac{d\bar{N}_j(s)}{\bar{Y}_j(s)},
```

and the Kaplan–Meier estimator is obtained via the product-limit
transformation:

``` math
\hat{S}_j(t) = \prod_{s \leq t} \left\{1 - \frac{d\bar{N}_j(s)}{\bar{Y}_j(s)}\right\}.
```

Under random censoring within each treatment group, both estimators are
consistent for their respective group-specific survival functions.

## Standard (Unweighted) Kaplan–Meier Estimator

### Definition

The Kaplan–Meier estimator (Kaplan and Meier 1958) is defined at the
ordered event times $`t_1 < t_2 < \cdots`$ as:

``` math
\hat{S}(t) = \prod_{t_k \leq t} \left(1 - \frac{d_k}{n_k}\right),
```

where $`d_k`$ denotes the number of events at time $`t_k`$ and $`n_k`$
is the number at risk just before $`t_k`$.

### Martingale Representation

Under the null hypothesis of a common hazard $`\lambda(t)`$ in group
$`j`$, the counting process martingale is defined as:

``` math
M_{i,j}(t) = N_{i,j}(t) - \int_0^t Y_{i,j}(s) \lambda_j(s)\, ds,
```

where $`Y_{i,j}(s) = I(T_i \geq s, A_i = j)`$ and
$`N_{i,j}(t) = I(T_i \leq t, \Delta_i = 1, A_i = j)`$.

The process $`M_{i,j}(t)`$ is a zero-mean martingale with respect to the
natural filtration. This representation is central to deriving the
variance of the KM estimator and forms the basis for the resampling
methods described below.

### Greenwood’s Variance Estimator

The classical variance estimator for $`\hat{S}_j(t)`$ is Greenwood’s
formula:

``` math
\widehat{\text{Var}}\bigl[\hat{S}_j(t)\bigr] = \hat{S}_j(t)^2 \sum_{t_k \leq t} \frac{d_k}{n_k(n_k - d_k)}.
```

## Causal Framework for Weighted Estimation

### Potential Outcomes and Causal Estimands

Following the potential outcomes framework, let $`\tilde{T}^a`$ denote
the potential survival time under treatment $`a \in \{0,1\}`$. The
population-level counterfactual survival function is:

``` math
S^a(t) = P(\tilde{T}^a > t),
```

which represents the survival probability that would be observed if the
entire population were assigned treatment $`a`$. This is the target
estimand for weighted KM estimation.

### Identifying Assumptions

Identification of $`S^a(t)`$ from observed data requires the following
assumptions, formalized for example by Deng and Wang (2025):

**Assumption 1 (Ignorability).**
$`A \perp\!\!\!\perp (\tilde{T}^a, \tilde{C}^a) \mid X`$, for
$`a = 0, 1`$. Treatment assignment is independent of the potential
outcomes given baseline covariates.

**Assumption 2 (Completely random censoring).**
$`\tilde{C}^a \perp\!\!\!\perp \tilde{T}^a`$, for $`a = 0, 1`$. The
censoring time is unconditionally independent of the event time. This
assumption is essential for nonparametric estimation. As discussed by
Deng and Wang (2025), it may be relaxed to conditional random censoring,
but at the cost of requiring additional modeling of the censoring
mechanism.

**Assumption 3 (Positivity).** $`0 < c < e(X) < 1 - c < 1`$ for a
constant $`c > 0`$, where $`e(X) = P(A = 1 \mid X)`$ is the propensity
score. Every subject has a positive probability of receiving either
treatment.

**Assumption 4 (Consistency).** $`T = T^A,\; \Delta = \Delta^A`$. The
observed outcome equals the potential outcome corresponding to the
treatment actually received.

## Weighted Kaplan–Meier Estimation via IPTW

### Propensity Score and IPW Weights

The propensity score $`e(X) = P(A = 1 \mid X)`$ summarizes the covariate
information relevant to treatment assignment (Rosenbaum and Rubin 1983).
The inverse probability weights are defined as:

``` math
w(a; A, X) = \frac{I(A = a)}{e(a; X)},
```

where $`e(1; X) = e(X)`$ and $`e(0; X) = 1 - e(X)`$.

**Stabilized weights.** To reduce variability, one commonly uses
stabilized weights:

``` math
w^{\text{stab}}(a; A, X) = \frac{P(A = a)}{e(a; X)} \cdot I(A = a),
```

where $`P(A = a)`$ is the marginal treatment probability. In practice,
both $`e(X)`$ and $`P(A = a)`$ are estimated from the data. When
stabilized weights are used, the weighted risk set at baseline equals
the original sample size, which facilitates interpretation (Cole and
Hernán 2004).

### Weighted Counting Processes

By applying IPW weights to the individual counting processes, we
construct weighted aggregate processes:

``` math
\bar{N}_j^w(t) = \sum_{i: A_i = j} w_i \cdot N_i(t), \qquad \bar{Y}_j^w(t) = \sum_{i: A_i = j} w_i \cdot Y_i(t),
```

where $`w_i`$ is the propensity score weight for subject $`i`$.

### Adjusted Nelson–Aalen Estimator

The adjusted Nelson–Aalen estimator for the counterfactual cumulative
hazard under treatment $`a`$ is defined as (Deng and Wang 2025; Winnett
and Sasieni 2002):

``` math
\hat{\Lambda}^a(t) = \int_0^t \frac{\sum_{i=1}^n w(a; A_i, X_i) \, dN_i(s)}{\sum_{i=1}^n w(a; A_i, X_i) \, Y_i(s)}.
```

This is equivalent to the standard Nelson–Aalen formula applied to the
weighted counting and at-risk processes.

### Weighted KM Estimator

The weighted Kaplan–Meier estimator is obtained via the product-limit:

``` math
\hat{S}^a(t) = \prod_{s \leq t} \left\{1 - \frac{d\bar{N}_a^w(s)}{\bar{Y}_a^w(s)}\right\}.
```

When there is a single terminal event (no competing risks), the weighted
KM estimator and the survival function derived from the adjusted
Nelson–Aalen estimator
$`\hat{S}^a_{\text{NA}}(t) = \exp\{-\hat{\Lambda}^a(t)\}`$ are
asymptotically equivalent (Deng and Wang 2025).

### Practical Weight Computation

In the `weightedsurv` package, propensity score weights are computed as
stabilized IPTW weights. A practical recipe, following Cole and Hernán
(2004) and as implemented in the `weightedsurv_examples` vignette, is:

1.  Fit a propensity score model:
    $`\hat{e}(X_i) = \text{expit}(\hat{\beta}^T X_i)`$ via logistic
    regression.
2.  Compute the marginal probability:
    $`\hat{P}(A = 1) = \hat{\pi}_{\text{null}}`$ from an intercept-only
    model.
3.  Calculate stabilized weights:
    ``` math
    w_i = \begin{cases} \hat{\pi}_{\text{null}} / \hat{e}(X_i) & \text{if } A_i = 1, \\ (1 - \hat{\pi}_{\text{null}}) / (1 - \hat{e}(X_i)) & \text{if } A_i = 0. \end{cases}
    ```
4.  (Optional) Truncate extreme weights to reduce variance, for example
    at the 5th and 95th percentiles.

## Asymptotic Theory

### Known Propensity Score

When the propensity score $`e(X)`$ is known (as in a randomized trial
with known allocation probabilities), Deng and Wang (2025) show that the
weighted counting process admits a martingale decomposition.
Specifically, for each subject $`i`$ in treatment group $`a`$:

``` math
w(a; A_i, X_i) \, dN_i(t) = w(a; A_i, X_i) \, Y_i(t) \, d\Lambda^a(t) + w(a; A_i, X_i) \, dM_i^a(t),
```

where $`M_i^a(t)`$ is a martingale. This yields the finite-sample
unbiasedness of $`d\hat{\Lambda}^a(t)`$ for $`d\Lambda^a(t)`$.

The variance of $`\hat{\Lambda}^a(t)`$ follows from martingale theory:

``` math
\text{Var}\bigl[\hat{\Lambda}^a(t)\bigr] \approx \int_0^t \frac{E\bigl[w(a; A, X)^2 \, Y(s) \, d\Lambda^a(s)\bigr]}{\bigl(E\bigl[w(a; A, X) \, Y(s)\bigr]\bigr)^2}.
```

### Estimated Propensity Score

In practice, the propensity score is unknown and must be estimated. A
key contribution of Deng and Wang (2025) is deriving the influence
function for $`\hat{\Lambda}^a(t)`$ when using an estimated propensity
score $`\hat{e}(X)`$. If the estimated propensity score is regular and
asymptotically linear (RAL), the additional variance term due to
propensity score estimation can be characterized.

The influence function of $`\hat{\Lambda}^a(t)`$ decomposes as:

``` math
\sqrt{n}\bigl(\hat{\Lambda}^a(t) - \Lambda^a(t)\bigr) = \frac{1}{\sqrt{n}} \sum_{i=1}^n \psi_i^a(t) + o_p(1),
```

where $`\psi_i^a(t)`$ is the influence function that incorporates both
the martingale term (from the counting process) and the correction term
(from estimating $`e(X)`$).

Importantly, Deng and Wang (2025) find through simulation that the
additional variance contribution from the estimated propensity score is
typically small. This is because the correction term takes the form of a
weighted martingale whose expectation is nearly zero.

### Variance Estimation in `weightedsurv`

In the `weightedsurv` package, standard errors for the weighted KM
estimator are computed via a robust variance formula analogous to
Greenwood’s formula, applied to the weighted counting processes. These
are validated against the `adjustedCurves` package and `survfit` with
case-weights (see the `weightedsurv_examples` vignette for
cross-checks).

For settings involving complex weighting and the simultaneous inference
procedures described below, variance estimation via martingale
resampling is preferred.

## Inference for Survival Differences

### Pointwise Confidence Intervals

For comparing two treatment groups, the survival difference
$`\Delta(t) = S^1(t) - S^0(t)`$ is estimated by
$`\hat{\Delta}(t) = \hat{S}^1(t) - \hat{S}^0(t)`$. Pointwise confidence
intervals are constructed as:

``` math
\hat{\Delta}(t) \pm z_{\alpha/2} \cdot \widehat{\text{se}}\bigl[\hat{\Delta}(t)\bigr],
```

where $`\widehat{\text{se}}[\hat{\Delta}(t)]`$ is obtained from the
(independent) variance estimates of $`\hat{S}^1(t)`$ and
$`\hat{S}^0(t)`$.

### Simultaneous Confidence Bands via Martingale Resampling

Pointwise confidence intervals provide valid inference at any single
time point, but the survival difference $`\hat{\Delta}(t)`$ is evaluated
across a range of time points $`t`$. Simultaneous confidence bands
guarantee coverage uniformly over an interval of time.

The `weightedsurv` package implements simultaneous confidence bands
using the martingale resampling approach of Dobler, Beyersmann, and
Pauly (2017).

#### Martingale Resampling Principle

The key idea is to approximate the distribution of the centered process:

``` math
\sqrt{n}\bigl(\hat{\Delta}(t) - \Delta(t)\bigr), \quad t \in [\tau_L, \tau_U],
```

by generating resampled versions of the process. Under the null (or
conditionally on the observed data), the martingale increments
$`dM_{i,j}(t)`$ are replaced by their observable counterparts multiplied
by independent standard normal variates.

The resampled process is:

``` math
W^\dagger = \sum_{i=1}^{n_0} \int_0^\infty \frac{K(t)}{\bar{Y}_0(t)} \, dN_{i,0}(t) \cdot G_{i,0} - \sum_{i=1}^{n_1} \int_0^\infty \frac{K(t)}{\bar{Y}_1(t)} \, dN_{i,1}(t) \cdot G_{i,1},
```

where $`G_{i,j}`$ are i.i.d. $`N(0,1)`$ random variables, independent of
the data.

The observable counting process increments $`dN_{i,j}(t)`$ serve as
proxies for the unobservable martingale increments, and the Gaussian
multipliers capture the stochastic variability. Unknown nuisance
parameters are replaced by consistent estimates (“plug-in”).

#### Algorithm

The procedure, as implemented in `plotKM.band_subgroups` and `KM_diff`,
proceeds as follows:

1.  Compute the observed survival difference $`\hat{\Delta}(t)`$ at each
    event time within the quantile-trimmed range $`[t_{q_L}, t_{q_U}]`$
    (controlled by the `qtau` parameter).
2.  For $`b = 1, \ldots, B`$ resamples:
    1.  Generate $`n`$ independent $`N(0,1)`$ variates
        $`\{G_i^{(b)}\}`$.
    2.  Compute the resampled survival difference process
        $`\hat{\Delta}^{(b)}(t)`$ by perturbing the counting process
        increments.
3.  For each resample $`b`$, compute
    $`\sup_t |\hat{\Delta}^{(b)}(t) - \hat{\Delta}(t)|`$.
4.  The $`(1-\alpha)`$ simultaneous band is:
    ``` math
    \hat{\Delta}(t) \pm c_{1-\alpha},
    ```
    where $`c_{1-\alpha}`$ is the $`(1-\alpha)`$ quantile of
    $`\sup_t |\hat{\Delta}^{(b)}(t) - \hat{\Delta}(t)|`$ over the $`B`$
    resamples.

#### Extension to Weighted Setting

The resampling framework naturally accommodates IPTW weights. In the
weighted case, the counting process increments are replaced by their
weighted counterparts, and the variance proxy accounts for the
propensity-score weights. This is operationally handled in
`weightedsurv` by passing a `weight.name` argument to the relevant
functions (e.g., `plotKM.band_subgroups`, `KM_diff`,
`cumulative_rmst_bands`).

## Cumulative RMST Inference

### RMST Definition

The restricted mean survival time (RMST) up to a truncation time
$`\tau`$ is defined as:

``` math
\mu^a(\tau) = \int_0^{\tau} S^a(t)\, dt = E\bigl[\min(\tilde{T}^a, \tau)\bigr].
```

The RMST difference $`\mu^1(\tau) - \mu^0(\tau)`$ is an interpretable,
model-free summary of the treatment effect (Uno et al. 2014).

### Cumulative RMST Curves

Rather than evaluating the RMST at a single $`\tau`$, the
`cumulative_rmst_bands` function in `weightedsurv` computes
$`\hat{\mu}^1(\tau) - \hat{\mu}^0(\tau)`$ as a function of $`\tau`$
across the entire follow-up range. This provides a “cumulative RMST
curve” that reveals how the treatment benefit (or harm) accumulates over
time.

Simultaneous confidence bands for the RMST difference curve are obtained
by the same martingale resampling approach, applied to the integrated
survival difference process.

## Connection to Related Approaches

### Relationship to the Adjusted Nelson–Aalen Estimator

As formalized by Deng and Wang (2025), the adjusted Nelson–Aalen
estimator and the weighted KM estimator are asymptotically equivalent
for a single terminal event. The Nelson–Aalen approach avoids product
limits in the estimator and connects directly to the hazard-based
identification results. The `weightedsurv` package works primarily with
the product-limit (KM) form, which is more commonly used in clinical
practice and has the natural interpretation as an estimated probability.

### Weighted Risk Set Estimators

An alternative approach to estimating treatment-specific survival
distributions from sequential randomization designs was proposed by Guo
and Tsiatis (2005). The weighted risk set (WRS) estimator uses inverse
probability weights within the risk sets of a Nelson–Aalen-type
estimator. Miyahara and Wahed (2010) proposed weighted KM estimators
(with both fixed and time-dependent weights) for two-stage treatment
regimes and showed that both forms are asymptotically unbiased.

These WRS and WKM estimators share the same fundamental structure as the
IPTW KM estimator: weight each subject’s contribution to the risk set
and event count by an inverse probability weight reflecting the
probability of following the treatment regime of interest.

### Doubly-Robust Estimators

Bai, Tsiatis, and O’Brien (2013) developed doubly-robust estimators for
treatment-specific survival distributions that are consistent if
*either* the propensity score model or a model for the survival
distribution as a function of covariates is correctly specified. These
augmented IPTW estimators can achieve greater efficiency than pure IPTW
methods. While the current `weightedsurv` implementation focuses on
IPTW, the counting-process infrastructure is amenable to such
extensions.

### Covariate-Adjusted Log-Rank Tests

In the randomized trial setting, Ye, Shao, and Yi (2024) developed
covariate-adjusted log-rank tests with guaranteed efficiency gains over
the unadjusted test. A key result is the universal applicability of
their method to different randomization schemes (simple randomization,
stratified permuted block, Pocock–Simon minimization). The weighted
log-rank statistics in `weightedsurv` are complementary: while the
package’s Fleming–Harrington class statistics target non-proportional
hazards alternatives, the covariate-adjusted framework of Ye, Shao, and
Yi (2024) targets efficiency gains through covariate adjustment under
any hazard structure.

### Adjusted Nelson–Aalen with Retrospective Matching

Winnett and Sasieni (2002) developed a weighted Nelson–Aalen estimator
for retrospectively matched data, with stratum-specific random effects.
Their work highlighted that ignoring heterogeneity from matched strata
leads to variance underestimation—a cautionary result relevant to any
weighted survival analysis where the source of weights introduces
additional variability.

### Covariate-Adjusted Logrank Tests and Working Models

Kong and Slud (1997) showed how to construct robust covariate-adjusted
logrank statistics by substituting maximum partial likelihood estimators
from various working models, providing valid hypothesis tests through
robust variance estimators even under model misspecification. This work
provides theoretical foundation for the approach taken in `weightedsurv`
where working models need not be exactly correctly specified.

### Sample Size Formulae for Weighted Designs

Li and Murphy (2011) developed sample size formulae for two-stage
randomized trials with survival outcomes based on both a weighted KM
estimator of survival probabilities at a fixed time-point and a weighted
version of the log-rank test. Their conservative formulae demonstrate
that the weighted KM framework extends naturally to the design phase of
clinical trials.

## Implementation in `weightedsurv`

The weighted KM methodology described above is implemented through
several interconnected functions in the `weightedsurv` package. The key
functions and their roles are:

- **`df_counting`** / **`get_dfcounting`**: Prepare the counting-process
  dataset from raw survival data. When a `weight.name` argument is
  supplied, these functions construct the weighted risk sets and event
  counts needed for weighted KM estimation.

- **`plot_weighted_km`**: Plots weighted (or unweighted) KM survival
  curves for two groups, with optional confidence intervals, risk
  tables, and summary statistics (log-rank p-value, Cox HR).

- **`KM_diff`**: Computes the KM survival difference
  $`\hat{S}^1(t) - \hat{S}^0(t)`$, with both pointwise and simultaneous
  confidence intervals. Accepts arbitrary non-negative weights for IPTW.

- **`plotKM.band_subgroups`**: Produces publication-ready plots of the
  survival difference with simultaneous confidence bands via martingale
  resampling.

- **`cumulative_rmst_bands`**: Computes cumulative RMST difference
  curves with simultaneous confidence bands, again using the martingale
  resampling framework.

For comprehensive examples demonstrating both the randomized trial
(GBSG) and observational (Rotterdam) settings, including
cross-validation of weighted SEs against the `adjustedCurves` and
`survfit` packages, see the `weightedsurv_examples` vignette.

## Summary

The weighted KM estimator provides a principled nonparametric approach
to estimating population-level survival distributions from observational
data via IPTW. The methodology rests on:

1.  **Causal identification** via the ignorability, completely random
    censoring, positivity, and consistency assumptions.
2.  **Counting process formulation** enabling martingale-based variance
    estimation and resampling inference.
3.  **Martingale resampling** for simultaneous confidence bands,
    extending naturally to the weighted setting and to cumulative RMST
    curves.

The `weightedsurv` package provides an integrated implementation of
these methods, with careful attention to numerical validation against
established packages.

## References

Austin, Peter C. 2014. “Use of the Propensity Score for Estimating
Comparative Effectiveness in Observational Studies.” *Statistics in
Medicine* 33 (13): 2202–11. <https://doi.org/10.1002/sim.6051>.

Bai, Xiaofei, Anastasios A. Tsiatis, and Sean M. O’Brien. 2013.
“Doubly-Robust Estimators of Treatment-Specific Survival Distributions
in Observational Studies with Stratified Sampling.” *Biometrics* 69 (4):
830–39. <https://doi.org/10.1111/biom.12076>.

Cole, Stephen R., and Miguel A. Hernán. 2004. “Adjusted Survival Curves
with Inverse Probability Weights.” *Computer Methods and Programs in
Biomedicine* 75 (1): 45–49.
<https://doi.org/10.1016/j.cmpb.2003.10.004>.

Deng, Yuhao, and Rui Wang. 2025. “Adjusted Nelson–Aalen Estimators by
Inverse Treatment Probability Weighting with an Estimated Propensity
Score.” *Statistics in Medicine* 44 (10–12): e70085.
<https://doi.org/10.1002/sim.70085>.

Dobler, D., J. Beyersmann, and M. Pauly. 2017. “Non-Strange Weird
Resampling for Complex Survival Data.” *Biometrika* 104 (3): 699–711.

Fleming, T. R., and D. P. Harrington. 1991. *Counting Processes and
Survival Analysis*. Wiley Series in Probability and Statistics. Wiley.

Guo, Xiang, and Anastasios Tsiatis. 2005. “A Weighted Risk Set Estimator
for Survival Distributions in Two-Stage Randomization Designs with
Censored Survival Data.” *The International Journal of Biostatistics* 1
(1): Article 1. <https://doi.org/10.2202/1557-4679.1000>.

Kaplan, E. L., and Paul Meier. 1958. “Nonparametric Estimation from
Incomplete Observations.” *Journal of the American Statistical
Association* 53 (282): 457–81.

Kong, Fan Hui, and Eric Slud. 1997. “Robust Covariate-Adjusted Logrank
Tests.” *Biometrika* 84 (4): 847–62.

Li, Zhiguo, and Susan A. Murphy. 2011. “Sample Size Formulae for
Two-Stage Randomized Trials with Survival Outcomes.” *Biometrika* 98
(3): 503–18. <https://doi.org/10.1093/biomet/asr019>.

Miyahara, Sachiko, and Abdus S. Wahed. 2010. “Weighted Kaplan–Meier
Estimators for Two-Stage Treatment Regimes.” *Statistics in Medicine*
29: 2581–91. <https://doi.org/10.1002/sim.4020>.

Rosenbaum, Paul R., and Donald B. Rubin. 1983. “The Central Role of the
Propensity Score in Observational Studies for Causal Effects.”
*Biometrika* 70 (1): 41–55.

Uno, Hajime, Brian Claggett, Lu Tian, Eiji Inoue, Peter Gallo, Toshio
Miyata, Deborah Schrag, et al. 2014. “Moving Beyond the Hazard Ratio in
Quantifying the Between-Group Difference in Survival Analysis.” *Journal
of Clinical Oncology* 32 (22): 2380–85.
<https://doi.org/10.1200/JCO.2014.55.2208>.

Winnett, Angela, and Peter Sasieni. 2002. “Adjusted Nelson–Aalen
Estimates with Retrospective Matching.” *Journal of the American
Statistical Association* 97 (457): 245–56.

Ye, Ting, Jun Shao, and Yanyao Yi. 2024. “Covariate-Adjusted Log-Rank
Test: Guaranteed Efficiency Gain and Universal Applicability.”
*Biometrika* 111 (2): 691–705. <https://doi.org/10.1093/biomet/asad045>.
