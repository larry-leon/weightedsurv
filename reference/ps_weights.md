# Propensity-score weights with the pieces needed to carry their estimation into the inference

Fits (or takes) a logistic regression for a binary treatment and returns
inverse-probability weights together with the influence of the fitted
coefficients and the derivative of each weight in those coefficients.
Passed as `weight_model` to
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`cumulative_rmst_bands()`](https://larry-leon.github.io/weightedsurv/reference/cumulative_rmst_bands.md)
or
[`plotKM.band_subgroups()`](https://larry-leon.github.io/weightedsurv/reference/plotKM.band_subgroups.md)
with `weight_variance = "estimated"`, the object lets standard errors,
intervals, bands and tests account for the fact that the weights were
estimated (see the section "Estimated weights" of
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)).

## Usage

``` r
ps_weights(
  formula = NULL,
  data = NULL,
  type = c("ipw", "stabilized"),
  model = NULL
)

# S3 method for class 'ps_weights'
print(x, digits = 4, ...)
```

## Arguments

- formula:

  A formula for the treatment model, treatment (coded 0/1) on the
  left-hand side. Used with `data`; fitted by
  `stats::glm(formula, family = stats::binomial(), data = data)` with
  the default control, so that weights computed by hand from the same
  call agree.

- data:

  Data frame holding the variables of `formula`: exactly the rows that
  will be analysed.

- type:

  `"ipw"` (default): `1 / e` for treated subjects and `1 / (1 - e)` for
  controls, with `e` the fitted probability of treatment.
  `"stabilized"`: those weights multiplied by the proportion treated
  (treated subjects) and the proportion not treated (controls).

- model:

  An existing logistic `glm` fit, instead of `formula` and `data`. Give
  either `formula` and `data`, or `model`, not both.

- x:

  A `"ps_weights"` object.

- digits:

  Number of significant digits to print.

- ...:

  Not used.

## Value

An object of class `"ps_weights"`, a list with

- weights:

  The weights, one per row.

- type:

  `"ipw"` or `"stabilized"`.

- treat:

  The treatment (0/1), one per row.

- ps:

  The fitted probabilities of treatment.

- IF:

  Matrix, one row per subject and one column per estimated coefficient
  (aliased columns dropped): the influence of each subject on the
  coefficients.

- dw:

  Matrix of the same shape: the derivative of each weight in the
  coefficients.

- model:

  The `glm` fit.

- call:

  The call.

[`print()`](https://rdrr.io/r/base/print.html) shows the type, the
number of rows and of coefficients, and by arm the size, the range of
the weights and the effective sample size \\(\sum w)^2 / \sum w^2\\.

## Details

With \\x_i\\ the row of the model matrix of subject \\i\\, \\A_i\\ the
treatment and \\e_i\\ the fitted probability, \$\$IF_i = \Big\\\sum_j
x_j x_j' e_j (1 - e_j)\Big\\^{-1} x_i (A_i - e_i)\$\$ is the influence
of subject \\i\\ on the fitted coefficients (the estimate minus its
limit is approximately \\\sum_i IF_i\\), and the derivative of the
weight of subject \\i\\ in the coefficients is \\dw_i = -\\(1 - e_i) /
e_i\\ x_i\\ for a treated subject and \\dw_i = \\e_i / (1 - e_i)\\ x_i\\
for a control. With `type = "stabilized"` the weights and `dw` are those
of `"ipw"` multiplied by the arm proportion, which is treated as a fixed
number; `IF` is the same.

The following are errors, each naming the problem: a `model` that is not
a `glm`; a family other than binomial or a link other than logit; prior
weights other than 1 or an offset; a fit that did not converge; a fitted
probability within `1e-10` of 0 or 1; rows dropped for missing values; a
response that does not take both values 0 and 1.

The weights are never truncated or rescaled. The rows of `weights`,
`treat`, `ps`, `IF` and `dw` are in the order of the rows of the data
the model was fitted on.

## See also

[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md)
(section "Estimated weights"),
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md)

Other survival_analysis:
[`KM_diff()`](https://larry-leon.github.io/weightedsurv/reference/KM_diff.md),
[`cox_rhogamma()`](https://larry-leon.github.io/weightedsurv/reference/cox_rhogamma.md),
[`df_counting()`](https://larry-leon.github.io/weightedsurv/reference/df_counting.md),
[`wt.rg.S()`](https://larry-leon.github.io/weightedsurv/reference/wt.rg.S.md)

## Examples

``` r
library(survival)
gbsg$tte <- gbsg$rfstime / 30.4375
pw <- ps_weights(hormon ~ age + meno + size + grade + nodes + pgr + er, data = gbsg,
                 type = "stabilized")
pw
#> Propensity-score weights (stabilized): n = 686, 8 coefficients in the logistic model
#>        arm   n weight min weight max   ess ess / n
#>  treat = 1 246     0.4493       2.63 203.0  0.8253
#>  treat = 0 440     0.7379       1.85 417.7  0.9493
gbsg$sw <- pw$weights

# weights treated as known (the default) and with their estimation carried through
known <- KM_diff(gbsg, "tte", "status", "hormon", weight.name = "sw")
est <- KM_diff(gbsg, "tte", "status", "hormon", weight.name = "sw",
               weight_variance = "estimated", weight_model = pw)
identical(known$dhat, est$dhat)                       # estimates do not change
#> [1] TRUE
summary(sqrt(est$sig2_dhat / known$sig2_dhat))        # standard errors do
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.     NAs 
#>  0.9231  0.9330  0.9445  0.9468  0.9556  0.9979      14 

dc <- df_counting(gbsg, "tte", "status", "hormon", weight.name = "sw",
                  weight_variance = "estimated", weight_model = pw)
dc$ess                                                # effective sample size by arm
#>     group   n      ess     ratio
#> 1   treat 246 203.0308 0.8253283
#> 2 control 440 417.7090 0.9493386
cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))$hr_ci_asy
#>         beta  sig_bhat       hr     lower     upper
#> 1 -0.3747545 0.1276936 0.687458 0.5352457 0.8829562
```
