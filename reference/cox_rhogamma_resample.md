# Resampling for Weighted Cox Model (rho, gamma)

Performs resampling to estimate uncertainty for the weighted Cox model
(rho, gamma).

## Usage

``` r
cox_rhogamma_resample(
  fit_rhogamma,
  i_bhat,
  K_wt_rg,
  i_zero,
  K_zero,
  G1.draws = NULL,
  G0.draws = NULL,
  draws = 100,
  seedstart = 8316951
)
```

## Arguments

- fit_rhogamma:

  List with fitted Cox model results.

- i_bhat:

  Information at estimated beta.

- K_wt_rg:

  Weights at estimated beta.

- i_zero:

  Information at beta=0.

- K_zero:

  Weights at beta=0.

- G1.draws:

  Optional: pre-generated random draws for groups.

- G0.draws:

  Optional: pre-generated random draws for groups.

- draws:

  Number of resampling iterations (default: 100).

- seedstart:

  Random seed for reproducibility (default: 8316951).

## Value

List with resampling results (score, beta, standard error, etc.).
