# weightedsurv — M1b: weighted-SE calibration with a weighted Greenwood estimator

Stage 2, Step 1 (kickoff: `M1b: GO`). Script `dev/audit/stage2/km_m1b_montecarlo.R`, a copy of `dev/audit/km_m1_montecarlo.R` with a fourth estimator added; results in `dev/audit/stage2/km_m1b_results.csv`, console output in `dev/audit/stage2/km_m1b_stdout.txt`. No `R/` file was involved.

Design as in M1, same seeds: 2,000 replicates per cell, n = 200, one arm, T ~ Exp(0.08), C ~ U(0, 30). Runtime 17 s. The three Stage-1 estimators reproduce the Stage-1 `km_m1_results.csv` exactly (all 54 rows identical).

Weighted Greenwood: `SE(t) = S(t) * sqrt( sum_{s<=t} [ sum_i w_i^2 dN_i(s) ] / [ Y_w(s) * (Y_w(s) - dN_w(s)) ] )`, a term being 0 where `Y_w(s) - dN_w(s) <= 0`.

Each entry is SE/SD (mean SE over the empirical SD of Ŝ), then coverage of Ŝ ± 1.96·SE for exp(−0.08t).

| ties | weights | estimator | t = 5 | t = 10 | t = 15 |
|---|---|---|---|---|---|
| continuous | none | package | 1.01 / 0.948 | 0.98 / 0.942 | 1.00 / 0.943 |
| continuous | none | survfit robust | 1.01 / 0.948 | 0.98 / 0.942 | 1.00 / 0.943 |
| continuous | none | documented | 1.00 / 0.947 | 0.97 / 0.942 | 0.99 / 0.939 |
| continuous | none | weighted Greenwood | 1.01 / 0.948 | 0.98 / 0.942 | 1.00 / 0.943 |
| continuous | U(0.3,3) | package | 0.99 / 0.941 | 0.96 / 0.940 | 0.97 / 0.937 |
| continuous | U(0.3,3) | survfit robust | 1.00 / 0.944 | 0.98 / 0.945 | 0.99 / 0.943 |
| continuous | U(0.3,3) | documented | 0.99 / 0.943 | 0.97 / 0.942 | 0.98 / 0.939 |
| continuous | U(0.3,3) | weighted Greenwood | 1.00 / 0.945 | 0.98 / 0.946 | 1.00 / 0.944 |
| continuous | U(0.3,3)/Σ | package | **0.45 / 0.628** | **1.29 / 0.974** | **3.12 / 1.000** |
| continuous | U(0.3,3)/Σ | survfit robust | 1.00 / 0.944 | 0.98 / 0.945 | 0.99 / 0.943 |
| continuous | U(0.3,3)/Σ | documented | 0.99 / 0.943 | 0.97 / 0.942 | 0.98 / 0.939 |
| continuous | U(0.3,3)/Σ | weighted Greenwood | 1.00 / 0.945 | 0.98 / 0.946 | 1.00 / 0.944 |
| ceiling 0.5 | none | package | 1.01 / 0.945 | 0.98 / 0.945 | 1.00 / 0.947 |
| ceiling 0.5 | none | survfit robust | 1.01 / 0.945 | 0.98 / 0.945 | 1.00 / 0.947 |
| ceiling 0.5 | none | documented | 0.98 / 0.940 | 0.95 / 0.937 | 0.97 / 0.938 |
| ceiling 0.5 | none | weighted Greenwood | 1.01 / 0.945 | 0.98 / 0.945 | 1.00 / 0.947 |
| ceiling 0.5 | U(0.3,3) | package | **2.31 / 0.999** | **1.97 / 1.000** | **1.71 / 1.000** |
| ceiling 0.5 | U(0.3,3) | survfit robust | 1.00 / 0.942 | 0.98 / 0.946 | 0.99 / 0.944 |
| ceiling 0.5 | U(0.3,3) | documented | 0.97 / 0.934 | 0.95 / 0.940 | 0.96 / 0.934 |
| ceiling 0.5 | U(0.3,3) | weighted Greenwood | 1.00 / 0.941 | 0.98 / 0.944 | 1.00 / 0.943 |
| ceiling 0.5 | U(0.3,3)/Σ | package | **0.88 / 0.886** | **2.18 / 0.993** | **4.32 / 1.000** |
| ceiling 0.5 | U(0.3,3)/Σ | survfit robust | 1.00 / 0.942 | 0.98 / 0.946 | 0.99 / 0.944 |
| ceiling 0.5 | U(0.3,3)/Σ | documented | 0.97 / 0.934 | 0.95 / 0.940 | 0.96 / 0.934 |
| ceiling 0.5 | U(0.3,3)/Σ | weighted Greenwood | 1.00 / 0.941 | 0.98 / 0.944 | 1.00 / 0.943 |

Reading:
- Weighted Greenwood is calibrated in all six cells: SE/SD 0.98–1.01, coverage 0.941–0.948. It is invariant to weight scale (the `U(0.3,3)` and `U(0.3,3)/Σ` rows are identical) and, unweighted, equals the package's Greenwood SE and `survfit` robust.
- Its mean SE differs from `survfit` robust by at most 1e-4 in every cell, including tied times; the documented formula runs slightly small with ties (SE/SD 0.95–0.97), as in Stage 1.
- The package estimator's rows are unchanged from Stage 1: too wide with ties and weights (SE/SD 1.71–2.31), and wrong in both directions with normalised weights (SE/SD 0.45–4.32).
