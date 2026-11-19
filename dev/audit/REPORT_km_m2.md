# weightedsurv: M2, weighted-SE evaluation when weights are related to the hazard

**Top line.** Branch `km-fixes-0.2.0`. START_SHA is `69398e972c2b60d01bfd2901bfdb945faa9845db`, and the package code is the Stage-2b code (no change under `R/`). The final SHA is the commit that adds this report; its parent is `52c064e`, the task-doc commit. Gates G0 and G1 passed. G2 (task section 5) is run after this report is committed. The full run took 0.4 min wall-clock on 12 cores, against a task estimate of 10–15 min and a sanity projection of 0.9 min. All replicates were run: 2,000 in each of the 10 arm-cells, and the first 500 of each for `pkg_resampling`. Non-finite SEs: 0. Files: `dev/audit/m2/km_m2_montecarlo.R`, `km_m2_results.csv` (150 rows), `km_m2_stdout.txt`, `km_m2_sanity.txt` (G1 run), and `km_m2_replicates.rds` (per-replicate values, used by `--check`).

Execution notes (findings):
- **G1, on 200 replicates of cell A.** Bias was at most 0.003. The largest difference from `survfit` was 3.3e-16. SE/SD was 0.98–1.11 for `pkg`, `robust`, `ij` and `xie_liu`. The projection was 0.9 min.
- **Formula cross-checks, on cell-D dataset r = 1 before the main run.**
  - `ij` vs `robust`: the largest relative difference is 4.4e-16 without ties and 5.7e-16 with times rounded up to 0.5.
  - `xie_liu` vs `adjustedCurves::adjustedsurv(method = "iptw_km", conf_int = TRUE)` (adjustedCurves 0.11.4, unstabilised weights): 3.1e-16 in arm 0 and 6.1e-16 in arm 1. They agree, so no separate `adjustedCurves` column was kept.
- **`survfit`'s time-fixing.** By default (`timefix = TRUE`, through `aeqSurv`), `survfit()` merges observed times that differ by less than about 1e-8 relative. The package does not merge them. This happened in 30 of the 20,000 arm-cell × replicate datasets; the count per arm-cell is the last column of Table 3.
  - In 17 arm-cell × replicate × time values, the point estimate differs from default `survfit(weights = w)` by up to 3.8e-05. In those replicates `ij` differs from `robust` by up to 1.1e-04 relative.
  - Against `survfit(weights = w, timefix = FALSE)`, the point estimate agrees in every replicate, arm-cell and time to 2.3e-15. In replicates without merged times, `ij` agrees with `robust` to 8.1e-15.
  - The script records both comparisons. `robust` is kept exactly as the task defines it (default `survfit`).
- **`--check`.** The first 20 replicates of cell B, re-run serially, are identical to the stored ones (largest difference 0).

---

## Table 1 — t = 10: SE/SD (coverage of S ± 1.96 SE)

| arm-cell | pkg | pkg_resampling | robust | ij | xie_liu |
|---|---|---|---|---|---|
| A | 1.011 (0.943) | 1.008 (0.964) | 1.010 (0.942) | 1.010 (0.942) | 1.011 (0.942) |
| B | 1.086 (0.970) | 1.079 (0.966) | 1.005 (0.954) | 1.005 (0.954) | 0.998 (0.951) |
| C | 0.967 (0.935) | 0.964 (0.940) | 1.012 (0.946) | 1.012 (0.946) | 1.050 (0.959) |
| Bt | 1.071 (0.965) | 1.038 (0.958) | 0.992 (0.950) | 0.992 (0.950) | 0.986 (0.949) |
| D0 | 1.077 (0.964) | 1.073 (0.968) | 1.002 (0.950) | 1.002 (0.950) | 0.947 (0.939) |
| D1 | 0.856 (0.906) | 0.852 (0.892) | 0.966 (0.938) | 0.966 (0.938) | 1.062 (0.961) |
| E0 | 1.159 (0.975) | 1.154 (0.968) | 1.078 (0.966) | 1.078 (0.966) | 1.017 (0.956) |
| E1 | 0.936 (0.932) | 0.932 (0.956) | 1.058 (0.960) | 1.058 (0.960) | 1.164 (0.978) |
| F0 | 1.104 (0.963) | 1.067 (0.958) | 0.946 (0.939) | 0.946 (0.939) | 0.822 (0.899) |
| F1 | 0.682 (0.817) | 0.691 (0.828) | 0.894 (0.925) | 0.894 (0.925) | 1.064 (0.947) |

**Monte Carlo error.** At 2,000 replicates it is about ±0.016 on SE/SD and ±0.005 on coverage (one SD).
- **`pkg_resampling`.** Its mean SE and coverage come from the first 500 replicates, divided by the SD over all 2,000. Its coverage error is therefore about ±0.010.
- **Larger errors for low ESS.** The error is larger where the effective-sample-size ratio is low: F1's is 0.34, and F0's is 0.47.

## Table 2 — t = 5 and t = 20: SE/SD

| arm-cell | pkg t=5 | pkg_resampling t=5 | robust t=5 | ij t=5 | xie_liu t=5 | pkg t=20 | pkg_resampling t=20 | robust t=20 | ij t=20 | xie_liu t=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 1.012 | 1.009 | 1.012 | 1.012 | 1.012 | 1.045 | 1.034 | 1.044 | 1.044 | 1.046 |
| B | 1.020 | 1.017 | 0.984 | 0.984 | 0.952 | 1.154 | 1.139 | 0.993 | 0.993 | 1.028 |
| C | 0.956 | 0.955 | 0.989 | 0.989 | 1.056 | 0.979 | 0.973 | 1.020 | 1.020 | 1.032 |
| Bt | 0.998 | 0.968 | 0.962 | 0.962 | 0.933 | 1.165 | 1.127 | 1.003 | 1.003 | 1.041 |
| D0 | 1.050 | 1.043 | 0.998 | 0.998 | 0.902 | 1.075 | 1.073 | 0.993 | 0.993 | 0.977 |
| D1 | 0.915 | 0.911 | 0.988 | 0.988 | 1.144 | 0.839 | 0.832 | 0.981 | 0.981 | 1.023 |
| E0 | 1.090 | 1.084 | 1.034 | 1.034 | 0.933 | 1.120 | 1.115 | 1.035 | 1.035 | 1.018 |
| E1 | 0.979 | 0.979 | 1.056 | 1.056 | 1.222 | 0.877 | 0.875 | 1.032 | 1.032 | 1.077 |
| F0 | 1.017 | 0.991 | 0.918 | 0.918 | 0.745 | 1.165 | 1.137 | 0.977 | 0.977 | 0.923 |
| F1 | 0.731 | 0.739 | 0.888 | 0.888 | 1.172 | 0.625 | 0.604 | 0.845 | 0.845 | 0.926 |

Same Monte Carlo error as Table 1. All coverages are in `km_m2_results.csv`. The lowest are F1 `pkg` at 0.860 / 0.778 (t = 5 / 20) and F1 `pkg_resampling` at 0.854 / 0.758.

## Table 3 — checks

| arm-cell | truth t = 5 / 10 / 20 | bias of point estimate t = 5 / 10 / 20 | ESS ratio | max rel. diff `ij` vs `robust`, replicates without merged times | same, all replicates | replicates with times merged by `survfit` | `adjustedCurves` cross-check |
|---|---|---|---|---|---|---|---|
| A | 0.6065 / 0.3679 / 0.1353 | +0.0002 / +0.0001 / −0.0003 | 0.800 | 5.6e-15 | 2.2e-05 | 3 | n/a (one group) |
| B | 0.5490 / 0.3190 / 0.1293 | +0.0009 / +0.0004 / +0.0003 | 0.800 | 5.4e-15 | 5.4e-15 | 6 | n/a |
| C | 0.7022 / 0.5107 / 0.2884 | −0.0003 / +0.0003 / −0.0001 | 0.800 | 5.0e-15 | 1.1e-04 | 3 | n/a |
| Bt | 0.5490 / 0.3190 / 0.1293 | +0.0023 / +0.0038 / +0.0028 | 0.800 | 6.4e-15 | 6.4e-15 | 0 | n/a |
| D0 | 0.6380 / 0.4408 / 0.2397 | +0.0008 / +0.0008 / +0.0005 | 0.835 | 5.9e-15 | 2.6e-05 | 7 | 3.1e-16 (dataset r = 1) |
| D1 | 0.7316 / 0.5587 / 0.3526 | −0.0010 / −0.0007 / −0.0016 | 0.686 | 4.1e-15 | 8.0e-06 | 3 | 6.1e-16 (dataset r = 1) |
| E0 | 0.6380 / 0.4408 / 0.2397 | +0.0009 / +0.0009 / +0.0002 | 0.831 | 8.1e-15 | 8.1e-15 | 3 | not run (cross-check is on D) |
| E1 | 0.7316 / 0.5587 / 0.3526 | −0.0008 / −0.0008 / −0.0012 | 0.683 | 4.0e-15 | 4.0e-15 | 1 | not run |
| F0 | 0.6380 / 0.4408 / 0.2397 | +0.0024 / +0.0019 / +0.0013 | 0.469 | 6.9e-15 | 1.8e-05 | 3 | not run |
| F1 | 0.7316 / 0.5587 / 0.3526 | −0.0023 / −0.0050 / −0.0069 | 0.342 | 5.0e-15 | 5.0e-15 | 1 | not run |

- **Bias scale.** The Monte Carlo SE of the bias is the SD of the estimate divided by √2000, about 0.0004–0.0016 here.
- **Real biases.** Two are beyond that: Bt at t = 10 (+0.0038), where observed times are rounded up, and F1 at t = 20 (−0.0069).
- **ESS ratio.** It does not depend on t.

## Reading

- **Weights unrelated to the hazard (A).** All five estimators agree and are calibrated: SE/SD 1.01–1.05, coverage 0.94–0.96.
- **`robust` and `ij`.** They are the same estimator to rounding. With true weights (A, B, C, Bt, D) they stay within 0.96–1.05 SE/SD at every time.
- **`pkg` with weights related to the hazard.** As expected, it separates from the others in the predicted direction.
  - Too large where high weights go with high hazard: B 1.09, Bt 1.07 and D0 1.08 at t = 10; up to 1.17 at t = 20.
  - Too small where high weights go with low hazard: C 0.97 and D1 0.86 at t = 10 (D1 0.84 at t = 20, coverage 0.90).
- **`pkg_resampling`.** It tracks `pkg` within about 0.04 everywhere.
- **`xie_liu`.** It is not consistently on either side. It is near 1 in B, but 1.05–1.06 in C and D1, 0.90–0.95 in D0 at t = 5 and 10, and 1.14 in D1 at t = 5.
- **Estimated weights (E against D).** Every estimator's SE/SD rises by 0.05–0.10, so `robust`/`ij` become conservative (1.06–1.08 at t = 10). `pkg` moves from 1.08 / 0.86 (D0 / D1) to 1.16 / 0.94.
- **Extreme weights (F against D; ESS ratio 0.34 in the treated arm).** All estimators fall in the treated arm at t = 10:
  - `pkg` 0.68 (coverage 0.82), `robust`/`ij` 0.89 (0.93), `xie_liu` 1.06 (0.95).
  - In the control arm, `xie_liu` drops to 0.82 (coverage 0.90). The point estimate also carries a bias of up to −0.007 in F1.
