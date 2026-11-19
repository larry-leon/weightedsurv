# weightedsurv 0.2.0: Stage 2c on branch `km-fixes-0.2.0`

**Top line.**
- **Branch and SHAs.** Branch `km-fixes-0.2.0`. START_SHA is `2d7a9b908579a766ec9b409d4431689c11393132` and MASTER_SHA is `166919ec7128f165fb380847f549fa26ca1db985`. The final SHA is the commit that adds this report together with the verification artifacts; its parent is `7ab4171` (docs), and no commit follows it.
- **Gates.** G0, G1 and G2 passed. The G2 output is in section H.
- **Items.** R1 and R2 are DONE; none is NOT DONE.
- **Evaluations.** All three ran (kickoff: `M2b: GO`, `M3: GO`, `M4: GO`), on 12 cores:

  | run | wall-clock | task estimate |
  |---|---|---|
  | M2b | 0.4 min | under 1 min |
  | M3 | 4.0 min (2.2 min "before" + 1.8 min "after") | 10–20 min |
  | M4 | 4.5 min | under 10 min |

  No run came near the 60-min stop.
- **Constraints kept.** Nothing was merged or pushed. Nothing under `vignettes/` was edited. `R/weightedCox.R` was not touched. Installs went only into temporary libraries under the session scratchpad.

Execution notes (findings, not package results):
- **Baseline.** `stage2c/golden_s2b.rds` (6,569 fields, 8 datasets, 164 recorded calls, none of which raised a warning) and `stage2c/api_s2b.rds` (46 exports) were made at `a09b5b9`, the task-doc commit. At that point `git diff --quiet START_SHA -- R/ NAMESPACE` was true. The "before" worktree of START_SHA lived in the scratchpad and was used by M3 and Step V.
- **Golden script.** `stage2c/golden.R` is a copy of `stage2b/golden.R` with the same datasets and calls. It changes three things:
  - Whitelists: R1 and R2.
  - The contract check: every moved analytic `sig2_surv0/1/P` is now compared with `survfit(Surv(aeqSurv) ~ 1, weights, robust = TRUE)`, at 1e-8 relative on the SE scale.
  - A new check that the RMST point estimate (`rmst_maxtau_ci[1]`) does not move.
- **API script.** `stage2c/api.R` accepts no API difference at all.
- **Implementation.** R1 and R2 share one internal helper set in `R/internal_helpers.R` (`ws_km_psi_parts()` and three small companions, all `@noRd`). It computes, in O(n + m) per draw:
  - V(t) = Σ w_i² ψ_i(t)²;
  - M(t, b) = Σ G_ib w_i ψ_i(t);

  from the decomposition ψ_i(g_j) = D_i 1(k_i ≤ j) a[k_i] − H[min(j, k_i)]. The scope switch is `length(unique(W)) > 1`, the test the code already used for the `draws > 0` path. An unsorted `at_points` is handled by computing on the sorted grid and mapping back.
- **R CMD check.** The first attempt stopped at `R CMD build` because pandoc was not on PATH. The second hung for 10 min at "checking use of S3 registration": the sub-process was loading `tcltk` with `DISPLAY` set to XQuartz. The check that counts was run as in Stage 2: `RSTUDIO_PANDOC` set to RStudio's pandoc and `DISPLAY` unset.
- **Golden calls with weights and `draws > 0`.** `KM_diff_weighted_draws` (draws = 200) moves on the five datasets whose `w` is non-constant. The R2 whitelist covers it.

---

## A. Items

| ID | status | commit | files touched | golden fields that moved | tests added or replaced |
|---|---|---|---|---|---|
| R1 | DONE | `6a325a2` | `R/df_counting_helpers.R` (`calculate_risk_event_counts()`, analytic path for unequal weights); `R/internal_helpers.R` (new internal helpers) | 159 fields, all whitelisted. `sig2_surv0`/`sig2_surv1`/`sig2_survP` in `df_counting_weighted` (5 non-constant-`w` datasets) and `df_counting_wv` (8). Weighted `quantile_results` `lower`/`upper` (8 fields). `sig2_surv0`/`sig2_surv1`/`sig2_dhat`/`lower`/`upper` in `KM_diff_wv` and `KM_diff_wv_band`. `Zdhat_star`/`c_alpha_band`/`sb_lower`/`sb_upper` in `KM_diff_wv_band` (the standardisation changed). All 71 moved analytic `sig2_surv*` equal `survfit` robust; the largest relative SE difference is 1.5e-14. Output: `stage2c/golden_compare_R1.txt` | **Added:** `test-r1-robust-variance.R` (7 blocks). **Replaced** (encoded the Stage-2b weighted Greenwood formula; removed from `test-m-weighted-greenwood.R`): "the multiplier equals the formula computed independently, with and without ties" and "df_counting() and KM_diff() report the formula; the survival estimate is unchanged". The helper `wg_ref()` went with them. |
| R2 | DONE | `5784aa1` | `R/df_counting_helpers.R` (`draws > 0` path); `R/km_wlr_calculations_helpers.R` (`resampling_survival()`) | 272 in total (R1 + R2), all whitelisted. R2 adds 113: `surv0_star`/`surv1_star`/`dhat_star` in `KM_diff_wv_band` (8 each); `sig2_rmst_time`, `rmst_time_lower/upper`, `rmst_maxtau_ci` (limits only; the point estimate is checked unchanged), `rmst_text`, `c_alpha_band`, `rmst_time_sb_lower/upper` in `rmst_band_wv` (8 each); `sig2_surv0/1`, `sig2_dhat`, `lower`, `upper` in `KM_diff_weighted_draws` (the 5 non-constant-`w` datasets). Nothing moved in any unweighted or constant-weight call, and no log-rank, Cox, count, risk-table, survival, `dhat` or RMST-estimate field moved. Output: `stage2c/golden_compare_R2.txt` | **Added:** `test-r2-weighted-resampling.R` (5 blocks). **Replaced:** "resampling path (draws > 0, unequal weights) is the one used before and is unchanged" in `test-m-weighted-greenwood.R`, which required the old weighted resampling to stay identical. Its equal-weights half is kept as the block "draws > 0 with equal weights takes the formula, as before". |
| docs | Step 3 | `7ab4171` | `NEWS.md` (M bullet rewritten as "M, R1"; new R2 bullet); roxygen in `R/df_counting_helpers.R`, `R/df_counting.R`, `R/km_wlr_calculations_helpers.R`, `R/rmst_calculations.R`; `man/calculate_risk_event_counts.Rd`, `KM_diff.Rd`, `df_counting.Rd`, `resampling_survival.Rd`, `cumulative_rmst_bands.Rd` | n/a | n/a |

Other commits: `a09b5b9` (task doc), `421caf7` (baseline and copied scripts), and the final commit after `7ab4171`, which adds the verification artifacts and this report.

**Suite.** 16 files, 71 `test_that` blocks, 613 expectations, 0 failures (inside `R CMD check`). Stage 2b had 61 blocks and 560 expectations.

**Tests unchanged.** The other M-test blocks still pass unchanged: unit weights equal Greenwood and `survfit`; the SE is scale-free; constant weights give the unit-weight SE; an empty group gives zeros.

**Error and warning messages.** No error or warning message was added.

## B. R1 and R2 in numbers

`stage2c/report_numbers.R` → `report_numbers.txt`. Golden datasets, weights `wv`.

**B1. Analytic SE.** `KM_diff(draws = 0)` at the grid points at the 25th / 50th / 75th percentile of `at_points`. Stage 2b values are from `golden_s2b.rds`. The ratio is to `survfit(Surv(aeqSurv) ~ 1, weights = wv, robust = TRUE)`.

| dataset | arm | t at 25 / 50 / 75% | SE 2b | SE 2c | ratio to robust 2b | ratio to robust 2c |
|---|---|---|---|---|---|---|
| veteran | 0 | 25 / 80 / 144 | 0.0465 / 0.0633 / 0.0647 | 0.0479 / 0.0644 / 0.0653 | 0.969 / 0.983 / 0.990 | 1.000 / 1.000 / 1.000 |
| veteran | 1 | 25 / 80 / 144 | 0.0560 / 0.0685 / 0.0621 | 0.0566 / 0.0683 / 0.0669 | 0.990 / 1.003 / 0.929 | 1.000 / 1.000 / 1.000 |
| gbsg | 0 | 18.6 / 35.5 / 55.4 | 0.0219 / 0.0275 / 0.0306 | 0.0219 / 0.0274 / 0.0306 | 0.997 / 1.006 / 0.999 | 1.000 / 1.000 / 1.000 |
| gbsg | 1 | 18.6 / 35.5 / 55.4 | 0.0247 / 0.0327 / 0.0363 | 0.0248 / 0.0328 / 0.0370 | 0.997 / 0.998 / 0.982 | 1.000 / 1.000 / 1.000 |
| gbsg_er_le0 | 0 | 14.7 / 25.6 / 48 | 0.0586 / 0.0643 / 0.0863 | 0.0579 / 0.0658 / 0.0823 | 1.012 / 0.977 / 1.049 | 1.000 / 1.000 / 1.000 |
| gbsg_er_le0 | 1 | 14.7 / 25.6 / 48 | 0.1013 / 0.1169 / 0.1091 | 0.1023 / 0.1166 / 0.1105 | 0.991 / 1.002 / 0.988 | 1.000 / 1.000 / 1.000 |
| gbsg_er_gt0 | 0 | 19.6 / 35.8 / 56 | 0.0253 / 0.0302 / 0.0337 | 0.0251 / 0.0299 / 0.0337 | 1.010 / 1.008 / 1.002 | 1.000 / 1.000 / 1.000 |
| gbsg_er_gt0 | 1 | 19.6 / 35.8 / 56 | 0.0262 / 0.0343 / 0.0402 | 0.0262 / 0.0341 / 0.0404 | 1.000 / 1.004 / 0.996 | 1.000 / 1.000 / 1.000 |
| sim_cont | 0 | 2.55 / 5.35 / 11.3 | 0.0392 / 0.0559 / 0.0613 | 0.0392 / 0.0559 / 0.0609 | 1.000 / 1.001 / 1.006 | 1.000 / 1.000 / 1.000 |
| sim_cont | 1 | 2.55 / 5.35 / 11.3 | 0.0478 / 0.0543 / 0.0625 | 0.0474 / 0.0544 / 0.0642 | 1.007 / 0.998 / 0.974 | 1.000 / 1.000 / 1.000 |
| sim_ties | 0 | 3.5 / 6.5 / 13 | 0.0505 / 0.0582 / 0.0666 | 0.0520 / 0.0602 / 0.0622 | 0.971 / 0.967 / 1.072 | 1.000 / 1.000 / 1.000 |
| sim_ties | 1 | 3.5 / 6.5 / 13 | 0.0409 / 0.0476 / 0.0623 | 0.0407 / 0.0485 / 0.0616 | 1.004 / 0.981 / 1.011 | 1.000 / 1.000 / 1.000 |
| sim_weighted | 0 | 3.37 / 6.72 / 13.2 | 0.0395 / 0.0464 / 0.0462 | 0.0396 / 0.0466 / 0.0484 | 0.998 / 0.995 / 0.955 | 1.000 / 1.000 / 1.000 |
| sim_weighted | 1 | 3.37 / 6.72 / 13.2 | 0.0271 / 0.0458 / 0.0525 | 0.0271 / 0.0453 / 0.0520 | 1.001 / 1.013 / 1.010 | 1.000 / 1.000 / 1.000 |
| sim_strata | 0 | 2.59 / 6.36 / 13 | 0.0359 / 0.0472 / 0.0523 | 0.0363 / 0.0470 / 0.0511 | 0.988 / 1.005 / 1.024 | 1.000 / 1.000 / 1.000 |
| sim_strata | 1 | 2.59 / 6.36 / 13 | 0.0304 / 0.0473 / 0.0503 | 0.0305 / 0.0465 / 0.0505 | 0.997 / 1.019 / 0.997 | 1.000 / 1.000 / 1.000 |

**B2. Bands.** `KM_diff(draws.band = 200)` and `cumulative_rmst_bands()` (defaults) on that fit. RMST values are at the last grid time.

| dataset | KM band c_alpha 2b → 2c | RMST (unchanged) | RMST pointwise limits 2b → 2c | RMST band c_alpha 2b → 2c | RMST band limits 2b → 2c |
|---|---|---|---|---|---|
| veteran | 3.193 → 3.118 | 40.535 | (−24.11, 105.17) → (−37.39, 118.46) | 2.598 → 2.534 | (−45.16, 126.23) → (−60.22, 141.29) |
| gbsg | 2.921 → 2.956 | 10.799 | (5.15, 16.45) → (5.11, 16.49) | 2.916 → 2.947 | (2.39, 19.21) → (2.25, 19.35) |
| gbsg_er_le0 | 2.553 → 2.647 | −14.542 | (−30.49, 1.41) → (−30.67, 1.59) | 2.607 → 2.655 | (−35.76, 6.68) → (−36.39, 7.30) |
| gbsg_er_gt0 | 3.043 → 3.011 | 11.232 | (5.32, 17.14) → (5.46, 17.01) | 2.826 → 2.890 | (2.71, 19.75) → (2.72, 19.75) |
| sim_cont | 2.867 → 3.124 | 2.575 | (−0.41, 5.56) → (−1.18, 6.33) | 2.808 → 2.634 | (−1.70, 6.85) → (−2.47, 7.62) |
| sim_ties | 2.706 → 2.839 | 3.983 | (0.66, 7.30) → (0.61, 7.36) | 2.574 → 2.557 | (−0.38, 8.34) → (−0.42, 8.38) |
| sim_weighted | 3.030 → 3.136 | 3.202 | (0.36, 6.05) → (0.54, 5.86) | 2.598 → 2.643 | (−0.57, 6.97) → (−0.39, 6.79) |
| sim_strata | 2.989 → 3.080 | 2.056 | (−0.50, 4.61) → (−0.63, 4.74) | 2.707 → 2.642 | (−1.47, 5.58) → (−1.57, 5.68) |

Reading:
- **Analytic SE.** The weighted SE now equals `survfit` robust to the printed precision in every row (golden contract: 1.5e-14). Stage 2b's weighted Greenwood ratios to it were 0.929–1.072 on these datasets.
- **Band critical values.** The KM band critical value moves by −0.08 to +0.26 because both the resampled curves (R2) and the SEs that standardise them (R1) changed.
- **RMST.** The RMST estimate is unchanged. The RMST pointwise limits widen or narrow with the new resampling variance; the largest change is veteran, whose last grid time is 999 days with a long sparse tail.

## C. M2b

`stage2c/km_m2b_montecarlo.R` is `m2/km_m2_montecarlo.R` with only its output path changed (`diff`: one line). It ran `--full` in 0.4 min. Results: `stage2c/km_m2_results.csv`, `km_m2_replicates.rds`, `km_m2b_stdout.txt`. The requirements are checked by `km_m2b_check.R` (→ `km_m2b_check.txt`, "M2b: PASS").
- **`pkg` equals `ij`.** The largest relative difference is 6.6e-14, over the 59,910 arm-cell × replicate × time values in which `survfit` merged no times (requirement: < 1e-8). Against `robust` it is 6.9e-14.
- **Resampling vs analytic.** On the first 500 replicates of each arm-cell, `mean(pkg_resampling SE) / mean(pkg SE)` is 0.9948–1.0034 over all 10 arm-cells × 3 times (requirement: 0.97–1.03).

**Table 1 of M2 again, t = 10: SE/SD (coverage).** M2 = Stage-2b code, M2b = Stage-2c code. `robust` is the same column in both runs.

| arm-cell | pkg M2 | pkg_resampling M2 | **pkg M2b** | **pkg_resampling M2b** | robust |
|---|---|---|---|---|---|
| A | 1.011 (0.943) | 1.008 (0.964) | 1.010 (0.942) | 1.015 (0.962) | 1.010 (0.942) |
| B | 1.086 (0.970) | 1.079 (0.966) | 1.005 (0.954) | 1.000 (0.948) | 1.005 (0.954) |
| C | 0.967 (0.935) | 0.964 (0.940) | 1.012 (0.946) | 1.012 (0.948) | 1.012 (0.946) |
| Bt | 1.071 (0.965) | 1.038 (0.958) | 0.992 (0.950) | 0.989 (0.948) | 0.992 (0.950) |
| D0 | 1.077 (0.964) | 1.073 (0.968) | 1.002 (0.950) | 1.003 (0.952) | 1.002 (0.950) |
| D1 | 0.856 (0.906) | 0.852 (0.892) | 0.966 (0.938) | 0.968 (0.934) | 0.966 (0.938) |
| E0 | 1.159 (0.975) | 1.154 (0.968) | 1.078 (0.966) | 1.080 (0.960) | 1.078 (0.966) |
| E1 | 0.936 (0.932) | 0.932 (0.956) | 1.058 (0.960) | 1.062 (0.978) | 1.058 (0.960) |
| F0 | 1.104 (0.963) | 1.067 (0.958) | 0.946 (0.939) | 0.939 (0.940) | 0.946 (0.939) |
| F1 | 0.682 (0.817) | 0.691 (0.828) | 0.894 (0.925) | 0.905 (0.932) | 0.894 (0.925) |

t = 5 and t = 20 are in `km_m2b_check.txt`. At t = 20 the M2b `pkg` SE/SD is 0.845 (F1) to 1.044 (A), against 0.625–1.165 in M2. The lowest coverage is F1: 0.889 for `pkg`, 0.874 for `pkg_resampling`.

## D. M3: simultaneous-band coverage

`stage2c/km_m3_bands.R`: 4 cells × 2 versions × 1,000 replicates, N = 800, `draws.band = draws_sb = 500`. Results: `km_m3_results.csv`, `km_m3_before.rds`, `km_m3_after.rds`, `km_m3_stdout.txt`. Truth is from numerical integration over X (step 0.01).

**How the table was computed:**
- **Pointwise read-off.** The estimate at t is the value at the last grid point ≤ t, compared with δ(t). The cumulative RMST at grid time t_j is compared with R(t_j).
- **Two kinds of simultaneous coverage.** "All" is over every grid point. "SE>0" is over the grid points where the SE is positive. Where the SE is 0 the band has width 0, and it cannot cover a true value that is not exactly the estimate. That happens on 0.26–0.32 points per replicate for δ, and on about 1.05 points for the RMST curve.
- **RMST all-grid coverage is 0 by construction.** The RMST curve starts at 0 with SE 0 at the first grid point, while R(t_1) > 0. So its all-grid coverage is 0.000 in every cell, and the RMST columns below show the SE>0 version only.
- **Monte Carlo error.** About ±0.007 on coverage.

| cell | version | δ pointwise t = 5 / 10 / 20 | δ band all | δ band all, SE>0 | δ band [2, 25] | RMST pointwise t = 20 | RMST band all, SE>0 | RMST band [2, 25] | mean c_alpha δ | mean c_alpha RMST | mean half-width δ band t = 10 | SE/SD of δ̂ t = 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U | before | 0.955 / 0.952 / 0.952 | 0.727 | 0.924 | 0.977 | 0.948 | 0.908 | 0.984 | 3.139 | 2.741 | 0.1186 | 1.001 |
| U | after | 0.955 / 0.952 / 0.952 | 0.727 | 0.924 | 0.977 | 0.948 | 0.908 | 0.984 | 3.139 | 2.741 | 0.1186 | 1.001 |
| D | before | 0.948 / 0.933 / 0.911 | 0.662 | 0.837 | 0.943 | 0.931 | 0.853 | 0.966 | 3.111 | 2.736 | 0.1298 | 0.930 |
| D | after | 0.950 / 0.947 / 0.937 | 0.688 | 0.867 | 0.960 | 0.942 | 0.860 | 0.973 | 3.111 | 2.748 | 0.1354 | 0.970 |
| E | before | 0.963 / 0.960 / 0.948 | 0.688 | 0.869 | 0.975 | 0.961 | 0.857 | 0.991 | 3.108 | 2.734 | 0.1301 | 1.073 |
| E | after | 0.969 / 0.970 / 0.965 | 0.705 | 0.892 | 0.988 | 0.966 | 0.867 | 0.994 | 3.107 | 2.746 | 0.1359 | 1.120 |
| F | before | 0.933 / 0.912 / 0.840 | 0.502 | 0.650 | 0.854 | 0.877 | 0.700 | 0.916 | 3.027 | 2.717 | 0.1730 | 0.855 |
| F | after | 0.936 / 0.932 / 0.913 | 0.569 | 0.737 | 0.923 | 0.911 | 0.745 | 0.949 | 3.044 | 2.722 | 0.1868 | 0.918 |

Required check: in cell U, "before" and "after" are identical in all 1,000 replicates. The digest of the full `KM_diff()` and `cumulative_rmst_bands()` output is identical, and so is every recorded number. In cells D, E and F, all 1,000 replicates differ.

Reading:
- **Cells D, E, F.** Coverage rises in every column. The largest gains are in F: δ pointwise at t = 20 goes 0.840 → 0.913, and the δ band on [2, 25] 0.854 → 0.923.
- **SE/SD of δ̂ at t = 10.** D 0.930 → 0.970, E 1.073 → 1.120, F 0.855 → 0.918.
- **Estimated weights (E).** E is above nominal in both versions, and more so after (band on [2, 25]: 0.988).
- **Extreme weights (F).** F stays below nominal: band on [2, 25] 0.923, RMST pointwise 0.911.
- **Whole grid.** Whole-grid band coverage (SE>0) is below nominal in every cell, including the unweighted U (0.924), because of the grid ends. On [2, 25], U is 0.977.
- **Critical values.** The mean critical values barely move (δ: 3.11 → 3.11 in D, 3.03 → 3.04 in F). The band widens through the SE: half-width at t = 10 is +4% in D and +8% in F.

## E. M4: weighted log-rank and Cox under the null (measurement only)

`stage2c/km_m4_null.R`, HEAD code, 4 cells × 2,000 replicates, N = 800, no treatment effect. Results: `km_m4_results.csv`, `km_m4_replicates.rds`, `km_m4_stdout.txt`.

**What the columns mean:**
- "z" is `df_counting()$z.score` (a non-finite value is set to `NA` by `df_counting()`, with a warning). Its rejection rate is `|z| > 1.96` over the finite values.
- "rhogamma" is `cox_rhogamma(rho = 0, gamma = 0, draws = 0)`: `fit$bhat` and `fit$sig_bhat_asy`.
- The `coxph` columns do not depend on the package scaling, since the reference is fitted once per replicate on the original weights.

| cell | weights | z mean | z SD | z reject | z non-finite | rhogamma: SD of b̂ / mean SE (ratio) | rhogamma reject | coxph robust: z SD, reject; SD of b̂ / mean SE | coxph model-based: z SD, reject; mean SE |
|---|---|---|---|---|---|---|---|---|---|
| U | none | 0.001 | 0.978 | 0.046 | 0 | 0.0843 / 0.0862 (1.023) | 0.045 | 0.978, 0.046; 0.0843 / 0.0861 | 0.977, 0.045; 0.0862 |
| D | w | −0.080 | 1.814 | 0.276 | 0 | 0.1105 / 0.0609 (0.551) | 0.276 | 1.012, 0.054; 0.1105 / 0.1095 | 1.808, 0.275; 0.0610 |
| D | 10·w | −0.253 | 5.741 | 0.739 | 0 | 0.1105 / 0.0193 (0.174) | 0.739 | (same) | (same) |
| D | w / Σw | ≈0 (finite ones) | 3.4e-08 | 0.000 | 1,286 | 0.1105 / 3.452 (31.2) | 0.000 | (same) | (same) |
| E | w | −0.099 | 1.589 | 0.215 | 0 | 0.0967 / 0.0609 (0.629) | 0.215 | 0.890, 0.030; 0.0967 / 0.1093 | 1.585, 0.214; 0.0610 |
| E | 10·w | −0.314 | 5.031 | 0.702 | 0 | 0.0967 / 0.0192 (0.199) | 0.702 | (same) | (same) |
| E | w / Σw | ≈0 | 2.7e-08 | 0.000 | 1,267 | 0.0967 / 3.446 (35.6) | 0.000 | (same) | (same) |
| F | w | −0.461 | 3.088 | 0.534 | 0 | 0.1875 / 0.0612 (0.326) | 0.534 | 1.185, 0.107; 0.1875 / 0.1603 | 3.054, 0.530; 0.0614 |
| F | 10·w | −1.460 | 9.773 | 0.841 | 0 | 0.1875 / 0.0193 (0.103) | 0.841 | (same) | (same) |
| F | w / Σw | ≈0 | 5.9e-08 | 0.000 | 1,308 | 0.1875 / 3.462 (18.5) | 0.000 | (same) | (same) |

Reading:
- **No weights (U).** The package z and the `cox_rhogamma()` SE are calibrated: z SD 0.978, SE/SD 1.023, rejection 0.045–0.046.
- **Propensity weights.** The z SD is 1.81 (D), 1.59 (E) and 3.09 (F); rejection at nominal 0.05 is 0.28, 0.22 and 0.53. The `cox_rhogamma()` SE is 0.55, 0.63 and 0.33 of the SD of its estimate.
- **Same as `coxph` model-based.** These numbers are essentially those of the model-based `coxph` SE (z SD 1.81 / 1.59 / 3.05). The robust `coxph` z has SD 1.01 / 0.89 / 1.18 and rejection 0.054 / 0.030 / 0.107.
- **Scaling the weights by 10.** This multiplies the z SD by about √10 (D 1.81 → 5.74) and divides the `cox_rhogamma()` SE by about √10. The estimate b̂ does not change.
- **Normalised weights (w / Σw).** The z-score is non-finite in 64–65% of replicates (1,267–1,308 of 2,000); `df_counting()` warns once per such replicate. The finite values are of order 1e-8. The `cox_rhogamma()` SE is about 3.45, 18–36 times the SD of b̂.
- **Point estimate.** b̂ from `cox_rhogamma()` has the same SD as `coxph`'s b̂ in every cell.

## F. Vignette impact inventory (no edits)

Method: `stage2c/vignette_impact.R`. Each version was installed into its own temporary library: START_SHA from the worktree, HEAD from `git archive`. The three files with weighted calls were knitted with `knitr::knit` under each version, so the vignettes' own code was used, nothing was edited, and figures were written as PNG.
- **What was captured.** After every chunk, a hook recorded the displayed quantities of every object that the chunk created or changed. Numbers: `stage2c/vignette_numbers.txt`. Figures: `stage2c/vignette_weighted_figures.pdf`, 15 pages, one per figure whose PNG differs, before | after.
- **Files with no weighted calls.** `weightedkm_methodology.Rmd` and `weightedrmst_methodology.Rmd` have only a hidden setup chunk; the weighted code in the RMST article is a ` ```r ` block that is not evaluated. `weightedcox_wald-vs-logrank_simulations.Rmd` uses only time-dependent weights. `README.md` passes no subject weights.

### F1. Calls that pass weights

| file:line (chunk header) | function(s) | what it shows |
|---|---|---|
| `weightedsurv_examples.Rmd:310` (call 313) | `get_dfcounting(weight.name = "sw.weights")` (Rotterdam) | feeds 320, 329, 624, 702 |
| `weightedsurv_examples.Rmd:320` | `plot_weighted_km()` on the weighted fit | weighted KM curves; median with CI (text); coxph HR (text) |
| `weightedsurv_examples.Rmd:329` | `plot_weighted_km(conf.int = TRUE, show.logrank = TRUE)` | pointwise band; median CI; **subject-weighted log-rank p**; coxph HR |
| `weightedsurv_examples.Rmd:346` (call 351) | `plotKM.band_subgroups(weight.name, draws.band = 1000, show_resamples = TRUE)` | 20 standardised resampled curves; pointwise CIs; simultaneous band |
| `weightedsurv_examples.Rmd:369` (calls 374, 386) | `plotKM.band_subgroups(weight.name)`; `cumulative_rmst_bands(weight.name)` | pointwise CIs and simultaneous band; RMST pointwise and simultaneous band; RMST(τ*) with CI (text) |
| `weightedsurv_examples.Rmd:476` (call 480) | `get_dfcounting(weight.name)` (cross-trial data) | feeds 485, 502 |
| `weightedsurv_examples.Rmd:485` | `plot_weighted_km(conf.int = TRUE, show.logrank = TRUE)` | pointwise band; **subject-weighted log-rank p**; coxph HR |
| `weightedsurv_examples.Rmd:502` | `wlr_cumulative(dfcount_CTC_wtd)` | **cumulative subject-weighted log-rank z** |
| `weightedsurv_examples.Rmd:520` (calls 546, 558) | `plotKM.band_subgroups(weight.name)`; `cumulative_rmst_bands(weight.name)` | as 369, cross-trial data |
| `weightedsurv_examples.Rmd:624` | `cox_rhogamma(dfcount_rotterdam_wtd, draws = 1000)` × 5 schemes | **Cox HR, resampling CI, log-rank z (gt table)** |
| `weightedsurv_examples.Rmd:702` (call 703) | `get_dfcounting(weight.name, check.seKM = TRUE, draws = 0)` | SE comparison figure, package vs `survfit(weights)` |
| `weightedsurv_examples.Rmd:711` (call 712) | same with `draws = 5000` | resampling-SE comparison figure |
| `weightedsurv_examples.Rmd:732` | `plot_weighted_km(dfcount_rotterdam_wtd, conf.int = TRUE)` (the `draws = 5000` fit) | pointwise band from the resampling SE |
| `weightedsurv_examples.Rmd:743` | SE plot of that fit against `adjustedCurves` | SE comparison figure |
| `martingale_resampling_applications.Rmd:196` | `df_counting(weight.name)`; `plot_weighted_km()` | weighted KM curves; median CI; coxph HR |
| `martingale_resampling_applications.Rmd:241` | `plot_weighted_km(conf.int = TRUE)` | pointwise band; median CI |
| `martingale_resampling_applications.Rmd:274` | `plotKM.band_subgroups(weight.name, show_resamples = TRUE)` | resampled curves; pointwise CIs; simultaneous band |
| `martingale_resampling_applications.Rmd:316` (bottom row) | `plotKM.band_subgroups(weight.name)`; `cumulative_rmst_bands(weight.name)` | simultaneous band; RMST bands; RMST(τ*) with CI |
| `martingale_resampling_applications.Rmd:421` | `eval = FALSE` workflow example | not run |
| `weightedcox_methodology.Rmd:315` | `get_dfcounting(weight.name)`; `cox_rhogamma(draws = 1000)` × 2 | **Cox HR and resampling CI (printed)** |
| `weightedrmst_methodology.Rmd:281` | ` ```r ` block, not evaluated | not run |

### F2. Numbers, before → after

From `vignette_numbers.txt`. Survival estimates, median point estimates, log-rank z (5.894 Rotterdam, 3.656 CTC), coxph HR text, the RMST estimate and every `cox_rhogamma()` result are identical before and after.

| chunk | quantity | before | after |
|---|---|---|---|
| examples 310 / 320 / 329 / 624 / 702; martingale 196 / 241; coxmeth 315 (Rotterdam, draws = 0) | treated median 80.39, lower limit | 60.25 | 59.33 |
| same | control median 40.84, CI | (38.21, 47.38) | (38.21, 47.38) |
| same | SE control arm at 24 / 60 / 120 months | 0.0141 / 0.0146 / 0.0140 | 0.0140 / 0.0144 / 0.0139 |
| same | SE treated arm at 24 / 60 / 120 | 0.0306 / 0.0447 / 0.0786 | 0.0318 / 0.0475 / 0.0741 |
| examples 711 / 732 / 743 (Rotterdam, draws = 5000) | control median CI | (38.21, 47.38) | (38.28, 46.75) |
| same | treated lower limit | 60.25 | 59.33 |
| same | SE control / treated at 24 / 60 / 120 | 0.0141 / 0.0145 / 0.0140; 0.0305 / 0.0441 / 0.0728 | 0.0139 / 0.0144 / 0.0139; 0.0317 / 0.0477 / 0.0743 |
| examples 476 / 485 (CTC) | control median 45.44, lower limit | 40.41 | 40.34 |
| same | SE control / treated at 24 / 60 / 120 | 0.0140 / 0.0150 / 0.0151; 0.0523 / 0.0489 / 0.0620 | 0.0141 / 0.0150 / 0.0149; 0.0535 / 0.0570 / 0.0738 |
| examples 346 / 369; martingale 274 / 316 (Rotterdam band) | band critical value | 3.166 | 3.539 |
| same | δ̂ 0.158 / 0.172 / 0.089 at 24 / 60 / 120: SE | 0.0338 / 0.0415 / 0.0783 | 0.0247 / 0.0512 / 0.0659 |
| same | pointwise CI at 60 | (0.090, 0.253) | (0.071, 0.272) |
| same | simultaneous band at 24 / 60 / 120 | (0.051, 0.265) / (0.040, 0.303) / (−0.159, 0.337) | (0.071, 0.246) / (−0.009, 0.353) / (−0.144, 0.322) |
| martingale 316 (and examples 369, same call) | RMST text | RMST(τ*) = 16.8 (8.6, 25) | RMST(τ*) = 16.8 (8, 25.7) |
| same | RMST band critical value; band at τ* = 118.4 | 2.678; (5.60, 28.05) | 2.684; (4.71, 28.94) |
| examples 520 (CTC band) | band critical value | 3.310 | 3.630 |
| same | SE of δ̂ at 24 / 60 / 120 | 0.0542 / 0.0511 / 0.0624 | 0.0626 / 0.0638 / 0.0858 |
| same | simultaneous band at 60 / 120 | (0.025, 0.363) / (−0.074, 0.339) | (−0.038, 0.426) / (−0.179, 0.444) |
| same | RMST text; band at τ* = 77.9 | RMST(τ*) = 8.1 (1.6, 14.7); (−0.62, 16.86) | RMST(τ*) = 8.1 (1, 15.3); (−1.40, 17.63) |

Notes on these numbers:
- **Band chunks use a 20-draw pointwise SE.** `plotKM.band_subgroups()` defaults to `draws = 20`. With weights, its pointwise SE (`sig2_dhat`, which also standardises the band) is therefore the `draws > 0` resampling variance from 20 draws, under both versions. On the Rotterdam call it is 0.70–1.07 times the analytic SE over the grid (median 0.93), and 0.71 at 24 months (`stage2c/band_draws20_check.txt`). That is why the band-chunk SE at 24 months falls (0.0338 → 0.0247) while the arm SEs at 24 months rise. Unweighted band calls are not affected: equal weights take the Greenwood path whatever `draws` is.
- **Sign changes.** In the Rotterdam band the simultaneous lower limit at 60 months changes sign (0.040 → −0.009). In the CTC band it does so at 60 months (0.025 → −0.038). Both pointwise CIs at 60 months still exclude 0.
- **Captured last value only.** In chunk `examples:369` the hook captured only the last value of `temp` and `get_bands`, which belong to the unweighted GBSG panels. The weighted Rotterdam panels of that chunk are the same calls as `martingale:316` (identical arguments), whose numbers are above. Its figure is on the PDF.
- **SE comparison chunks** (`examples:702–768`, `vignette_adjcurves_check.R` → `vignette_adjcurves_check.txt`):
  - Package / `survfit` robust, `draws = 0`: medians 1.0081 / 0.9576 (control / treated) → 1.0000 / 1.0000.
  - The same ratio for the `draws = 5000` fit that the chunks at 732 and 743 plot: 1.0074 / 0.9524 → 0.9928 / 1.0023.
  - Package / `adjustedCurves` (`iptw_km`) at `draws = 5000`: medians 1.026 / 0.827 → 1.008 / 0.877; treated-arm range 0.41–1.02 in both versions.
  - `check.seKM` compares with `survfit(weights = w)`, which already uses the robust variance by default when case weights are given (checked: identical to `robust = TRUE`).

### F3. Figures

`stage2c/vignette_weighted_figures.pdf`, one page per figure that differs, before on the left and after on the right.

**15 figures differ:**
- `weightedsurv_examples.Rmd` chunks at 320, 329, 346, 369, 485, 520, 702, 711, 732, 743;
- `martingale_resampling_applications.Rmd` chunks at 196, 241, 274 (two figures), 316.

**Byte-identical** (not in the PDF):
- examples 70, 203, 213, 233 (unweighted);
- examples 502 (weighted cumulative log-rank);
- examples 578, 606 (unweighted band plots);
- examples 737 (`adjustedCurves` own plot);
- martingale 132;
- coxmeth 305.

The 320 and 196 figures differ only in the median-CI text, because `plot_weighted_km()` defaults to `conf.int = FALSE`.

### F4. Sentences or formulas no longer accurate

| file:line | text | why |
|---|---|---|
| `weightedkm_methodology.Rmd:203–207` (edited in Stage 2b) | "the standard error of the weighted KM estimator is the weighted Greenwood estimator (used whenever no resampling is requested …)", the weighted Greenwood formula, and "With all weights equal this is Greenwood's formula" | for unequal weights the SE is now the robust (infinitesimal-jackknife) variance, equal to `survfit(robust = TRUE)`; the formula shown is no longer the one computed. The sentences on fixed/known weights and scale-invariance remain true |
| `weightedkm_methodology.Rmd:209` | "For settings involving complex weighting … variance estimation via martingale resampling is preferred." | for weighted data the resampling now estimates the same robust variance as the analytic formula (M2b: ratio 0.995–1.003); the contrast the sentence draws no longer exists |
| `weightedkm_methodology.Rmd:234` | "the martingale increments dM_{i,j}(t) are replaced by their observable counterparts multiplied by independent standard normal variates" | does not hold for weighted data: the multipliers now enter through each subject's full influence contribution (event term and compensator, censored subjects included) |
| `weightedkm_methodology.Rmd:242` | "The observable counting process increments dN_{i,j}(t) serve as proxies for the unobservable martingale increments" | as above, for weighted data |
| `weightedkm_methodology.Rmd:251` | "b. Compute the resampled survival difference process … by perturbing the counting process increments." | as above, for weighted data |
| `weightedkm_methodology.Rmd:259` | "In the weighted case, the counting process increments are replaced by their weighted counterparts, and the variance proxy accounts for the propensity-score weights." | describes the replaced events-only weighted perturbation |
| `weightedkm_methodology.Rmd:276` | "Simultaneous confidence bands for the RMST difference curve are obtained by the same martingale resampling approach" | for weighted data the resampled curves are influence-based, not the dN-based process described at 234–255 |
| `weightedrmst_methodology.Rmd:219–223` | "The weighted version of the resampled statistic for group j replaces unweighted counting process …" and the weighted W† formula (dN_{i,j} G_{i,j} / Ȳ^w_j) | not what is computed for weighted KM / RMST curves any more |
| `weightedrmst_methodology.Rmd:225` | "The `weightedsurv` package implements this weighted resampling directly in the `cumulative_rmst_bands` function when a `weight.name` argument is supplied" | it now implements the influence-based resampling |
| `martingale_resampling_applications.Rmd:303` | "The grey traces show individual resampled processes ($W^{\dagger}$)" | W† (line 85) is the events-only weighted log-rank process; for the weighted Rotterdam plot the traces are standardised influence-based KM-difference resamples |
| `martingale_resampling_applications.Rmd:419` | "The resampling distribution is modified accordingly, with the at-risk process Ȳ_j(t) replaced by the weighted version Ȳ_j^w(t)" | for weighted data the resampling is no longer the unweighted one with Ȳ replaced by Ȳ^w |
| `martingale_resampling_applications.Rmd:448` | "The same resampling framework applies naturally to propensity-score weighted analyses" | for KM and RMST curves the weighted resampling is a different construction from the unweighted one |
| `weightedkm_methodology.Rmd:326`; `README.md:231` | "cross-validation of weighted SEs against the `adjustedCurves` and `survfit` packages"; "Validation against the `survival` and `adjustedCurves` packages" | the `survfit` comparison now agrees exactly (`draws = 0`). The `adjustedCurves` comparison (`iptw_km`, a different variance estimator) still does not agree: treated-arm ratios 0.41–1.02 |

The Stage-2b RMST passage (`weightedrmst_methodology.Rmd:90–95`, integration rule) was checked and is still accurate. Its last sentence, "The same rule is applied to every resampled curve", still holds.

Checked and not listed:
- `weightedsurv_examples.Rmd:341`, the resampled approximations: holds, and is unchanged for unweighted data.
- `martingale_resampling_applications.Rmd:263, 408`: qualitative statements about point estimates, which do not move.
- `weightedrmst_methodology.Rmd:112`: says the package uses the delta-method RMST variance. It used the resampling variance in both versions, so this is not caused by this task.
- `weightedkm_methodology.Rmd:252–255`: describe an unstandardised sup and band Δ̂ ± c. The code uses the standardised version in both versions, so this is not caused by this task either.

### F5. Pending scope (subject-weighted log-rank z or `cox_rhogamma()` with weights; not changed by this task)

| file:line | what is displayed | before = after? |
|---|---|---|
| `weightedsurv_examples.Rmd:329` | `plot_weighted_km(show.logrank = TRUE)`: "logrank (1-sided) p = <0.001" from the subject-weighted z = 5.894 | identical |
| `weightedsurv_examples.Rmd:485` | same, cross-trial data, z = 3.656 | identical |
| `weightedsurv_examples.Rmd:502` | `wlr_cumulative(dfcount_CTC_wtd)`: cumulative subject-weighted log-rank z curve | identical (figure byte-identical) |
| `weightedsurv_examples.Rmd:624` | `cox_rhogamma()` with sw.weights, 5 schemes, `draws = 1000`: HR, resampling CI, z (gt table) | identical |
| `weightedcox_methodology.Rmd:315` | `cox_rhogamma()` with sw.weights, FH(0,0) and FH(0,1): HR and resampling CI (printed) | identical |

## G. Seed probe, pressure test, API, CRAN check, forestsearch

**Probe** (`dev/audit/km_audit_probe.R --no-write` → `stage2c/probe_2c.txt`):
- Statuses are unchanged from Stage 2b: OK1–OK3 PASS, I3 the only REPRODUCED, P5 INFO, everything else NOT_REPRODUCED.
- Three lines changed text:
  - OK3: pkg / survfit robust SE was in [0.943, 1.022] and is now in [1.000, 1.000].
  - K4: tied times, ratio at 4 times was 1.00 / 1.02 / 0.97 / 1.01 and is now 1.00 / 1.00 / 1.00 / 1.00.
  - P5: timing only, 0.3 s → 0.4 s.

**Pressure test** (`stage2c/km_pressure_test_s2c.R`, a copy of `stage2b/km_pressure_test_s2b.R` that differs only in its default output directory; 92 s; `km_pressure_results.csv`, `km_pressure_stdout.txt`). Failing cells by family:

| family | checked cells | failing Stage 2b | failing Stage 2c |
|---|---|---|---|
| REF | 1 | 0 | 0 |
| T1 | 111 | 0 | 0 |
| T2 | 17 | 0 | 0 |
| T3 | 288 | 0 | 0 |
| T4 | 21 | 0 | 0 |
| T5 | 9 | 0 | 0 |
| T6 | 24 | 0 | 0 |
| T7 | 18 | 0 | 0 |
| T8 | 73 | 0 | 0 |
| T9 | 1 | 0 | 0 |
| T10 | 5 | 3 | 3 |

No cell changed `n_checks` or `n_fail`. Of the 612 rows, 32 report-only rows changed their value, each explained:
- **T1, "weighted SE ratio pkg/survfit robust (report)"** (15 rows, weights int / unif / ipw / norm / zero10 × three tie levels). Over the design grid (n 2–500, censoring 0–0.95), every row is now median 1, range [1, 1]. Stage 2b had ranges up to [0.683, 1.7] (int) and [0.226, 5.69] (ipw). The package SE is now `survfit` robust's.
- **T1, "weighted SE ratio pkg/documented formula (report)"** (15 rows). The reference is the formula in the 0.1.0 help page, which is not the package's method in either stage. The worst ratio fell, for example from 22.3 to 9.1 for zero10 with heavy ties, and medians stayed at 1.02–1.15. These rows compare against a formula that no longer describes any version.
- **T2.8, SE with draws = 2000 / draws = 0, weighted (report).** The quartiles were 0.921 / 0.978 / 0.987 and are now 1.00 / 1.01 / 1.02 (worst 0.725 → 1.04). This is R2: the resampling SE now estimates the same robust variance.
- **T5, var(resampled surv0) / sig2_surv0, weights = unif (report).** At 10/25/50/75/90% of the curve the ratio was 0.967, 0.94, 0.97, 0.933, 0.931 and is now 0.987, 0.966, 1.01, 0.998, 0.997. Both numerator (R2) and denominator (R1) changed.
- **T10.** The same three cells fail as in Stage 2 and 2b, for the reasons given there: S5 checks the old linear log-CI formula, and R1 and S4 count the informative errors the Stage-2 contract requires.
- **T2.4 (duplicated rows vs weight 2).** Unchanged: all-equal weights keep the Greenwood path.

**API** (`Rscript dev/audit/stage2c/api.R compare` against `api_s2b.rds`), verbatim:

```
API compare: 46 exports before, 46 now
DIFFERENCES from api_s2b.rds:
  (none)
API: identical to api_s2b.rds
```

**CRAN check.** `R CMD check --as-cran --timings weightedsurv_0.2.0.tar.gz`. The tarball came from `R CMD build` of a clean `git clone` at `7ab4171`, the last commit touching `R/`, `man/` or `NEWS.md`. The check ran with `RSTUDIO_PANDOC` set and `DISPLAY` unset. Log: `stage2c/00check_2c.log`.

```
Status: 2 NOTEs
```

- **NOTEs.** These are the same two local-tooling NOTEs as in Stage 2b: pandoc not on PATH for README/NEWS, and old HTML Tidy.
- **Tests.** `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 613 ]`. Examples, `--run-donttest` examples and the vignette rebuild are OK.

**forestsearch 0.5.0.** `R CMD check --no-manual --ignore-vignettes` with `R_LIBS` set to a temporary library holding the HEAD build. It ran for 50 s. Logs: `stage2c/forestsearch_00check_2c.log`, `forestsearch_testthat_2c.Rout`.

```
Status: OK
[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]
```

As before, no forestsearch test or example calls weightedsurv.

## H. Gate G2 output

Checks 1–9 were run by `stage2c/g2.sh` (the block from task section 9, verbatim) after the docs commit `7ab4171`, with `START_SHA` and `MASTER_SHA` set (`stage2c/g2_output.txt`). The golden compare lines are abridged here; the full list of moved fields is in `golden_compare_R2.txt` and `g2_output.txt`.

```
OK 1a
OK 1b
API compare: 46 exports before, 46 now
DIFFERENCES from api_s2b.rds:
  (none)
API: identical to api_s2b.rds
OK 2a
OK 2b
final whitelist = union over: R1, R2
golden compare, whitelist 'final': 6569 fields in golden, 6569 in new run, 0 added
WHITELISTED fields that moved (272):
  [... 272 lines, see g2_output.txt ...]
R1 check: 71 moved analytic sig2_surv* fields checked against survfit robust; max relative SE difference 1.48e-14
Warnings: no call raises a warning it did not raise in golden_s2b.rds
GOLDEN: identical outside the whitelist
OK 3
OK 4
OK 5
OK 6a
OK 6b
OK 7a
OK 7b
OK 8
7ab4171 stage2c docs: robust variance (R1) and influence-based resampling (R2) in NEWS and help pages
5784aa1 R2: resampled weighted KM curves use per-subject influence contributions
6a325a2 R1: analytic weighted KM variance is the robust (infinitesimal-jackknife) variance
421caf7 stage2c: baseline (golden_s2b.rds, api_s2b.rds) and copied scripts
a09b5b9 stage2c: add task doc
```

Check 9 at that run printed no `FAIL: missing` line, and did not print `OK 9` because this report did not yet exist. Re-run of check 9 after this report was written:

```
OK 9 (if no FAIL above)
```

Check 10 is run after the final commit; its output is given with the reporting step.
