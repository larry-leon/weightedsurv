# weightedsurv 0.2.0: Stage 2b on branch `km-fixes-0.2.0`

**Top line.** Branch `km-fixes-0.2.0`. START_SHA `6860cb43b41cefa1de5fc710bacb27ecb4682000`, MASTER_SHA `166919ec7128f165fb380847f549fa26ca1db985`. The final SHA is the commit that adds this report; its parent is `b5b8f67`, and no commit follows it. Gates G0, G1 and G2 passed; the G2 output is in section H. All five items are DONE: T0, B9, B11, B8 and M. None is NOT DONE. M1c was run (kickoff: `M1c: GO`). Nothing was merged or pushed. Installs went only into temporary libraries under the session scratchpad.

Execution notes (findings, not package results):
- **Baseline.** `golden_s2.rds` and `api_s2.rds` were made at `d7c3900`, after the task-doc commit and before any edit to `R/`. At that point `git diff --quiet START_SHA -- R/ NAMESPACE` was true. The snapshot has 6,553 fields from 8 datasets and 164 recorded calls. None of the 164 calls raised a warning at baseline.
- **Weights in the golden set.** Every dataset has a non-constant `wv` column. Where the dataset's `w` is already non-constant (`veteran`, the three GBSG sets, `sim_weighted`), `wv` is `w`. Otherwise `wv` is seeded U(0.3, 3). The old call `df_counting_weighted` uses `w`, which is constant (= 1) in `sim_cont`, `sim_ties` and `sim_strata`. Those three calls therefore cover the constant-weights case, and they did not move.
- **Golden compare script.** The B8 commit changed `stage2b/golden.R` so that a list whose names grow by appended elements counts as an allowed addition. Before the change, the 16 `rmst_band*/@names` entries, which gained `mean_rule`, were reported as differences. Only lists that gain elements are affected; old names must be an unchanged prefix of the new ones.
- **Commit layout.** The two article edits (B8, M) and the M help-page prose went into the single Step-3 documentation commit `932ba3a`. B8's `@param mean_rule` went with the B8 code, because the new argument has to be documented for `R CMD check`.
- **Docs commit amended.** The Step-3 commit was amended once, before any later commit, to correct a NEWS sentence. The draft said that results with all weights equal do not change. They do change when the constant is not 1: such weights now give the unweighted SE, where before a constant weight c gave the SE of the data with every row repeated c times. The contract requires this ("Constant weights: all weights equal to any c > 0 give the unit-weight SE"). No golden call uses a constant weight other than 1, so the golden compare does not show it. Pressure-test T2.4 does show it (section D).
- **Which commit the checks ran on.** `R CMD check` and forestsearch ran on a clean clone at `932ba3a`. The probe, pressure test, M1c and plot sheets ran on the same `R/` code (before the amend, which changed only `NEWS.md`). Later commits touch only `dev/audit/`.

---

## A. Items

| ID | class | status | commit | files touched | golden fields that moved | tests added |
|---|---|---|---|---|---|---|
| T0 | housekeeping | DONE | `769605b` | `DESCRIPTION` | none | none (suite runs under testthat 3.3.2) |
| B9 | changes behaviour | DONE | `7834e43` | `R/km_wlr_calculations_helpers.R` (default, `@param`, example adds `show_resamples = TRUE`), `man/KM_diff.Rd`, `README.md` (the band call adds `show_resamples = TRUE`) | none | `test-b9-show-resamples.R` (2 blocks) |
| B11 | changes behaviour (drawing) | DONE | `87ca3c8` | `R/kmplotting_helpers.R`, `DESCRIPTION` (`grDevices` in Imports), `NAMESPACE` (`importFrom(grDevices,adjustcolor)`), `man/plot_km_curves_counting.Rd` | none | `test-b11-band-transparency.R` (1) |
| B8 | changes behaviour (adds option) | DONE | `fcb47ac` | `R/rmst_calculations.R`, `R/internal_helpers.R` (new `ws_cumulative_mean()`), `man/cumulative_rmst_bands.Rd` | none; 16 `rmst_band*` lists gained the element `mean_rule` | `test-b8-rmst-mean-rule.R` (4) |
| M | changes the method | DONE | `470c40b` | `R/df_counting_helpers.R` (`calculate_risk_event_counts()`, analytic path) | 165 fields, all whitelisted: `sig2_surv0` 29, `sig2_surv1` 29, `sig2_survP` 13, `sig2_dhat` 16, `lower`/`upper` 16 each, `Zdhat_star` 8, `c_alpha_band` 8, `sb_lower`/`sb_upper` 8 each, weighted `quantile_results` `lower` 6 and `upper` 8. By call: `df_counting_weighted` 20 (only the 5 datasets with non-constant `w`), `df_counting_wv` 33, `KM_diff_wv` 40, `KM_diff_wv_band` 72. All 71 moved `sig2_surv*` fields equal the approved formula (1e-9 relative). No field of an unweighted, constant-weight or `draws > 0` call moved, and neither did any RMST output | `test-m-weighted-greenwood.R` (7) |
| docs | Step 3 | done | `932ba3a` | `NEWS.md`; help-page sections on `calculate_risk_event_counts()`, `KM_diff()`, `df_counting()` (3 Rd); `vignettes/articles/weightedkm_methodology.Rmd` (the "Variance Estimation in `weightedsurv`" paragraph); `vignettes/articles/weightedrmst_methodology.Rmd` (a paragraph after the step-function RMST formula) | n/a | n/a |

Other commits: `d7c3900` (task doc), `0525f3b` (baseline and copied scripts), `b5b8f67` (verification artifacts). No existing test had to change. Suite: 13 files, 61 `test_that` blocks, 560 expectations, 0 failures (1.9 s locally; 2.4 s inside `R CMD check`). No test compares a weighted SE with `survfit`.

## B. M in numbers

`KM_diff(weight.name = "wv", draws = 0)` on each golden dataset. The SE is taken at the grid points at the 25th, 50th and 75th percentile of `at_points`. "Before" is `golden_s2.rds` and "after" is the current code. The ratio is to `survfit(weights = wv, robust = TRUE)` for the same arm. Script: `stage2b/report_numbers.R`; output: `report_numbers.txt`.

| dataset | arm | t (grid pct) | SE before | SE after | ratio before | ratio after |
|---|---|---|---|---|---|---|
| veteran | 0 | 25 (25%) | 0.05177 | 0.04646 | 1.080 | 0.969 |
| veteran | 0 | 80 (50%) | 0.06602 | 0.06326 | 1.026 | 0.983 |
| veteran | 0 | 144 (75%) | 0.06760 | 0.06470 | 1.035 | 0.990 |
| veteran | 1 | 25 (25%) | 0.06402 | 0.05599 | 1.132 | 0.990 |
| veteran | 1 | 80 (50%) | 0.07618 | 0.06846 | 1.116 | 1.003 |
| veteran | 1 | 144 (75%) | 0.06690 | 0.06214 | 1.000 | 0.929 |
| gbsg | 0 | 18.6 (25%) | 0.02363 | 0.02185 | 1.078 | 0.997 |
| gbsg | 0 | 35.5 (50%) | 0.02863 | 0.02754 | 1.046 | 1.006 |
| gbsg | 0 | 55.4 (75%) | 0.03150 | 0.03060 | 1.028 | 0.999 |
| gbsg | 1 | 18.6 (25%) | 0.02574 | 0.02474 | 1.038 | 0.997 |
| gbsg | 1 | 35.5 (50%) | 0.03312 | 0.03272 | 1.010 | 0.998 |
| gbsg | 1 | 55.4 (75%) | 0.03644 | 0.03631 | 0.985 | 0.982 |
| gbsg_er_le0 | 0 | 14.7 (25%) | 0.05702 | 0.05856 | 0.986 | 1.012 |
| gbsg_er_le0 | 0 | 25.6 (50%) | 0.06255 | 0.06427 | 0.951 | 0.977 |
| gbsg_er_le0 | 0 | 48 (75%) | 0.08129 | 0.08626 | 0.988 | 1.049 |
| gbsg_er_le0 | 1 | 14.7 (25%) | 0.09520 | 0.10131 | 0.931 | 0.991 |
| gbsg_er_le0 | 1 | 25.6 (50%) | 0.10626 | 0.11685 | 0.911 | 1.002 |
| gbsg_er_le0 | 1 | 48 (75%) | 0.09332 | 0.10910 | 0.845 | 0.988 |
| gbsg_er_gt0 | 0 | 19.6 (25%) | 0.02812 | 0.02533 | 1.121 | 1.010 |
| gbsg_er_gt0 | 0 | 35.8 (50%) | 0.03166 | 0.03016 | 1.058 | 1.008 |
| gbsg_er_gt0 | 0 | 56 (75%) | 0.03475 | 0.03374 | 1.032 | 1.002 |
| gbsg_er_gt0 | 1 | 19.6 (25%) | 0.02615 | 0.02618 | 0.999 | 1.000 |
| gbsg_er_gt0 | 1 | 35.8 (50%) | 0.03407 | 0.03425 | 0.999 | 1.004 |
| gbsg_er_gt0 | 1 | 56 (75%) | 0.03982 | 0.04018 | 0.987 | 0.996 |
| sim_cont | 0 | 2.55 (25%) | 0.03866 | 0.03921 | 0.986 | 1.000 |
| sim_cont | 0 | 5.35 (50%) | 0.05490 | 0.05593 | 0.983 | 1.001 |
| sim_cont | 0 | 11.3 (75%) | 0.05938 | 0.06131 | 0.974 | 1.006 |
| sim_cont | 1 | 2.55 (25%) | 0.04694 | 0.04777 | 0.990 | 1.007 |
| sim_cont | 1 | 5.35 (50%) | 0.05326 | 0.05428 | 0.979 | 0.998 |
| sim_cont | 1 | 11.3 (75%) | 0.06100 | 0.06253 | 0.951 | 0.974 |
| sim_ties | 0 | 3.5 (25%) | 0.09302 | 0.05049 | 1.788 | 0.971 |
| sim_ties | 0 | 6.5 (50%) | 0.09790 | 0.05822 | 1.625 | 0.967 |
| sim_ties | 0 | 13 (75%) | 0.09308 | 0.06664 | 1.498 | 1.072 |
| sim_ties | 1 | 3.5 (25%) | 0.06111 | 0.04087 | 1.502 | 1.004 |
| sim_ties | 1 | 6.5 (50%) | 0.06893 | 0.04756 | 1.422 | 0.981 |
| sim_ties | 1 | 13 (75%) | 0.08423 | 0.06231 | 1.367 | 1.011 |
| sim_weighted | 0 | 3.37 (25%) | 0.03908 | 0.03949 | 0.988 | 0.998 |
| sim_weighted | 0 | 6.72 (50%) | 0.04574 | 0.04635 | 0.982 | 0.995 |
| sim_weighted | 0 | 13.2 (75%) | 0.04527 | 0.04624 | 0.935 | 0.955 |
| sim_weighted | 1 | 3.37 (25%) | 0.02683 | 0.02709 | 0.991 | 1.001 |
| sim_weighted | 1 | 6.72 (50%) | 0.04525 | 0.04583 | 1.000 | 1.013 |
| sim_weighted | 1 | 13.2 (75%) | 0.05155 | 0.05251 | 0.992 | 1.010 |
| sim_strata | 0 | 2.59 (25%) | 0.03554 | 0.03590 | 0.978 | 0.988 |
| sim_strata | 0 | 6.36 (50%) | 0.04661 | 0.04723 | 0.992 | 1.005 |
| sim_strata | 0 | 13 (75%) | 0.05107 | 0.05232 | 0.999 | 1.024 |
| sim_strata | 1 | 2.59 (25%) | 0.03016 | 0.03044 | 0.988 | 0.997 |
| sim_strata | 1 | 6.36 (50%) | 0.04671 | 0.04733 | 1.005 | 1.019 |
| sim_strata | 1 | 13 (75%) | 0.04942 | 0.05034 | 0.979 | 0.997 |

**M1c.** `stage2b/km_m1c_montecarlo.R` is a copy of the M1b script that differs only in its output path. It ran for 17 s: 6 cells, 2,000 replicates each, n = 200. Results are in `km_m1c_results.csv`. The "package" rows equal the "weighted Greenwood" rows in every cell and at every time: mean SE, SE/SD and coverage are identical to the printed precision.

| ties | weights | t | package SE/SD | wG SE/SD | package coverage | wG coverage |
|---|---|---|---|---|---|---|
| ceiling 0.5 | none | 5 | 1.0058 | 1.0058 | 0.9445 | 0.9445 |
| ceiling 0.5 | none | 10 | 0.9776 | 0.9776 | 0.9450 | 0.9450 |
| ceiling 0.5 | none | 15 | 1.0008 | 1.0008 | 0.9465 | 0.9465 |
| ceiling 0.5 | U(0.3,3) | 5 | 0.9969 | 0.9969 | 0.9405 | 0.9405 |
| ceiling 0.5 | U(0.3,3) | 10 | 0.9779 | 0.9779 | 0.9440 | 0.9440 |
| ceiling 0.5 | U(0.3,3) | 15 | 0.9956 | 0.9956 | 0.9425 | 0.9425 |
| ceiling 0.5 | U(0.3,3)/sum | 5 | 0.9969 | 0.9969 | 0.9405 | 0.9405 |
| ceiling 0.5 | U(0.3,3)/sum | 10 | 0.9779 | 0.9779 | 0.9440 | 0.9440 |
| ceiling 0.5 | U(0.3,3)/sum | 15 | 0.9956 | 0.9956 | 0.9425 | 0.9425 |
| continuous | none | 5 | 1.0060 | 1.0060 | 0.9475 | 0.9475 |
| continuous | none | 10 | 0.9777 | 0.9777 | 0.9415 | 0.9415 |
| continuous | none | 15 | 1.0008 | 1.0008 | 0.9425 | 0.9425 |
| continuous | U(0.3,3) | 5 | 0.9970 | 0.9970 | 0.9450 | 0.9450 |
| continuous | U(0.3,3) | 10 | 0.9781 | 0.9781 | 0.9455 | 0.9455 |
| continuous | U(0.3,3) | 15 | 0.9958 | 0.9958 | 0.9435 | 0.9435 |
| continuous | U(0.3,3)/sum | 5 | 0.9970 | 0.9970 | 0.9450 | 0.9450 |
| continuous | U(0.3,3)/sum | 10 | 0.9781 | 0.9781 | 0.9455 | 0.9455 |
| continuous | U(0.3,3)/sum | 15 | 0.9958 | 0.9958 | 0.9435 | 0.9435 |

Reading:
- **Unweighted cells.** The package rows are unchanged from M1b. In the weighted cells, Stage 2 had SE/SD from 0.45 to 4.32 and coverage from 0.63 to 1.00 (`stage2/km_m1b_results.csv`). It now has SE/SD from 0.978 to 0.997 and coverage from 0.941 to 0.946, and the U(0.3,3) and U(0.3,3)/sum cells are identical.
- **Golden datasets.** The ratio of the weighted SE to `survfit` robust was 0.85–1.79 before and is 0.93–1.07 now. The largest shift is in `sim_ties`, where the old weighted formula was furthest off.
- **What does not move.** The resampling path (`draws > 0` with unequal weights) is byte-identical to the old code (test `m-weighted-greenwood`, last block, `expect_identical`). The golden `KM_diff_weighted_draws` and `cox_*_d200` calls did not move.

## C. B8

`KM_diff(draws.band = 1000)` on the default grid, then `cumulative_rmst_bands()` with defaults (`draws_sb = 1000`) under each rule. "survfit" is `summary(survfit(Surv ~ treat), rmean = tau)` treated minus control, with tau the last grid time.

| dataset | RMST trapezoid | RMST rectangle | survfit rmean diff | rectangle − survfit | c_alpha trapezoid | c_alpha rectangle |
|---|---|---|---|---|---|---|
| veteran (tau = 999) | 15.139392 | 18.133115 | 18.133115 | 6e-14 | 2.6306 | 2.6732 |
| gbsg (tau = 87.36) | 8.468821 | 8.434356 | 8.434356 | 6.4e-14 | 2.7253 | 2.7287 |
| sim_cont (tau = 28.39) | 2.578505 | 2.594907 | 2.594907 | -4.4e-15 | 2.6748 | 2.6756 |

Reading:
- **The rectangle rule.** It reproduces `survfit`'s restricted-mean difference to rounding. The tests check this at the last grid point and three interior ones, on veteran, GBSG and a simulated set with ties, both on the default `KM_diff()` grid and on `fit_itt` from `plotKM.band_subgroups()`. The trapezoid default is unchanged. On veteran it differs from the exact area by 3.0 days of 999, because that grid has long gaps between events.
- **Band critical values.** They differ slightly between rules (2.63 vs 2.67 on veteran) because the same resampled curves are integrated differently. The draws, seeds, the 0.95 quantile and the 1.96 are the same under both rules (test `b8`, last block).

## D. Seed probe and pressure test, Stage 2 → Stage 2b

**Probe** (`dev/audit/km_audit_probe.R`; output in `stage2/probe_after.txt` → `stage2b/probe_2b.txt`):

| ID | Stage 2 | Stage 2b |
|---|---|---|
| OK1, OK2, OK3 | PASS | PASS (OK3: pkg/survfit robust SE in [0.943, 1.022], was [0.923, 1.005]) |
| K3 (SE not scale-free) | REPRODUCED 0.18/0.64/1.60/6.31 | NOT_REPRODUCED 1.00/1.00/1.00/1.00 (M) |
| K4 (tied times vs survfit robust) | REPRODUCED 2.06/1.77/1.52/1.25 | NOT_REPRODUCED 1.00/1.02/0.97/1.01 (M) |
| K5 (weights 2 vs nearly 2) | REPRODUCED 0.0297 vs 0.0416 | NOT_REPRODUCED 0.0420 / 0.0420 / 0.0420 (M) |
| I3 (trapezoid vs exact area) | REPRODUCED 4.2836 vs 4.2688 | REPRODUCED, unchanged: the probe uses the default rule, which stays the trapezoid (B8 adds `"rectangle"`) |
| all other items | NOT_REPRODUCED | NOT_REPRODUCED, unchanged |

**Pressure test** (`stage2b/km_pressure_test_s2b.R`, a copy that changes only the default output directory; 96 s). Failing cells by family:

| family | cells | failing Stage 2 | failing Stage 2b |
|---|---|---|---|
| REF | 1 | 0 | 0 |
| T1 | 111 | 0 | 0 |
| T2 | 17 | 3 | 0 |
| T3 | 288 | 0 | 0 |
| T4 | 21 | 0 | 0 |
| T5 | 9 | 0 | 0 |
| T6 | 24 | 0 | 0 |
| T7 | 18 | 0 | 0 |
| T8 | 73 | 0 | 0 |
| T9 | 1 | 0 | 0 |
| T10 | 5 | 3 | 3 |

Changes, each explained:
- **T2.2, weights w → c·w** (c = 10, 0.01, 1/sum w). As expected, the three "SE ratio (K3)" cells now pass: the worst ratio is 1, where Stage 2 had 1.41, 50.5 and 168.
- **T2.4, duplicated rows vs weight 2.** The SE part is the expected change by design. The script records it as a report-only note, not as a failing check, so the cell still counts as passing on its survival check. That note was "SE difference max 0" in Stage 1 and Stage 2 and is now "max 0.0181". The old code put all-equal weights through Greenwood with the weighted counts, so weight 2 behaved like duplicated rows. Weight 2 now gives the unweighted SE, so weights are no longer frequency weights.
- **T2.8, SE with draws = 2000 / draws = 0, weighted (report).** Quartiles were 0.361, 0.875 and 1.01, and are now 0.921, 0.978 and 0.987. The resampling SE is unchanged; the analytic SE it is divided by changed, and the two now agree more closely.
- **T5, var(resampled surv0) / sig2_surv0, weights = unif (report).** The ratios at 10/25/50/75/90% of the curve were 0.967, 0.94, 0.97, 0.933 and 0.931, and are now 0.994, 0.968, 1.01, 1.01 and 1.05. The resampled curves are unchanged; the denominator is the new weighted SE. The unweighted row is unchanged.
- **T10.** The same three cells fail as in Stage 2, for the same reasons. In S5, the test computes the old linear log-CI formula itself rather than reading what is drawn. R1 and S4 count the informative errors the Stage-2 contract requires as failures.
- No other cell changed `n_checks` or `n_fail`.

**T1: weighted SE / `survfit(robust = TRUE)`.** Each entry is the median and range of the per-arm medians over the T1 design grid (n 2 to 500, censoring 0 to 0.95, last observation as-is / event / censored).

| weights | ties | median before | range before | median after | range after |
|---|---|---|---|---|---|
| int | none | 0.964 | 0.458 to 1.13 | 0.999 | 0.683 to 1.7 |
| int | moderate | 0.975 | 0.354 to 2.67 | 1 | 0.745 to 1.68 |
| int | heavy | 1.42 | 0.358 to 8.5 | 1 | 0.75 to 1.58 |
| ipw | none | 0.898 | 0.222 to 2.53 | 0.98 | 0.226 to 5.69 |
| ipw | moderate | 0.958 | 0.398 to 2.43 | 0.986 | 0.352 to 4.34 |
| ipw | heavy | 1.39 | 0.525 to 6.59 | 0.984 | 0.397 to 4.83 |
| norm | none | 1.84 | 0.205 to 1.01e+03 | 0.999 | 0.749 to 2.63 |
| norm | moderate | 2.21 | 0.313 to 975 | 0.999 | 0.735 to 2.82 |
| norm | heavy | 2.47 | 0.282 to 600 | 0.999 | 0.71 to 2.41 |
| unif | none | 0.941 | 0.28 to 1.15 | 0.999 | 0.715 to 2.98 |
| unif | moderate | 0.934 | 0.00542 to 2.84 | 0.999 | 0.726 to 2.41 |
| unif | heavy | 1.41 | 0.16 to 8.24 | 0.999 | 0.762 to 1.82 |
| zero10 | none | 0.933 | 0.233 to 1.19 | 1 | 0.685 to 2.76 |
| zero10 | moderate | 0.968 | 0.00101 to 2.68 | 0.999 | 0.699 to 2.67 |
| zero10 | heavy | 1.46 | 0.187 to 7.86 | 1 | 0.724 to 2.45 |

- **Medians.** Every median is now 0.98–1.00. The scale-dependent extremes are gone; they came from normalised weights (`norm`) and from small n with near-zero ratios.
- **Remaining spread.** It reflects what weighted Greenwood and the robust (sandwich) estimator are, not scale. It is widest for IPW weights (to 5.69), and the worst cells are at n = 100–500 with no censoring.

**Plot sheets.** `stage2b/km_plot_compare.R after2` writes `km_plots_after2.pdf` and `km_plot_compare_after2.csv`. The CSV is `identical()` to `stage2/km_plot_compare_after.csv` (74 rows, 13 columns), so B11 left the geometry unchanged. On the sheets, the two bands are drawn at alpha 0.5 and the overlap shows as a darker region.

## E. API

`Rscript dev/audit/stage2b/api.R compare` against `api_s2.rds`, printed verbatim:

```
API compare: 46 exports before, 46 now
DIFFERENCES from api_s2.rds:
  cumulative_rmst_bands: appended mean_rule = c("trapezoid", "rectangle")
  KM_diff: default of 'show_resamples' changed: TRUE -> FALSE
API: only the two allowed changes
```

## F. CRAN check and forestsearch

`R CMD check --as-cran --timings` on `weightedsurv_0.2.0.tar.gz`, built from a clean clone at `932ba3a` (`stage2b/00check_2b.log`):

```
Status: 2 NOTEs
```

The NOTEs, verbatim:

```
* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
```

- **Notes.** Both come from local tooling. Stage 2's third NOTE ("unable to verify current time") did not appear this time.
- **Tests.** `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 560 ]`, elapsed 2.4 s. Examples, `--run-donttest` examples and the vignette rebuild (11 s) are OK.

forestsearch 0.5.0: `R CMD check --no-manual --ignore-vignettes` against a temporary library holding the new weightedsurv build (`stage2b/forestsearch_00check_2b.log`, `forestsearch_testthat_2b.Rout`):

```
Status: OK
[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]
```

The check took 45 s. As in Stage 2, no forestsearch test or example calls weightedsurv.

## G. Sentences in `vignettes/` or `README.md` made inaccurate and not edited

None found. Searched every file in `vignettes/`, `vignettes/articles/` and `README.md` for variance and SE wording, survfit/adjustedCurves comparisons, frequency/case weights, the integration rule, `show_resamples` and band colours. Checked and not listed:
- **B9.** Every vignette call that draws resamples sets `show_resamples` explicitly (`weightedsurv_examples.Rmd:356, 379, 403, 531, 552, 613`; `martingale_resampling_applications.Rmd:294, 337, 380`). These are `plotKM.band_subgroups()` calls, whose own default was already `FALSE`. No vignette calls `KM_diff()` with `draws.band > 0`. The one README call that relied on the old default was edited in B9.
- **Weighted-SE comparison.** `weightedkm_methodology.Rmd:326` says the examples vignette includes "cross-validation of weighted SEs against the `adjustedCurves` and `survfit` packages". The vignette still makes those comparisons (`weightedsurv_examples.Rmd:698–768`). I re-ran its Rotterdam SE comparison at START_SHA and now:
  - pkg/adjustedCurves SE ratio, control arm: median 1.227 → 1.028.
  - pkg/adjustedCurves SE ratio, treated arm: median 0.840 → 0.834.
  - pkg/`survfit` robust: medians 1.203 → 1.008 (control) and 0.959 → 0.958 (treated).

  The agreement was not exact before and is not exact now, so M does not make the sentence false.
- **Exact-area wording.** `weightedrmst_methodology.Rmd:80–88` and `:211` define the RMST as the exact area under the KM curve. The default trapezoid rule already differed from that before this task (probe I3). The edited paragraph at `:90–96` now says so.
- **Resampling recommendation.** `weightedkm_methodology.Rmd:209` says that "variance estimation via martingale resampling is preferred" for complex weighting. That is a recommendation, not a statement M falsifies.

## H. Gate G2 output

Run after the last item commit, with this report in place (checks 1–8). Check 9 was run after the commit that adds this report.

```
OK 1a
OK 1b
OK 2a
OK 2b
API compare: 46 exports before, 46 now
DIFFERENCES from api_s2.rds:
  cumulative_rmst_bands: appended mean_rule = c("trapezoid", "rectangle")
  KM_diff: default of 'show_resamples' changed: TRUE -> FALSE
API: only the two allowed changes
OK 2c
final whitelist = union over: T0, B9, B11, B8, M
golden compare, whitelist 'final': 6553 fields in golden, 6569 in new run, 16 added
ADDED (allowed):
  veteran/rmst_band/mean_rule
  veteran/rmst_band_wv/mean_rule
  gbsg/rmst_band/mean_rule
  gbsg/rmst_band_wv/mean_rule
  gbsg_er_le0/rmst_band/mean_rule
  gbsg_er_le0/rmst_band_wv/mean_rule
  gbsg_er_gt0/rmst_band/mean_rule
  gbsg_er_gt0/rmst_band_wv/mean_rule
  sim_cont/rmst_band/mean_rule
  sim_cont/rmst_band_wv/mean_rule
  sim_ties/rmst_band/mean_rule
  sim_ties/rmst_band_wv/mean_rule
  sim_weighted/rmst_band/mean_rule
  sim_weighted/rmst_band_wv/mean_rule
  sim_strata/rmst_band/mean_rule
  sim_strata/rmst_band_wv/mean_rule
LISTS THAT GAINED ELEMENTS (allowed, 16):
  veteran/rmst_band@names: + mean_rule
  veteran/rmst_band_wv@names: + mean_rule
  gbsg/rmst_band@names: + mean_rule
  gbsg/rmst_band_wv@names: + mean_rule
  gbsg_er_le0/rmst_band@names: + mean_rule
  gbsg_er_le0/rmst_band_wv@names: + mean_rule
  gbsg_er_gt0/rmst_band@names: + mean_rule
  gbsg_er_gt0/rmst_band_wv@names: + mean_rule
  sim_cont/rmst_band@names: + mean_rule
  sim_cont/rmst_band_wv@names: + mean_rule
  sim_ties/rmst_band@names: + mean_rule
  sim_ties/rmst_band_wv@names: + mean_rule
  sim_weighted/rmst_band@names: + mean_rule
  sim_weighted/rmst_band_wv@names: + mean_rule
  sim_strata/rmst_band@names: + mean_rule
  sim_strata/rmst_band_wv@names: + mean_rule
WHITELISTED fields that moved (165):
  veteran/df_counting_weighted/quantile_results/upper: max abs diff 7
  veteran/df_counting_weighted/sig2_surv0: max abs diff 0.000878
  veteran/df_counting_weighted/sig2_surv1: max abs diff 0.00122
  veteran/df_counting_weighted/sig2_survP: max abs diff 0.00151
  veteran/df_counting_wv/quantile_results/upper: max abs diff 7
  veteran/df_counting_wv/sig2_surv0: max abs diff 0.000878
  veteran/df_counting_wv/sig2_surv1: max abs diff 0.00122
  veteran/df_counting_wv/sig2_survP: max abs diff 0.00151
  veteran/KM_diff_wv/sig2_surv0: max abs diff 0.000878
  veteran/KM_diff_wv/sig2_surv1: max abs diff 0.00122
  veteran/KM_diff_wv/sig2_dhat: max abs diff 0.00169
  veteran/KM_diff_wv/lower: max abs diff 0.0188
  veteran/KM_diff_wv/upper: max abs diff 0.0188
  veteran/KM_diff_wv_band/sig2_surv0: max abs diff 0.000878
  veteran/KM_diff_wv_band/sig2_surv1: max abs diff 0.00122
  veteran/KM_diff_wv_band/sig2_dhat: max abs diff 0.00169
  veteran/KM_diff_wv_band/lower: max abs diff 0.0188
  veteran/KM_diff_wv_band/upper: max abs diff 0.0188
  veteran/KM_diff_wv_band/Zdhat_star: max abs diff 0.925
  veteran/KM_diff_wv_band/c_alpha_band: max abs diff 0.0504
  veteran/KM_diff_wv_band/sb_lower: max abs diff 0.0273
  veteran/KM_diff_wv_band/sb_upper: max abs diff 0.0273
  gbsg/df_counting_weighted/quantile_results/lower: max abs diff 0.0329
  gbsg/df_counting_weighted/quantile_results/upper: max abs diff 5.78
  gbsg/df_counting_weighted/sig2_surv0: max abs diff 0.00544
  gbsg/df_counting_weighted/sig2_surv1: max abs diff 0.00035
  gbsg/df_counting_weighted/sig2_survP: max abs diff 0.000703
  gbsg/df_counting_wv/quantile_results/lower: max abs diff 0.0329
  gbsg/df_counting_wv/quantile_results/upper: max abs diff 5.78
  gbsg/df_counting_wv/sig2_surv0: max abs diff 0.00544
  gbsg/df_counting_wv/sig2_surv1: max abs diff 0.00035
  gbsg/df_counting_wv/sig2_survP: max abs diff 0.000703
  gbsg/KM_diff_wv/sig2_surv0: max abs diff 0.00544
  gbsg/KM_diff_wv/sig2_surv1: max abs diff 0.00035
  gbsg/KM_diff_wv/sig2_dhat: max abs diff 0.00579
  gbsg/KM_diff_wv/lower: max abs diff 0.0532
  gbsg/KM_diff_wv/upper: max abs diff 0.0532
  gbsg/KM_diff_wv_band/sig2_surv0: max abs diff 0.00544
  gbsg/KM_diff_wv_band/sig2_surv1: max abs diff 0.00035
  gbsg/KM_diff_wv_band/sig2_dhat: max abs diff 0.00579
  gbsg/KM_diff_wv_band/lower: max abs diff 0.0532
  gbsg/KM_diff_wv_band/upper: max abs diff 0.0532
  gbsg/KM_diff_wv_band/Zdhat_star: max abs diff 0.804
  gbsg/KM_diff_wv_band/c_alpha_band: max abs diff 0.0266
  gbsg/KM_diff_wv_band/sb_lower: max abs diff 0.0818
  gbsg/KM_diff_wv_band/sb_upper: max abs diff 0.0818
  gbsg_er_le0/df_counting_weighted/sig2_surv0: max abs diff 0.00109
  gbsg_er_le0/df_counting_weighted/sig2_surv1: max abs diff 0.00319
  gbsg_er_le0/df_counting_weighted/sig2_survP: max abs diff 0.000493
  gbsg_er_le0/df_counting_wv/sig2_surv0: max abs diff 0.00109
  gbsg_er_le0/df_counting_wv/sig2_surv1: max abs diff 0.00319
  gbsg_er_le0/df_counting_wv/sig2_survP: max abs diff 0.000493
  gbsg_er_le0/KM_diff_wv/sig2_surv0: max abs diff 0.00109
  gbsg_er_le0/KM_diff_wv/sig2_surv1: max abs diff 0.00319
  gbsg_er_le0/KM_diff_wv/sig2_dhat: max abs diff 0.00428
  gbsg_er_le0/KM_diff_wv/lower: max abs diff 0.0315
  gbsg_er_le0/KM_diff_wv/upper: max abs diff 0.0315
  gbsg_er_le0/KM_diff_wv_band/sig2_surv0: max abs diff 0.00109
  gbsg_er_le0/KM_diff_wv_band/sig2_surv1: max abs diff 0.00319
  gbsg_er_le0/KM_diff_wv_band/sig2_dhat: max abs diff 0.00428
  gbsg_er_le0/KM_diff_wv_band/lower: max abs diff 0.0315
  gbsg_er_le0/KM_diff_wv_band/upper: max abs diff 0.0315
  gbsg_er_le0/KM_diff_wv_band/Zdhat_star: max abs diff 0.368
  gbsg_er_le0/KM_diff_wv_band/c_alpha_band: max abs diff 0.151
  gbsg_er_le0/KM_diff_wv_band/sb_lower: max abs diff 0.0221
  gbsg_er_le0/KM_diff_wv_band/sb_upper: max abs diff 0.0221
  gbsg_er_gt0/df_counting_weighted/quantile_results/lower: max abs diff 0.0657
  gbsg_er_gt0/df_counting_weighted/sig2_surv0: max abs diff 0.0149
  gbsg_er_gt0/df_counting_weighted/sig2_surv1: max abs diff 0.0002
  gbsg_er_gt0/df_counting_weighted/sig2_survP: max abs diff 0.001
  gbsg_er_gt0/df_counting_wv/quantile_results/lower: max abs diff 0.0657
  gbsg_er_gt0/df_counting_wv/sig2_surv0: max abs diff 0.0149
  gbsg_er_gt0/df_counting_wv/sig2_surv1: max abs diff 0.0002
  gbsg_er_gt0/df_counting_wv/sig2_survP: max abs diff 0.001
  gbsg_er_gt0/KM_diff_wv/sig2_surv0: max abs diff 0.0149
  gbsg_er_gt0/KM_diff_wv/sig2_surv1: max abs diff 0.0002
  gbsg_er_gt0/KM_diff_wv/sig2_dhat: max abs diff 0.0151
  gbsg_er_gt0/KM_diff_wv/lower: max abs diff 0.116
  gbsg_er_gt0/KM_diff_wv/upper: max abs diff 0.116
  gbsg_er_gt0/KM_diff_wv_band/sig2_surv0: max abs diff 0.0149
  gbsg_er_gt0/KM_diff_wv_band/sig2_surv1: max abs diff 0.0002
  gbsg_er_gt0/KM_diff_wv_band/sig2_dhat: max abs diff 0.0151
  gbsg_er_gt0/KM_diff_wv_band/lower: max abs diff 0.116
  gbsg_er_gt0/KM_diff_wv_band/upper: max abs diff 0.116
  gbsg_er_gt0/KM_diff_wv_band/Zdhat_star: max abs diff 1.29
  gbsg_er_gt0/KM_diff_wv_band/c_alpha_band: max abs diff 0.0158
  gbsg_er_gt0/KM_diff_wv_band/sb_lower: max abs diff 0.182
  gbsg_er_gt0/KM_diff_wv_band/sb_upper: max abs diff 0.182
  sim_cont/df_counting_wv/quantile_results/upper: max abs diff 2.96
  sim_cont/df_counting_wv/sig2_surv0: max abs diff 0.000894
  sim_cont/df_counting_wv/sig2_surv1: max abs diff 0.00276
  sim_cont/df_counting_wv/sig2_survP: max abs diff 0.000472
  sim_cont/KM_diff_wv/sig2_surv0: max abs diff 0.000894
  sim_cont/KM_diff_wv/sig2_surv1: max abs diff 0.00276
  sim_cont/KM_diff_wv/sig2_dhat: max abs diff 0.00348
  sim_cont/KM_diff_wv/lower: max abs diff 0.04
  sim_cont/KM_diff_wv/upper: max abs diff 0.04
  sim_cont/KM_diff_wv_band/sig2_surv0: max abs diff 0.000894
  sim_cont/KM_diff_wv_band/sig2_surv1: max abs diff 0.00276
  sim_cont/KM_diff_wv_band/sig2_dhat: max abs diff 0.00348
  sim_cont/KM_diff_wv_band/lower: max abs diff 0.04
  sim_cont/KM_diff_wv_band/upper: max abs diff 0.04
  sim_cont/KM_diff_wv_band/Zdhat_star: max abs diff 0.6
  sim_cont/KM_diff_wv_band/c_alpha_band: max abs diff 0.163
  sim_cont/KM_diff_wv_band/sb_lower: max abs diff 0.0463
  sim_cont/KM_diff_wv_band/sb_upper: max abs diff 0.0463
  sim_ties/df_counting_wv/quantile_results/lower: max abs diff 2
  sim_ties/df_counting_wv/quantile_results/upper: max abs diff 1.5
  sim_ties/df_counting_wv/sig2_surv0: max abs diff 0.00866
  sim_ties/df_counting_wv/sig2_surv1: max abs diff 0.00503
  sim_ties/df_counting_wv/sig2_survP: max abs diff 0.00644
  sim_ties/KM_diff_wv/sig2_surv0: max abs diff 0.00866
  sim_ties/KM_diff_wv/sig2_surv1: max abs diff 0.00503
  sim_ties/KM_diff_wv/sig2_dhat: max abs diff 0.0135
  sim_ties/KM_diff_wv/lower: max abs diff 0.114
  sim_ties/KM_diff_wv/upper: max abs diff 0.114
  sim_ties/KM_diff_wv_band/sig2_surv0: max abs diff 0.00866
  sim_ties/KM_diff_wv_band/sig2_surv1: max abs diff 0.00503
  sim_ties/KM_diff_wv_band/sig2_dhat: max abs diff 0.0135
  sim_ties/KM_diff_wv_band/lower: max abs diff 0.114
  sim_ties/KM_diff_wv_band/upper: max abs diff 0.114
  sim_ties/KM_diff_wv_band/Zdhat_star: max abs diff 1.58
  sim_ties/KM_diff_wv_band/c_alpha_band: max abs diff 0.511
  sim_ties/KM_diff_wv_band/sb_lower: max abs diff 0.083
  sim_ties/KM_diff_wv_band/sb_upper: max abs diff 0.083
  sim_weighted/df_counting_weighted/quantile_results/upper: max abs diff 0.0139
  sim_weighted/df_counting_weighted/sig2_surv0: max abs diff 0.000387
  sim_weighted/df_counting_weighted/sig2_surv1: max abs diff 0.000756
  sim_weighted/df_counting_weighted/sig2_survP: max abs diff 0.000153
  sim_weighted/df_counting_wv/quantile_results/upper: max abs diff 0.0139
  sim_weighted/df_counting_wv/sig2_surv0: max abs diff 0.000387
  sim_weighted/df_counting_wv/sig2_surv1: max abs diff 0.000756
  sim_weighted/df_counting_wv/sig2_survP: max abs diff 0.000153
  sim_weighted/KM_diff_wv/sig2_surv0: max abs diff 0.000387
  sim_weighted/KM_diff_wv/sig2_surv1: max abs diff 0.000756
  sim_weighted/KM_diff_wv/sig2_dhat: max abs diff 0.00114
  sim_weighted/KM_diff_wv/lower: max abs diff 0.0165
  sim_weighted/KM_diff_wv/upper: max abs diff 0.0165
  sim_weighted/KM_diff_wv_band/sig2_surv0: max abs diff 0.000387
  sim_weighted/KM_diff_wv_band/sig2_surv1: max abs diff 0.000756
  sim_weighted/KM_diff_wv_band/sig2_dhat: max abs diff 0.00114
  sim_weighted/KM_diff_wv_band/lower: max abs diff 0.0165
  sim_weighted/KM_diff_wv_band/upper: max abs diff 0.0165
  sim_weighted/KM_diff_wv_band/Zdhat_star: max abs diff 0.366
  sim_weighted/KM_diff_wv_band/c_alpha_band: max abs diff 0.09
  sim_weighted/KM_diff_wv_band/sb_lower: max abs diff 0.0199
  sim_weighted/KM_diff_wv_band/sb_upper: max abs diff 0.0199
  sim_strata/df_counting_wv/quantile_results/lower: max abs diff 0.0558
  sim_strata/df_counting_wv/sig2_surv0: max abs diff 0.000647
  sim_strata/df_counting_wv/sig2_surv1: max abs diff 0.00419
  sim_strata/df_counting_wv/sig2_survP: max abs diff 0.000413
  sim_strata/KM_diff_wv/sig2_surv0: max abs diff 0.000647
  sim_strata/KM_diff_wv/sig2_surv1: max abs diff 0.00419
  sim_strata/KM_diff_wv/sig2_dhat: max abs diff 0.00483
  sim_strata/KM_diff_wv/lower: max abs diff 0.0609
  sim_strata/KM_diff_wv/upper: max abs diff 0.0609
  sim_strata/KM_diff_wv_band/sig2_surv0: max abs diff 0.000647
  sim_strata/KM_diff_wv_band/sig2_surv1: max abs diff 0.00419
  sim_strata/KM_diff_wv_band/sig2_dhat: max abs diff 0.00483
  sim_strata/KM_diff_wv_band/lower: max abs diff 0.0609
  sim_strata/KM_diff_wv_band/upper: max abs diff 0.0609
  sim_strata/KM_diff_wv_band/Zdhat_star: max abs diff 1.34
  sim_strata/KM_diff_wv_band/c_alpha_band: max abs diff 0.161
  sim_strata/KM_diff_wv_band/sb_lower: max abs diff 0.0829
  sim_strata/KM_diff_wv_band/sb_upper: max abs diff 0.0829
M formula check: 71 moved sig2_surv* fields checked against the weighted Greenwood formula
Warnings: no call raises a warning it did not raise in golden_s2.rds
GOLDEN: identical outside the whitelist
OK 3
OK 4
OK 5a
OK 5b
OK 6a
OK 6b
OK 6c
OK 7
OK 8 (if no FAIL above)
b5b8f67 stage2b: verification artifacts (pressure test, M1c, plot sheets after2, R CMD check, forestsearch, report numbers)
932ba3a stage2b docs: NEWS for M, B8, B9, B11; help pages for M; article passages for B8 and M
470c40b M: weighted Kaplan-Meier variance is weighted Greenwood
fcb47ac B8: cumulative_rmst_bands(mean_rule = c("trapezoid", "rectangle"))
87ca3c8 B11: semi-transparent KM confidence bands (alpha 0.5, same hues)
7834e43 B9: KM_diff(show_resamples = FALSE) by default
769605b T0: testthat (>= 3.1.7) in Suggests
0525f3b stage2b: baseline (golden_s2.rds, api_s2.rds) and copied scripts
d7c3900 stage2b: add task doc
```
