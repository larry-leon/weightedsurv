# weightedsurv 0.2.0: Stage 2d on branch `km-fixes-0.2.0`

**Top line.**
- **Branch and SHAs.** Branch `km-fixes-0.2.0`. START_SHA is `b5843ecf365d123c26f67c1bf625a012988c307d` and MASTER_SHA is `166919ec7128f165fb380847f549fa26ca1db985`. The final SHA is the commit that adds this report together with the verification artifacts; its parent is `8f820bb` (docs), and no commit follows it.
- **Gates.** G0, G1 and G2 passed. The G2 output is in section G.
- **Items.** L1, L2, L3, L4 and B12 are DONE; none is NOT DONE.
- **M6.** Ran (kickoff: `M6: GO`): 11 cells × 2,000 replicates on 12 cores in 2.9 min (task estimate 3–8 min). Both required checks passed.
- **Constraints kept.** This task ran no merge and no push. Nothing under `vignettes/` was edited. Installs went only into temporary libraries under the session scratchpad. No package dependency was added. No error or warning message was added.

Execution notes (findings, not package results):
- **Baseline.** `stage2d/golden_s2c.rds` (16,829 fields, 8 datasets, 504 recorded calls, none of which raised a warning) and `stage2d/api_s2c.rds` (46 exports) were made at `868bda6`, the task-doc commit. At that point `git diff --quiet START_SHA -- R/ NAMESPACE` was true. The "before" worktree of START_SHA lived in the scratchpad and was used by Step V.
- **Golden script.** `stage2d/golden.R` is a copy of `stage2c/golden.R` with the same datasets and the same Stage-2c calls. It changes four things:
  - Calls added, each with weights `wv`, with no weights and with a constant weight column `wc` (= 2): `df_counting()` FH(0,0) and FH(0,1), and both with `strata.name` on the four datasets with strata; `wlr_cumulative(<df_counting result>, return_cumulative = TRUE)` for both schemes; `wlr_cumulative()` on a bare data frame of counts (cumulative and not); `wlr_dhat_estimates(tzero = median of the grid)` for both schemes; `cox_rhogamma(draws = 0)` for FH(0,0), FH(0,1), MB and schemper, and `draws = 200` for FH(0,0); `plotKM.band_subgroups(draws.band = 200)` without `draws`.
  - Whitelists: L1, L2, L3, L4, B12.
  - "Identical" is exact: a numeric field outside the whitelist fails on any difference at all (Stage 2c allowed 1e-12).
  - Contract check: every moved FH(0,0) `z.score` (unstratified, stratified, and `cox_rhogamma()`) and `sig_bhat_asy` is compared with `coxph(weights = , robust = TRUE, ties = "breslow")`.
- **One whitelist entry is placed differently from the task document.** `fit$wald_debiased1`, `z.score_debiased` and `fit_resamples$Z_bstar` are listed under L3 in the task document. Each is divided by an L2 quantity (`sig_bhat_asy` or `sig2_score`), so they already moved at the L2 commit. `golden.R` whitelists them from L2 on. The final whitelist (the union) is the one in the task document.
- **API script.** `stage2d/api.R compare` exits 0 only when the B12 default is the single difference; `compare pre` (used before B12) accepts none.
- **Scope switch.** The switch is `length(unique(w)) > 1` on the subject weights of the call, as the task document writes it and as Stage 2c did. In `df_counting()` it is evaluated once on the whole sample and then applied in every stratum. Zero weights count as a distinct value: 0/1 weights therefore take the new variance. On a 300-row example with 30 zero weights, z is −0.7569 with 0/1 weights and −0.7554 when the same 270 rows are analysed without weights (the score is the same, −4.677; `coxph` robust on the 270 rows gives 0.7569 in absolute value).
- **How `df_counting()` hands subject-level data to `wlr_cumulative()`.** As an attribute (`ws_subject`) on the data frame of counts, set only when the weights are not all equal. No signature changed.
- **Implementation.** One internal helper, `ws_score_psi()` in `R/internal_helpers.R`, computes the per-subject contributions for all four L items; two small companions find the subject-level data and compute the cumulative V(t). All are `@noRd`.
- **M6 commit.** M6, the pressure test and the vignette "after" build ran on the `R/` code of `53b5748` / `8f820bb`. The diff between those two commits under `R/` is roxygen comments only (checked).
- **Remote branch.** The local reflog of `origin/km-fixes-0.2.0` records "update by push" to `8f820bb` at 19:27:21 on 2026-11-19, while this task was running, and `git ls-remote` confirms the remote branch is at `8f820bb`. No command of this task pushed; GitHub Desktop was open on the machine. The final commit is local only. `master` is at MASTER_SHA locally and on the remote.
- **Constraints kept.** This task ran no merge and no push. Nothing under `vignettes/` was edited. at `R CMD build` because `RSTUDIO_PANDOC` pointed at the directory above the binary. The check that counts ran with `RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64` and `DISPLAY` unset.

---

## A. Items

| ID | status | commit | files touched | golden fields that moved | tests added or replaced |
|---|---|---|---|---|---|
| L1 | DONE | `33b3dfd` | `R/internal_helpers.R` (new `ws_score_psi()`, `ws_lr_subject_data()`, `ws_lr_cumvar()`); `R/km_wlr_calculations_helpers.R` (`wlr_cumulative()`); `R/df_counting.R` (attaches the subject-level data, unstratified and per stratum) | 187, all whitelisted: `z.score`, `sig2_lr`, `get_score/sig2.score`, `get_score/z.score`, `zlogrank_text` (where the printed p changes) in `df_counting_wv`, `df_counting_fh01_wv`, `df_counting_weighted` (5 datasets); the same plus `sig2_lr_stratified`, `z.score_stratified` in the two stratified `_wv` calls (4 datasets); `sig2.score`, `z.score` of `wlr_cum_fh00_wv` and `wlr_cum_fh01_wv`. All 17 moved FH(0,0) z-scores (8 `wv`, 5 `w`, 4 stratified) equal `coxph` robust; largest relative difference of z² 5.4e-14. `stage2d/golden_compare_L1.txt` | **Added:** `test-l1-weighted-logrank.R` (6 blocks) and `helper-score-resid.R` (the reference implementation of the task document, as given). **Replaced:** none; no existing test encoded the earlier weighted log-rank variance |
| L2 | DONE | `20a3b5f` | `R/weightedCox.R` (`cox_score_rhogamma(score_only = FALSE)`) | 385 more (572 in total): `fit/sig2_score`, `fit/sig_bhat_asy`, `fit/z.score`, `z.score`, `hr_ci_asy/sig_bhat`, `lower`, `upper`, `zlogrank_text`, `cox_text_asy` in the five `cox_*_wv` calls; `fit/wald`, `fit/wald_debiased1`, `z.score_debiased`, `fit_resamples/Z_bstar` in `cox_fh00_d200_wv` (see the execution note). All 16 moved FH(0,0) `sig_bhat_asy` equal `coxph` robust (largest relative difference 1.7e-09) and all 16 `fit/z.score` do (5.4e-14). `golden_compare_L2.txt` | **Added:** `test-l2-weighted-cox-se.R` (4 blocks). **Replaced:** none |
| L3 | DONE | `4c755b5` | `R/weightedCox.R` (`cox_rhogamma_resample()`) | 112 more (684): `fit_resamples/score_star`, `score_star_null`, `bhat_center_star`, `sig_bhat_star`; `fit/sig_bhat_star`, `fit/bhat_debiased`, `fit/wald_debiased2`; `hr_ci_star` (all five columns); `sig2_score_star`; `cox_text_star`, in `cox_fh00_d200_wv`. `golden_compare_L3.txt` | **Added:** `test-l3-weighted-cox-resampling.R` (4 blocks). **Replaced:** none |
| L4 | DONE | `b3057f3` | `R/km_wlr_calculations_helpers.R` (`wlr_dhat_estimates()`) | 80 more (764): `sig2.score`, `z.score`, `sig2_dhat`, `cov_wlr_dhat`, `cor_wlr_dhat` in `wlr_dhat_grid_wv` and `wlr_dhat_fh01_grid_wv` (8 datasets each). `golden_compare_L4.txt` | **Added:** `test-l4-wlr-dhat-weighted.R` (4 blocks). **Replaced:** none |
| B12 | DONE | `53b5748` | `R/kmplotting_helpers.R` (default and help text); `man/plotKM.band_subgroups.Rd` | 88 more (852): `fit_itt/sig2_surv0`, `sig2_surv1`, `sig2_dhat`, `lower`, `upper`, `Zdhat_star`, `c_alpha_band`, `sb_lower`, `sb_upper`, and top-level `lower`, `upper`, in `pkb_wv` (8 datasets). `golden_compare_B12.txt` | **Added:** `test-b12-band-default-draws.R` (3 blocks). **Replaced:** none |
| docs | Step 3 | `8f820bb` | `NEWS.md` (bullets L1, L2, L3, L4, B12; one sentence of the Stage-2c "M, R1" bullet, which said log-rank and Cox do not change, now points to them); roxygen in `R/df_counting.R`, `R/km_wlr_calculations_helpers.R`, `R/weightedCox.R`, `R/kmplotting_helpers.R`; `man/df_counting.Rd`, `wlr_cumulative.Rd`, `wlr_dhat_estimates.Rd`, `cox_rhogamma.Rd`, `cox_score_rhogamma.Rd`, `cox_rhogamma_resample.Rd`, `plotKM.band_subgroups.Rd` | n/a | n/a |

Other commits: `868bda6` (task doc), `857848e` (baseline and copied scripts), and the final commit after `8f820bb`, which adds the verification artifacts and this report.

**What did not move.** In every one of the 852 moved fields the call passes non-constant weights. No field moved in any call without weights or with the constant weight column, and none in `wlr_cum_bare_*` or `wlr_bare_*` (bare tables of counts, including those taken from the weighted fit). `lr`, `score`, `lr_stratified`, `bhat`, `u.beta`, `u.zero`, the time weights, every Kaplan–Meier and RMST field outside `pkb_wv`, counts and risk tables are identical, exactly. No call raises a warning it did not raise before.

**Suite.** 21 test files, 92 `test_that` blocks, 1,036 expectations, 0 failures (inside `R CMD check`). Stage 2c had 16 files, 71 blocks and 613 expectations. No existing test was changed or weakened.

**Checks in the tests that go beyond the contracts.**
- L4: the covariance is checked against numerical derivatives. With FH(0,0) weights, a_i × ∂score/∂a_i and a_i × ∂d̂/∂a_i are computed by central differences for every subject; `sig2.score`, `sig2_dhat` and `cov_wlr_dhat` equal the sums of their squares and products within 1e-6.
- L1: the cumulative V(t) at grid time t equals the `coxph` robust score variance of the same data censored at t (`stage2d/cumulative_check.txt`; 10 time points with and without ties, largest relative difference of z² 5e-14).

## B. In numbers

`stage2d/report_numbers.R` → `report_numbers.txt`. Golden datasets, weights `wv` (uniform on 0.3–3, unrelated to outcome). "2c" is from `golden_s2c.rds`. The value in brackets is `coxph(weights = wv, robust = TRUE, ties = "breslow")`: the signed square root of its robust score test for z, and its robust SE for the SE.

| dataset | n | log-rank z FH(0,0): 2c → 2d (coxph robust) | log-rank z FH(0,1): 2c → 2d | `cox_rhogamma()` SE FH(0,0): 2c → 2d (coxph robust) | SE FH(0,1): 2c → 2d | stratified z FH(0,0): 2c → 2d (coxph robust) | resampling SE FH(0,0), draws = 200: 2c → 2d |
|---|---|---|---|---|---|---|---|
| veteran | 137 | 0.713 → 0.505 (0.505) | 2.172 → 1.620 | 0.1421 → 0.2027 (0.2027) | 0.1680 → 0.2332 | −0.504 → −0.377 (−0.377) | 0.2019 → 0.2085 |
| gbsg | 686 | 4.721 → 3.461 (3.461) | 3.832 → 2.751 | 0.0994 → 0.1371 (0.1371) | 0.1161 → 0.1622 | 4.494 → 3.278 (3.278) | 0.1347 → 0.1326 |
| gbsg_er_le0 | 82 | −2.664 → −1.661 (−1.661) | −1.109 → −0.701 | 0.2375 → 0.3634 (0.3634) | 0.2948 → 0.4605 | – | 0.3315 → 0.3448 |
| gbsg_er_gt0 | 604 | 4.691 → 3.392 (3.392) | 3.158 → 2.241 | 0.1071 → 0.1500 (0.1500) | 0.1231 → 0.1744 | 4.489 → 3.268 (3.268) | 0.1473 → 0.1353 |
| sim_cont | 200 | 2.007 → 1.389 (1.389) | 2.333 → 1.628 | 0.1475 → 0.2146 (0.2146) | 0.1770 → 0.2520 | – | 0.1930 → 0.2321 |
| sim_ties | 200 | 3.304 → 2.413 (2.413) | 2.380 → 1.744 | 0.1406 → 0.1954 (0.1954) | 0.1707 → 0.2379 | – | 0.2026 → 0.1978 |
| sim_weighted | 300 | 3.502 → 2.418 (2.418) | 1.737 → 1.183 | 0.1140 → 0.1643 (0.1643) | 0.1335 → 0.1966 | – | 0.1709 → 0.1580 |
| sim_strata | 300 | 2.201 → 1.579 (1.579) | 1.871 → 1.354 | 0.1168 → 0.1631 (0.1631) | 0.1387 → 0.1914 | 2.359 → 1.614 (1.614) | 0.1575 → 0.1629 |

Reading:
- **z and SE.** Every weighted z-score is smaller in absolute value (ratio 2d / 2c 0.62–0.75) and every weighted SE is larger (ratio 1.38–1.56). The FH(0,0) values equal `coxph` robust to the printed precision in every row, stratified included.
- **Resampling SE.** With 200 draws it was already near the robust SE in Stage 2c (ratio to the 2d analytic SE 0.90–1.04) and is 0.90–1.08 of it in 2d. In Stage 2c the analytic SE was 0.67–0.76 of the resampling SE; the two now estimate the same quantity (contract check with 4,000 draws: ratio within 0.95–1.05).
- **Point estimates.** The score, the stratified score and b̂ for both schemes are identical to 10 digits in every dataset (listed in `report_numbers.txt`).

## C. M6

`stage2d/m6_sim.R` (data model, cells, seeds and `cox_robust()` copied verbatim into `m6_common.R`; nothing under `oc/` was edited) and `m6_summarise.R`. Results: `m6_results.csv`, `m6_l4_results.csv`, `m6_tables.md`, `m6_checks.txt`, `m6_stdout.txt`, `m6_raw/`. N = 400, 2,000 replicates per cell, no failed replicate and no non-finite z in any cell. The resampling columns use the first 500 replicates (`draws = 200`), with the interval centred on the de-biased estimate as `hr_ci_star` is. Monte Carlo error on coverage: ±0.005 at 2,000 and ±0.010 at 500.

| k | cell | scheme | bias (MCSE) | SE/SD | coverage | mean length | resampling, first 500: SE/SD, coverage, length | reject (null cells) | evaluation `proposed`: bias, SE/SD, coverage, length, reject |
|---|---|---|---|---|---|---|---|---|---|
| 1 | U-400 | FH(0,0) | −0.0017 (0.0028) | 1.026 | 0.955 | 0.506 | 0.988, 0.948, 0.499 | – | −0.0017, 1.023, 0.954, 0.505, – |
| 1 | U-400 | FH(0,1) | −0.0005 (0.0034) | 1.021 | 0.958 | 0.602 | – | – | −0.0005, 1.012, 0.955, 0.597, – |
| 2 | R-400 | FH(0,0) | −0.0022 (0.0031) | 0.981 | 0.950 | 0.533 | 1.016, 0.956, 0.536 | – | −0.0022, 0.981, 0.950, 0.533, – |
| 2 | R-400 | FH(0,1) | −0.0031 (0.0036) | 0.990 | 0.951 | 0.629 | – | – | −0.0031, 0.990, 0.951, 0.629, – |
| 5 | D-400 | FH(0,0) | 0.0020 (0.0035) | 1.007 | 0.952 | 0.611 | 0.987, 0.946, 0.609 | – | 0.0020, 1.007, 0.952, 0.611, – |
| 5 | D-400 | FH(0,1) | 0.0004 (0.0043) | 0.976 | 0.945 | 0.731 | – | – | 0.0004, 0.976, 0.945, 0.731, – |
| 8 | E-400 | FH(0,0) | 0.0050 (0.0030) | 1.156 | 0.976 | 0.612 | 1.182, 0.980, 0.609 | – | 0.0050, 1.156, 0.976, 0.612, – |
| 8 | E-400 | FH(0,1) | 0.0099 (0.0040) | 1.027 | 0.953 | 0.728 | – | – | 0.0099, 1.027, 0.953, 0.728, – |
| 10 | F-400 | FH(0,0) | 0.0380 (0.0054) | 0.861 | 0.899 | 0.815 | 0.856, 0.896, 0.809 | – | 0.0380, 0.861, 0.899, 0.815, – |
| 10 | F-400 | FH(0,1) | 0.0441 (0.0063) | 0.829 | 0.875 | 0.917 | – | – | 0.0441, 0.829, 0.875, 0.917, – |
| 14 | D-400-null | FH(0,0) | 0.0103 (0.0034) | 0.992 | 0.944 | 0.598 | 0.970, 0.940, 0.595 | 0.059 | 0.0103, 0.992, 0.944, 0.598, 0.059 |
| 14 | D-400-null | FH(0,1) | 0.0141 (0.0042) | 0.974 | 0.936 | 0.718 | – | 0.068 | 0.0141, 0.974, 0.936, 0.718, 0.068 |
| 16 | F-400-null | FH(0,0) | 0.0386 (0.0061) | 0.782 | 0.877 | 0.832 | 0.808, 0.856, 0.802 | 0.135 | 0.0386, 0.782, 0.877, 0.832, 0.135 |
| 16 | F-400-null | FH(0,1) | 0.0494 (0.0071) | 0.738 | 0.850 | 0.918 | – | 0.166 | 0.0494, 0.738, 0.850, 0.918, 0.166 |
| 12 | Dt-400 | FH(0,0) | 0.0105 (0.0035) | 0.980 | 0.943 | 0.601 | 0.989, 0.948, 0.591 | – | not run |
| 12 | Dt-400 | FH(0,1) | 0.0116 (0.0043) | 0.969 | 0.942 | 0.727 | – | – | not run |
| 21 | Dt-400-null | FH(0,0) | 0.0100 (0.0035) | 0.956 | 0.933 | 0.589 | 0.955, 0.916, 0.579 | 0.069 | not run |
| 21 | Dt-400-null | FH(0,1) | 0.0088 (0.0044) | 0.925 | 0.923 | 0.718 | – | 0.082 | not run |
| 22 | D-400-null-strata | FH(0,0) | – | – | – | – | – | 0.052 | not run |
| 22 | D-400-null-strata | FH(0,1) | – | – | – | – | – | 0.064 | not run |
| 23 | Dt-400-null-strata | FH(0,0) | – | – | – | – | – | 0.059 | not run |
| 23 | Dt-400-null-strata | FH(0,1) | – | – | – | – | – | 0.069 | not run |

**L4, cells 5 and 14** (`wlr_dhat_estimates(dc, tzero = 18)`): mean of each estimate over 2,000 replicates against the empirical value.

| k | cell | scheme | correlation: mean estimated / empirical | var(score): mean estimated / empirical | var(d̂ at 18): mean estimated / empirical | covariance: mean estimated / empirical | largest relative difference of `sig2_dhat` from `KM_diff()` | largest abs. correlation |
|---|---|---|---|---|---|---|---|---|
| 5 | D-400 | FH(0,0) | 0.792 / 0.794 | 348.9 / 335.6 | 0.00356 / 0.00346 | 0.8816 / 0.8562 | 2.2e-15 | 0.878 |
| 5 | D-400 | FH(0,1) | 0.508 / 0.493 | 57.9 / 55.0 | 0.00356 / 0.00346 | 0.2314 / 0.2153 | 2.2e-15 | 0.757 |
| 14 | D-400-null | FH(0,0) | 0.805 / 0.807 | 407.8 / 409.1 | 0.00407 / 0.00407 | 1.0350 / 1.0418 | 1.8e-15 | 0.932 |
| 14 | D-400-null | FH(0,1) | 0.543 / 0.545 | 83.5 / 86.0 | 0.00407 / 0.00407 | 0.3170 / 0.3228 | 1.8e-15 | 0.850 |

**Required checks** (`m6_checks.txt`): both PASS.
- Cells 2, 5, 8, 10, 14, 16, replicates without a tied time after the time fix (1,994–2,000 of 2,000 per cell), both schemes: the package's `sig_bhat_asy` and |z| equal `cox_robust()` within 2.9e-15 (SE) and 1.2e-11 (|z|, from `df_counting()` and from `cox_rhogamma()`); limit 1e-8.
- Cells 12 and 21, all 2,000 replicates (every one has tied times), FH(0,0): `sig_bhat_asy` equals `coxph(robust = TRUE, ties = "breslow")` within 1.4e-09; limit 1e-6.

Reading:
- In cells 2, 5, 8, 10, 14 and 16 the package reproduces the evaluation's `proposed` columns to every printed digit, for both schemes.
- Cell 1 (no weights) takes the unchanged code path; it differs from `proposed` in the third decimal because `proposed` applied the score-residual variance to unweighted data too.
- Tied times (cells 12, 21): coverage 0.943 and 0.933 for FH(0,0), rejection 0.069 under the null; the same design without ties (cells 5, 14) gave 0.952, 0.944 and 0.059.
- Stratified log-rank (cells 22, 23): rejection 0.052 and 0.059 for FH(0,0), 0.064 and 0.069 for FH(0,1). z² equals the stratified `coxph` robust score test within 4.6e-10 in every replicate.
- Estimated weights (cell 8): coverage 0.976 for FH(0,0), SE/SD 1.156. Extreme weights (cells 10, 16): coverage 0.850–0.899, rejection 0.135 and 0.166.
- The resampling SE/SD and coverage on the first 500 replicates are within about two Monte Carlo standard errors of the asymptotic ones on 2,000 in every cell (largest coverage gaps: 0.933 against 0.916 in cell 21, 0.877 against 0.856 in cell 16).
- L4: the estimated correlation is within 0.015 of the empirical one in all four rows, and `sig2_dhat` equals `KM_diff(draws = 0)` to rounding in every replicate.

## D. Vignette impact inventory (no edits)

Method: `stage2d/vignette_impact.R`, a copy of the Stage-2c script. Each version was installed into its own temporary library (START_SHA from the worktree, HEAD from `git archive`), and the three files with weighted calls were knitted with `knitr::knit` under each. Numbers: `stage2d/vignette_numbers_2d.txt`. Figures: `stage2d/vignette_figures_2d.pdf`, 8 pages. The script also re-runs, in the knitted environment, the five `cox_rhogamma()` calls of the chunk at `weightedsurv_examples.Rmd:624`, which overwrites its result object each time.

Files with no evaluated weighted call, as in Stage 2c: `weightedkm_methodology.Rmd`, `weightedrmst_methodology.Rmd`, `weightedcox_wald-vs-logrank_simulations.Rmd` (time weights only), `README.md`.

### D1. Calls

| file:line (chunk header) | call | shows |
|---|---|---|
| `weightedsurv_examples.Rmd:329` | `plot_weighted_km(dfcount_rotterdam_wtd, show.logrank = TRUE)` | subject-weighted log-rank p (text on the plot) |
| `weightedsurv_examples.Rmd:485` | `plot_weighted_km(dfcount_CTC_wtd, show.logrank = TRUE)` | subject-weighted log-rank p (text on the plot) |
| `weightedsurv_examples.Rmd:502` | `wlr_cumulative(dfcount_CTC_wtd, return_cumulative = TRUE)` | cumulative subject-weighted log-rank z (figure) |
| `weightedsurv_examples.Rmd:624` | `cox_rhogamma(dfcount_rotterdam_wtd, draws = 1000)` × 5 schemes | z, HR and resampling interval (gt table) |
| `weightedcox_methodology.Rmd:315` | `cox_rhogamma(dfcount_wtd, draws = 1000)`, FH(0,0) and FH(0,1) | HR and resampling interval (printed) |
| `weightedsurv_examples.Rmd:346` (call 351) | `plotKM.band_subgroups(weight.name = "sw.weights", draws.band = 1000, show_resamples = TRUE)`, no `draws` | resampled curves, pointwise CIs, band (Rotterdam) |
| `weightedsurv_examples.Rmd:369` (call 374) | `plotKM.band_subgroups(weight.name, draws.band = 1000)`, no `draws` | pointwise CIs, band (Rotterdam panel) |
| `weightedsurv_examples.Rmd:520` (call 546) | `plotKM.band_subgroups(weight.name, draws.band = draws)`, no `draws` | pointwise CIs, band (cross-trial data) |
| `martingale_resampling_applications.Rmd:274` | `plotKM.band_subgroups(weight.name, show_resamples = TRUE)`, no `draws` | resampled curves, pointwise CIs, band (two figures) |
| `martingale_resampling_applications.Rmd:316` | `plotKM.band_subgroups(weight.name)`, no `draws` | band (Rotterdam panel) |

No vignette or README call to `plotKM.band_subgroups()` passes `draws`. The weighted `df_counting()` results that feed these chunks are created at `weightedsurv_examples.Rmd:310`, `476`, `702`, `711`; `martingale_resampling_applications.Rmd:196`; `weightedcox_methodology.Rmd:315`.

### D2. Numbers, before → after

Survival estimates, standard errors of the curves, medians and their limits, the `coxph` HR text, the unweighted `survdiff` text and the RMST results are identical before and after in every chunk. b̂ and the HR are unchanged in every `cox_rhogamma()` call.

**Log-rank.**

| chunk | quantity | before | after |
|---|---|---|---|
| examples 310 / 329 / 624 / 702 / 711; martingale 196; coxmeth 315 (Rotterdam) | score (unchanged) 85.04; variance | 208.2 | 599.2 |
| same | z | 5.894 | 3.474 |
| same | text on the plot | "logrank (1-sided) p = <0.001" | "logrank (1-sided) p = <0.001" |
| examples 476 / 485 (cross-trial) | score (unchanged) 37.25; variance | 103.8 | 311.5 |
| same | z | 3.656 | 2.110 |
| same | text on the plot | "logrank (1-sided) p = <0.001" | "logrank (1-sided) p = 0.02" |
| examples 502 (cumulative, cross-trial) | z at 24 / 60 / 120 months | 1.095 / 3.985 / 3.656 | 0.690 / 2.303 / 2.110 |
| same | largest z and its time | 4.088 at 64.8 months | 5.509 at 5.5 months |
| same | grid points with z > 5 (the plot's `ymax`) | 0 | 10 |
| same | first time z exceeds 1.96 | 3.84 months | 2.53 months |
| same | range of z from 80 months on | 3.653–3.656 | 2.108–2.110 |

**`cox_rhogamma()` on the Rotterdam weighted fit, `draws = 1000`** (examples 624; the first two rows are also coxmeth 315).

| scheme | z | asymptotic SE | asymptotic 95% CI of the HR (not displayed) | resampling SE | displayed HR (resampling CI) |
|---|---|---|---|---|---|
| FH(0,0) | 5.894 → 3.474 | 0.0786 → 0.1330 | (0.542, 0.737) → (0.487, 0.820) | 0.1307 → 0.1306 | 0.632 (0.4938, 0.8244) → 0.632 (0.4945, 0.8251) |
| FH(0,1) | 4.188 → 2.308 | 0.0878 → 0.1638 | (0.584, 0.824) → (0.503, 0.956) | 0.1584 → 0.1572 | 0.694 (0.5129, 0.9541) → 0.694 (0.5154, 0.9545) |
| FH_exp2 | 5.502 → 3.152 | 0.0788 → 0.1382 | (0.557, 0.759) → (0.496, 0.853) | 0.1348 → 0.1348 | 0.650 (0.5041, 0.8553) → 0.650 (0.5050, 0.8566) |
| MB(12) | 5.893 → 3.468 | 0.0786 → 0.1333 | (0.542, 0.737) → (0.487, 0.821) | 0.1310 → 0.1308 | 0.632 (0.4936, 0.8249) → 0.632 (0.4944, 0.8255) |
| zero_one(6) | 5.928 → 3.547 | 0.0818 → 0.1372 | (0.527, 0.726) → (0.473, 0.809) | 0.1342 → 0.1329 | 0.619 (0.4807, 0.8134) → 0.619 (0.4825, 0.8125) |

**Band chunks (B12; `draws` not supplied).**

| chunk | quantity | before | after |
|---|---|---|---|
| examples 346 / 369; martingale 274 / 316 (Rotterdam) | band critical value | 3.539 | 3.002 |
| same | SE of δ̂ at 24 / 60 / 120 months | 0.0247 / 0.0512 / 0.0659 | 0.0348 / 0.0496 / 0.0754 |
| same | pointwise CI at 24 | (0.110, 0.207) | (0.090, 0.227) |
| same | simultaneous band at 24 / 60 / 120 | (0.071, 0.246) / (−0.009, 0.353) / (−0.144, 0.322) | (0.054, 0.263) / (0.023, 0.321) / (−0.137, 0.315) |
| examples 520 (cross-trial) | band critical value | 3.630 | 3.051 |
| same | SE of δ̂ at 24 / 60 / 120 | 0.0626 / 0.0638 / 0.0858 | 0.0553 / 0.0589 / 0.0753 |
| same | simultaneous band at 24 / 60 / 120 | (−0.194, 0.260) / (−0.038, 0.426) / (−0.179, 0.444) | (−0.136, 0.202) / (0.014, 0.373) / (−0.097, 0.362) |

δ̂ is unchanged (0.158 / 0.172 / 0.089 Rotterdam; 0.033 / 0.194 / 0.132 cross-trial). The RMST objects of chunks 369, 520 and martingale 316 are identical before and after.

Notes on these numbers:
- **Displayed and not displayed.** The gt table at examples 624 and the printout at coxmeth 315 show z and the resampling interval. z changes by a factor 0.55–0.60; the displayed interval limits change in the third or fourth decimal. The asymptotic interval, which changes most, is not displayed.
- **Asymptotic and resampling SE now agree.** Before, the asymptotic SE was 0.55–0.61 of the resampling SE on this fit; after, 1.02–1.04.
- **Rotterdam log-rank p.** The printed text is the same because p stays below 0.001 at z = 3.474.
- **Cumulative z (examples 502).** Up to 5.52 months the cross-trial data have 44 control-arm events and no treated-arm event. With no event in one arm the sum of squared contributions is small, and the cumulative z reaches 5.51 there; ten grid points lie above the figure's upper limit of 5. The cumulative variance is the `coxph` robust score variance of the data censored at each time (section A).
- **Band sign changes.** The simultaneous lower limit at 60 months changes sign in both datasets: Rotterdam −0.009 → 0.023, cross-trial −0.038 → 0.014.
- **Captured last value only.** As in Stage 2c, the hook captured only the last value of `temp` in chunk examples 369 (the unweighted GBSG panel). Its weighted Rotterdam panel is the same call as martingale 316.

### D3. Figures

`stage2d/vignette_figures_2d.pdf`, one page per figure whose PNG differs, before on the left and after on the right.

**8 figures differ:**
- `weightedsurv_examples.Rmd` chunks at 346, 369, 520 (weighted band plots, B12), 485 (log-rank p text, L1), 502 (cumulative z curve, L1);
- `martingale_resampling_applications.Rmd` chunks at 274 (two figures) and 316 (B12).

**Byte-identical** (not in the PDF): examples 70, 203, 213, 233, 320, 329, 578, 606, 702, 711, 732, 737, 743; martingale 132, 196, 241; coxmeth 305.

The figure at examples 329 is identical although its weighted z changed, because the printed p is "<0.001" in both.

### D4. Sentences or formulas no longer accurate

| file:line | text | why |
|---|---|---|
| `weightedcox_methodology.Rmd:101–103` | "the score variance σ̂²(U_w(β)) is estimated … following the martingale-based approach of @FH1991", with the formula in Ȳ₀, Ȳ₁ and dN̄ | for fits with subject weights that are not all equal (the section at line 311) the score variance is now Σ ψ_i(β)²; the formula shown is the one computed only without subject weights or with equal weights |
| `weightedcox_methodology.Rmd:105` | "at β = 0 the score statistic and its variance reduce to the weighted log-rank statistic and its standard variance estimator" | with subject weights the variance at β = 0 is Σ ψ_i(0)², not the standard estimator |
| `weightedcox_methodology.Rmd:112–114` | "replace the unobservable martingale increments dM_{i,j} with independent standard normal realizations multiplied by observable counting-process increments", with the formula for U†_w | with subject weights each resampled score is −Σ G_i ψ_i(β̂): the multiplier applies to the subject's whole contribution (event term and compensator, censored subjects included) |
| `weightedcox_methodology.Rmd:383` | "computed via the same synthetic-martingale resampling" | as above, when the fits carry subject weights |
| `weightedsurv_examples.Rmd:500` | "Notice how it is flat for time points ≈ ≥ 80." | still holds (range 2.108–2.110 from 80 months); the level it is flat at changed from 3.65 to 2.11, and the curve now leaves the plotted range (`ymax <- 5`, line 505) near 5.5 months |

Checked and not listed:
- `weightedsurv_examples.Rmd:632`, the code comment "with case-weights (here propensity-score weighting) recommend asymptotic versions of SEs", above code that displays the resampling interval. Not a statement about a formula. Before, the two SEs differed by a factor 1.6–1.8 on this fit; they now agree within 4%.
- `weightedcox_methodology.Rmd:107`, "The standard error of β̂_w is then σ̂(U_w(β̂_w)) / I_w(β̂_w)": holds for both variances.
- `weightedcox_methodology.Rmd:313`, `411`, `413`: say that subject weights can be combined with time weights and that both asymptotic and resampling inference exist; still true.
- `martingale_resampling_applications.Rmd:85–87` (W† built from dN·G): describes the unweighted log-rank resampling, which is unchanged; the article runs no weighted log-rank or Cox call. Its lines 303, 419 and 448 were listed in Stage 2c for the Kaplan–Meier resampling and are not changed by this task.
- `weightedsurv_examples.Rmd:341`, "the first graph displays 20 (1st 20) resampled approximations": refers to `show_resamples`, not to `draws`; still true.
- No vignette text mentions the default `draws = 20`.
- The sentences listed in Stage 2c (F4 of that report) are unaffected by this task.

### D5. Weights

Every weighted example uses stabilised inverse-probability-of-treatment weights from a logistic propensity model fitted on the analysis rows, entered as a known weight column.

| example | file:lines | model | form of the weight |
|---|---|---|---|
| Rotterdam (node-positive, hormonal treatment) | `weightedsurv_examples.Rmd:299–304`; `martingale_resampling_applications.Rmd:174–185`; `weightedcox_methodology.Rmd:329–335` | `glm(treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er, family = "binomial", data = df_rotterdam)`; marginal model `glm(treat ~ 1, …)` | `sw.weights` = P̂(treat = 1) / π̂(X) for treated, (1 − P̂(treat = 1)) / (1 − π̂(X)) for controls. Range 0.50–12.88; n = 1,207 / 339; effective sample size 1,163 / 105 (control / treated) |
| Cross-trial (GBSG treated against Rotterdam controls) | `weightedsurv_examples.Rmd:444–450` | `glm(treat ~ age + meno + size + grade3 + nodes + pgr + er, family = "binomial", data = df_CTC)`; marginal model `glm(treat ~ 1, …)` | the same form. Range 0.22–12.98; n = 1,207 / 246; effective sample size 1,099 / 73 |

- Truncated weights (`sw.weights_trunc`, the 5th and 95th percentiles; `weightedsurv_examples.Rmd:306`, `452`; `martingale_resampling_applications.Rmd:188`) are created and not passed to any call.
- `weightedsurv_examples.Rmd:724` fits `glm(hormon2 ~ age + meno + size + grade + nodes + pgr + chemo + er, …)` for `adjustedCurves::adjustedsurv()`, not for a weightedsurv call.

## E. Recorded, not changed

**`draws.band` and forestsearch.**
- The default of `draws.band` in `plotKM.band_subgroups()` is 20 (unchanged). In `KM_diff()` the defaults of `draws` and `draws.band` are 0.
- forestsearch 0.5.0 (CRAN source) calls `weightedsurv::plotKM.band_subgroups()` once, in `R/plot_km_band_forestsearch.R:349`. It passes `draws.band = draws_band` (the wrapper's default is 1000). It does not pass `draws`, and it does not pass `weight.name`. The local development tree (`0.4.1.9003`) has the same call. Without weights, `draws` has no numerical effect.

**Outputs that, for weighted data, are still computed without the weights or compared with an unweighted reference.** Example values are from the golden dataset `sim_weighted` with weights `wv` (`report_numbers.txt`).

| output | what it is for weighted data | example |
|---|---|---|
| `df_counting()$logrank_results` (`chisq`, `pval`, `logrank_text`) | `survival::survdiff()` on the positive-weight rows, without the weights | χ² 10.995, "Logrank p = <0.001"; the weighted z² is 5.849, "logrank (1-sided) p = 0.01" |
| `check_results()` | compares `lr² / sig2_lr` and `z.score²` (the same number) with the unweighted `survdiff` χ², and warns when they differ by more than 1% | 5.849, 10.995, 5.849; warning "Statistics differ by more than 1%" |
| `df_counting()$cox_results$score` | `summary(coxph)$sctest` of the weighted fit: the model-based score test, not the robust one | 12.184; the robust score test of the same fit is 5.849 |
| `df_counting()$cox_results` (`hr`, `hr_ci`, `pval`, `cox_text`) | `coxph(weights = , robust = TRUE)` with the default (Efron) ties; `cox_rhogamma()` uses Breslow sums | HR 0.6723, robust SE 0.1643 in both on this dataset (no ties) |
| the `survdiff` comparison inside `df_counting()` ("Discrepancy with log-rank and survdiff") | runs only when `weight.name` is `NULL` | not run with weights |
| `check.KM` quantile comparison inside `df_counting()` | point estimates always; limits only without weights | unchanged from Stage 2c |
| `plotKM.band_subgroups()` risk table | unweighted numbers of subjects with time ≥ t | 300 / 166 / 88 / 40 at t = 0 / 6 / 12 / 18; weighted numbers at risk 494 / 288.6 / 152 / 67.9 |
| `plotKM.band_subgroups()` `max_tau` | the smaller of the two arms' last event times, from the unweighted data | 23.77 |
| `KM_diff(modify_tau = TRUE)` time range | unweighted quantiles (`qtau`, `1 − qtau`) of the event times | not evaluated here |
| `wlr_cumulative()` on a bare table of counts taken from a weighted fit | the earlier variance formula (by contract: no subject-level data) | golden: `wlr_cum_bare_wv`, `wlr_bare_wv` identical to Stage 2c |

## F. Seed probe, pressure test, API, CRAN check, forestsearch

**Probe** (`dev/audit/km_audit_probe.R --no-write` → `stage2d/probe_2d.txt`). Statuses are unchanged from Stage 2c: OK1–OK3 PASS, I3 the only REPRODUCED, P5 INFO, everything else NOT_REPRODUCED. One line changed text: P5, timing and memory only (0.4 s → 0.3 s).

**Pressure test** (`stage2d/km_pressure_test_s2d.R`, a copy of `stage2c/km_pressure_test_s2c.R` that differs only in its header and default output directory; 85 s; `km_pressure_results.csv`, `km_pressure_stdout.txt`). Failing cells by family:

| family | checked cells | failing Stage 2c | failing Stage 2d |
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

No cell changed `n_checks` or `n_fail`. Of the 612 rows, 2 changed their `worst` value, both in T2 and both at rounding level:
- **T2.1, row permutation: S, SE, quantiles, z (weighted).** Worst difference 4.44e-15 → 3.11e-15. The weighted z is now computed from the per-subject contributions; it is still invariant to the order of the rows.
- **T2.5, swap arm labels: dhat, z, quantile rows.** Worst difference 0 → 8.88e-16. Under the earlier formula the weighted z changed sign exactly when the arm labels were swapped; the sum of squared contributions is now accumulated in a different order for the two labelings.
- **T10.** The same three cells fail as in Stages 2, 2b and 2c, for the reasons given there.
- The pressure test has no check of the weighted log-rank or Cox variance against a reference, so no other row moved. Those references are in the new tests, the golden contract check and M6.

**API** (`Rscript dev/audit/stage2d/api.R compare` against `api_s2c.rds`), verbatim:

```
API compare: 46 exports before, 46 now
DIFFERENCES from api_s2c.rds:
  plotKM.band_subgroups: default of 'draws' changed: 20 -> 0
API: the B12 default is the only difference from api_s2c.rds
```

**CRAN check.** `R CMD check --as-cran --timings weightedsurv_0.2.0.tar.gz`. The tarball came from `R CMD build` of a clean `git clone` at `8f820bb`, the last commit touching `R/`, `man/` or `NEWS.md`. Log: `stage2d/00check_2d.log`.

```
Status: 2 NOTEs
```

- **NOTEs.** The same two local-tooling NOTEs as in Stage 2c: pandoc not on PATH for README/NEWS (top-level files), and old HTML Tidy (HTML version of the manual).
- **Tests.** `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1036 ]`. Examples, `--run-donttest` examples and the vignette rebuild are OK.

**forestsearch 0.5.0** (CRAN source). `R CMD check --no-manual --ignore-vignettes` with `R_LIBS` set to a temporary library holding the HEAD build (`8f820bb`). It ran for 47 s. Logs: `stage2d/forestsearch_00check_2d.log`, `forestsearch_testthat_2d.Rout`.

```
Status: OK
[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]
```

As before, no forestsearch test or example calls weightedsurv.

## G. Gate G2 output

Checks 1–8 were run by `stage2d/g2.sh` (the block from task section 8, verbatim) after the docs commit `8f820bb`, with `START_SHA` and `MASTER_SHA` set and this report in place (`stage2d/g2_output.txt`). The golden compare lines are abridged here; the full list of moved fields is in `golden_compare_B12.txt` and `g2_output.txt`.

```
OK 1a
OK 1b
API compare: 46 exports before, 46 now
DIFFERENCES from api_s2c.rds:
  plotKM.band_subgroups: default of 'draws' changed: 20 -> 0
API: the B12 default is the only difference from api_s2c.rds
OK 2a
OK 2b
final whitelist = union over: L1, L2, L3, L4, B12
golden compare, whitelist 'final': 16829 fields in golden, 16829 in new run, 0 added
WHITELISTED fields that moved (852):
  [... 852 lines, see g2_output.txt ...]
coxph check: 33 moved FH(0,0) z-scores against coxph robust score test (max relative difference of z^2 5.38e-14); 16 moved sig_bhat_asy against coxph robust SE (max relative difference 1.71e-09)
Warnings: no call raises a warning it did not raise in golden_s2c.rds
GOLDEN: identical outside the whitelist
OK 3
OK 4
OK 5
OK 6a
OK 6b
OK 7a
OK 7b
OK 8 (if no FAIL above)
8f820bb stage2d docs: score-residual variance for subject-weighted log-rank and Cox (L1-L4), B12 default, in NEWS and help pages
53b5748 B12: plotKM.band_subgroups() default draws = 0
b3057f3 L4: wlr_dhat_estimates() consistent with L1 for subject weights
4c755b5 L3: cox_rhogamma() resampling from the per-subject score contributions with subject weights
20a3b5f L2: cox_rhogamma() SE and interval from the score-residual variance with subject weights
33b3dfd L1: subject-weighted log-rank uses the score-residual variance
857848e stage2d: baseline (golden_s2c.rds, api_s2c.rds) and copied golden/api scripts
868bda6 stage2d: add task doc
```

Check 9 is run after the final commit; its output is in the session summary and is not part of this file.
