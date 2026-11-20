# weightedsurv 0.2.0: Stage 2e on branch `km-fixes-0.2.0`

**Top line.**
- **Branch and SHAs.** Branch `km-fixes-0.2.0`. START_SHA is `222878dd8fc08050f54d9f7ac1fbb86bd269008b` and MASTER_SHA is `166919ec7128f165fb380847f549fa26ca1db985`. The final SHA is the commit that adds this report together with the verification artifacts; its parent is `c195953` (docs), and no commit follows it.
- **Gates.** G0 and G1 passed. **G2 did not pass: check 9 (M7 requirements) prints FAIL**; checks 1–8 and 10–13 print OK. The output is in section G.
- **Why check 9 fails.** M7 requirement 1 (the package under `"estimated"` against the evaluation's script, cells 7, 8, 9, 15, 19) fails in 3 of 2,500 replicates, in the Cox / log-rank part only (largest relative difference 4.2e-04). In each of the three, the package's time fix turns two raw times about 1e-07 apart into a tie; the script has no handling of tied times and its alignment confirmation does not detect the change. In those replicates the package equals the `coxph()` residual routes within 7e-13. Requirements 2, 3 and 4 are OK. Details in section C.
- **Items.** Rd, Rb, Rc, Ra and W1, W2, W3, W4 are DONE; none is NOT DONE. Each met its contract of §5 at its commit (tests, golden compare, Gate G2 checks 4–8).
- **M7.** Ran (kickoff: `M7: GO`): 15 cells, 22,500 replicates, and the bootstrap benchmark, on 12 cores. Pilot projection 18.5 min; wall-clock 32.1 min (task estimate 35–70 min).
- **Deliverables.** This report; `stage2e/m7_summary.csv`; `stage2e/vignette_known_vs_estimated.pdf`.

Execution notes (findings, not package results):
- **START_SHA.** The kickoff asked for the task document to be copied and committed first. START_SHA is the HEAD found at intake, `222878d`; the task document was committed as `fd2e815`. `R/`, `NAMESPACE`, `vignettes/` and `README.md` are identical at the two commits, and `git diff --quiet 8f820bb fd2e815 -- R/ NAMESPACE` is true. The baselines were made at `fd2e815`; the Step V "before" worktree is of `fd2e815`.
- **Golden script.** `stage2e/golden.R` is a copy of `stage2d/golden.R` with the same datasets and the earlier calls. It adds two weight columns to every dataset (`w01`: 0/1, about 30% zeros; `wz`: `wv` with the same rows set to 0) and these calls: on `gbsg` (tied times, strata), `sim_weighted` and `sim_strata`, with `w01` and with `wz`, `df_counting()` FH(0,0), the same with `strata.name` where the dataset has strata, `KM_diff(draws.band = 200)`, and `cox_rhogamma()` with `draws = 0` and `draws = 200`; on every dataset, `check_results(verbose = FALSE)` without weights and with `wv`, and `plotKM.band_subgroups()` with neither `draws` nor `draws.band`, without weights and with `wv`. Baseline: 18,946 fields from 564 recorded calls, none of which raised a warning.
- **Golden rules changed from Stage 2d.** A list may gain elements anywhere as long as its earlier names keep their order (Stage 2d required the old names to be a prefix; `fit_itt` of `plotKM.band_subgroups()` gains `weight_variance` before `max_tau`). Every new field and every new list name must match the task's list of allowed new elements for the items done, otherwise the compare fails. The Stage 2d `coxph` contract check is dropped from the script (the references are in the tests and in Gate G2).
- **API script.** `stage2e/api.R compare` exits 0 only when the differences from `api_s2d.rds` are exactly those of §5 "API" for the items in `stage2e/items_done.txt`.
- **No existing test encoded the earlier zero-weight behaviour (Rd).** Three existing expectations encoded behaviour this task changes and were replaced; they are listed in section A.
- **Where the code is.** `ps_weights()` and its print method: `R/ps_weights.R`. The internal helpers of the option (validation, the group-level Kaplan–Meier pieces, the adjusted score contributions, the cumulative variance): `R/estimated_weights.R`, all `@noRd`. The Kaplan–Meier part follows the implementation note of the task (only `m × p` matrices, from `ws_km_psi_parts()`, `ws_km_ij_draws()` and `ws_km_ij_var()`).
- **How the option travels.** `df_counting()` returns `weight_variance` and, with `"estimated"`, `ew_parts = list(IF, dw)` in the order of `time`. `wlr_cumulative()`, `wlr_dhat_estimates()` and `cox_rhogamma()` read those two elements. Inside `df_counting()` the pieces go to `wlr_cumulative()` with the subject-level attribute Stage 2d introduced. In a stratified analysis the per-stratum contributions are collected and the variance is one sum over all subjects.
- **One element outside the task's list, present only under `"estimated"`.** `cox_rhogamma()$fit` carries `ew_parts` when `dfcount` was made with `"estimated"`; `cox_rhogamma_resample()` needs it, since it receives only `fit`. No call that does not ask for `"estimated"` has it (golden: no such field).
- **Argument documentation.** The `@param` entries of the two new arguments were added with W2 so that each commit documents its signature; the shared section and the rest of the help text are in the documentation commit.
- **M7 code.** M7, the pressure test, Step V and the checks ran on the `R/` code of `c195953`. `git diff 4035b67 c195953 -- R/` is roxygen comments only (checked: no changed line outside `#'`).
- **Remote branch.** At intake the local reflog showed `origin/km-fixes-0.2.0` at `222878d` ("update by push", before this task). No command of this task pushed. At the end of the task `git ls-remote` still shows the remote branch at `222878d`; every commit of this task is local only. `master` is at MASTER_SHA.
- **Constraints kept.** No merge, no push. Nothing under `vignettes/` and not `README.md` was edited. Installs went only into temporary libraries under the session scratchpad. No package dependency was added. The error messages added (in `ps_weights()` and in the validation of the two new arguments) contain neither `TRUE/FALSE` nor `subscript out of bounds`. No warning or message depends on the effective sample size. Weights are never truncated or rescaled.

---

## A. Items

| ID | status | commit | files touched | golden fields that moved | tests added or replaced |
|---|---|---|---|---|---|
| Rd | DONE | `7a68390` | `R/internal_helpers.R` (new `ws_weights_unequal()`), `R/df_counting.R`, `R/df_counting_helpers.R`, `R/km_wlr_calculations_helpers.R`, `R/weightedCox.R` (the seven tests) | 189, all in the `w01` calls (variances, z-scores, texts, limits, bands, resampled quantities). None in the `wz` calls, none in a call without a zero weight | added `test-rd-zero-weights.R` (4 tests, 77 expectations) |
| Rb | DONE | `d244eb7` | `R/df_counting.R` | 139: `cox_results$score`, `logrank_results$chisq`, `$pval`, `$logrank_text` in the calls with unequal positive weights, and `check_results()` `logrank_chisq`. New: `cox_results$se` and `logrank_results$method` in all 106 `df_counting()` results | added `test-rb-weighted-survival-outputs.R` (4 tests, 40 expectations). **Replaced** in `test-b4-validation.R` (line 63): "weighted `logrank_results$chisq` equals the unweighted `survdiff` chi-square" → equals `coxph(..., robust = TRUE, ties = "breslow")$rscore` on the positive-weight rows, `method = "coxph_robust_score"` |
| Rc | DONE | `e08a27d` | `R/kmplotting_helpers.R`, `man/plotKM.band_subgroups.Rd` | 176, all in the 16 `plotKM.band_subgroups()` calls without `draws.band`: `surv0_star`, `surv1_star`, `dhat_star`, `Zdhat_star` (and their dimensions), `c_alpha_band`, `sb_lower`, `sb_upper` | added `test-rc-band-default-draws-band.R` (2 tests, 6 expectations). **Replaced** in `test-b12-band-default-draws.R`: `formals(...)$draws.band` is 20 → is 1000 |
| Ra | DONE | `556e51a` | `R/df_counting.R` | 0. New: `ess` (5 fields) in all 106 `df_counting()` results | added `test-ra-ess.R` (3 tests, 16 expectations) |
| W1 | DONE | `704c156` | `R/ps_weights.R` (new), `NAMESPACE`, `man/ps_weights.Rd` (new; three other Rd files gain a see-also link) | 0 | added `test-w1-ps-weights.R` (6 tests, 52 expectations) and `helper-ew-ref.R` (the reference implementation of §2, as given) |
| W2 | DONE | `f406510` | `R/estimated_weights.R` (new), `R/km_wlr_calculations_helpers.R` (`KM_diff()`), `R/df_counting.R`, `R/rmst_calculations.R`, `R/kmplotting_helpers.R`, four Rd files | 0. New: `weight_variance` in the `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()` results and in `fit_itt` | added `test-w2-estimated-km.R` (9 tests, 145 expectations). **Replaced** in `test-b8-rmst-mean-rule.R` (line 27): "`mean_rule` is the last argument" → the last three arguments are `mean_rule`, `weight_variance`, `weight_model` |
| W3 | DONE | `0977d68` | `R/df_counting.R`, `R/estimated_weights.R`, `R/internal_helpers.R`, `R/km_wlr_calculations_helpers.R` (`wlr_cumulative()`, `wlr_dhat_estimates()`) | 0 | added `test-w3-estimated-logrank.R` (6 tests, 208 expectations) |
| W4 | DONE | `4035b67` | `R/weightedCox.R` | 0. New: `weight_variance` in the 164 `cox_rhogamma()` results | added `test-w4-estimated-cox.R` (4 tests, 110 expectations) |
| docs | — | `c195953` | `NEWS.md`, `_pkgdown.yml`, roxygen comments in seven `R/` files, eleven Rd files | — | — |

- **Golden, final compare** (`golden.R compare final`, union of the whitelists): 18,946 baseline fields; 504 moved, all whitelisted; 1,122 fields added, all of them listed new elements; no warning that the baseline did not have; nothing outside the whitelist differs at tolerance 0. Per-item compare outputs: `stage2e/golden_compare_<ID>.txt`.
- **Rc, examples and tests that rely on the default.** None. No help-page example calls `plotKM.band_subgroups()`. Every test call passes `draws.band` (`test-b4` 0, `test-b8` 40, `test-b12` 100); the only reliance on the default was the `formals()` expectation above. Run time of the default call on `gbsg` (n = 686): 0.04 s with 20 draws, 0.21 s with 1,000. No example needed an explicit `draws.band`.
- **W2, `plot_weighted_km()`.** Confirmed by reading: `KM_plot_2sample_weighted_counting()` takes its standard errors from the result (`sqrt(dfcount$sig2_surv0)`, `sqrt(dfcount$sig2_surv1)`), so a plot of a result made with `"estimated"` draws the estimated-weight limits. Its log-rank legend is `zlogrank_text` and its Cox legend `cox_results$cox_text`.
- **Tests.** 1,691 expectations, 0 failures, 0 warnings, 0 skips (`stage2e/tests_2e.txt`); 1,036 at the end of Stage 2d.

## B. In numbers

`stage2e/report_numbers.R` → `report_numbers.txt`.

### B1. Riders, on the golden datasets

**Zero weights (Rd).** before (`golden_s2d.rds`) → after | the positive-weight rows analysed with the same weights (no weights for `w01`).

| dataset | weights | zero-weight rows | log-rank z | stratified z | `cox_rhogamma()` SE FH(0,0) | SE of the survival difference, median grid point | band critical value (200 draws) | resampling SE (200 draws) |
|---|---|---|---|---|---|---|---|---|
| gbsg | `w01` | 223 of 686 | 2.915839 → 2.894557 \| 2.894557 | 2.727412 → 2.710705 \| 2.710705 | 0.159712 → 0.156658 \| 0.156658 | 0.049366 → 0.049366 \| 0.049366 | 2.9731 → 2.9286 \| 3.1094 | 0.15524 → 0.16227 \| 0.15545 |
| gbsg | `wz` | 223 of 686 | 2.880339 → 2.880339 \| 2.880339 | 2.687019 → 2.687019 \| 2.687019 | 0.170257 → 0.170257 \| 0.170257 | 0.053951 → 0.053951 \| 0.053951 | 3.0958 → 3.0958 \| 2.9685 | 0.17394 → 0.17394 \| 0.16428 |
| sim_weighted | `w01` | 91 of 300 | 3.011930 → 3.025530 \| 3.025530 | – | 0.180321 → 0.179121 \| 0.179121 | 0.071261 → 0.071261 \| 0.071261 | 3.1138 → 3.0380 \| 2.8901 | 0.18019 → 0.18144 \| 0.16932 |
| sim_weighted | `wz` | 91 of 300 | 2.193960 → 2.193960 \| 2.193960 | – | 0.197632 → 0.197632 \| 0.197632 | 0.079298 → 0.079298 \| 0.079298 | 3.0584 → 3.0584 \| 2.9568 | 0.19064 → 0.19064 \| 0.18428 |
| sim_strata | `w01` | 96 of 300 | 0.956096 → 0.946984 \| 0.946984 | 0.999049 → 0.979497 \| 0.979497 | 0.175923 → 0.177675 \| 0.177675 | 0.072453 → 0.072453 \| 0.072453 | 3.0446 → 2.9905 \| 3.0074 | 0.18392 → 0.17095 \| 0.17159 |
| sim_strata | `wz` | 96 of 300 | 0.503564 → 0.503564 \| 0.503564 | 0.335167 → 0.335167 \| 0.335167 | 0.195778 → 0.195778 \| 0.195778 | 0.080543 → 0.080543 \| 0.080543 | 2.8871 → 2.8871 \| 3.0508 | 0.19850 → 0.19850 \| 0.19287 |

**Weighted datasets (Rb).** `z.score^2` is unchanged by Rb. In brackets: `coxph(weights = , robust = TRUE)`, `$rscore` with `ties = "breslow"` beside the log-rank chi-square and `summary()$robscore` (default ties) beside the Cox score.

| dataset | weights | `z.score^2` | `logrank_results$chisq`: before → after (coxph) | `cox_results$score`: before → after (coxph) | `check_results()` largest relative difference: before → after | ess / n control | ess / n treated |
|---|---|---|---|---|---|---|---|
| veteran | `wv` | 0.2551 | 0.0082 → 0.2551 (0.2551) | 0.4921 → 0.2500 (0.2500) | 187.5% → 1.6e-12% | 0.842 | 0.786 |
| gbsg | `wv` | 11.9752 | 8.5648 → 11.9752 (11.9752) | 22.2401 → 11.9763 (11.9763) | 33.2% → 3.3e-13% | 0.816 | 0.830 |
| gbsg_er_le0 | `wv` | 2.7577 | 4.6935 → 2.7577 (2.7577) | 6.9865 → 2.7577 (2.7577) | 52.0% → 1.6e-14% | 0.830 | 0.816 |
| gbsg_er_gt0 | `wv` | 11.5062 | 12.8826 → 11.5062 (11.5062) | 21.9588 → 11.5060 (11.5060) | 11.3% → 1.1e-13% | 0.795 | 0.826 |
| sim_cont | `wv` | 1.9299 | 2.6571 → 1.9299 (1.9299) | 3.9792 → 1.9299 (1.9299) | 31.7% → 2.3e-12% | 0.834 | 0.795 |
| sim_ties | `wv` | 5.8211 | 6.7954 → 5.8211 (5.8211) | 10.9067 → 5.8089 (5.8089) | 15.4% → 6.1e-13% | 0.804 | 0.837 |
| sim_weighted | `w` (= `wv`) | 5.8488 | 10.9951 → 5.8488 (5.8488) | 12.1843 → 5.8488 (5.8488) | 61.1% → 1.4e-13% | 0.822 | 0.820 |
| sim_strata | `wv` | 2.4930 | 4.0959 → 2.4930 (2.4930) | 4.8108 → 2.4930 (2.4930) | 48.7% → 5.4e-12% | 0.835 | 0.811 |

`check_results()` on `sim_weighted` (weights `w`), printed output:

```
before (START_SHA build):
   zlr_sq logrank_chisq  zCox_sq
 5.848805      10.99509 5.848805
Relative difference: 61.11%
Warning: Statistics differ by more than 1%. Check for computational issues.

after:
   zlr_sq logrank_chisq  zCox_sq
 5.848805      5.848805 5.848805
Relative difference: 1.367e-13%
(no warning)
```

`ess` on `sim_weighted` (weights `w`): treat n = 148, ess 121.31, ratio 0.820; control n = 152, ess 125.01, ratio 0.822.

Reading:
- **Rd.** With a 0/1 column every analytic quantity now equals that of the positive-weight rows without weights to the printed precision (in the tests and Gate G2 check 4: within 1e-10 relative, with and without tied times). The log-rank z moved by up to 0.021 and the Cox SE by up to 2%; the Kaplan–Meier standard error did not move at the printed precision. With zeros among unequal positive weights nothing moved, and the values already equalled those of the subset.
- **Rd, resampled quantities.** They are not equal to the subset's draw by draw, because the multiplier matrices keep a row per subject: band critical value 2.93 against 3.11 (gbsg, `w01`), resampling SE 0.162 against 0.155.
- **Rb.** The weighted `logrank_results$chisq` now equals `z.score^2` and the `coxph` robust score test (Breslow) in every dataset; before it was the unweighted `survdiff` value, between 0.03 and 1.9 times `z.score^2`. `cox_results$score` equals `summary(coxph(robust = TRUE))$robscore`; before it was 1.7–2.5 times that. `check_results()` no longer warns on any weighted dataset.
- **Ra.** The golden weights (uniform on 0.3–3) give ratios 0.79–0.84 in every arm.

### B2. The option

n = 400, propensity model `treat ~ x1 + x2`, data of Gate G2 checks 6 and 7 (`set.seed(11)`). `"known"` → `"estimated"`; in brackets, under `"estimated"`, the `coxph()`-route value (score residuals for z, dfbeta residuals with Breslow ties for the SE) or, for the survival difference, the `n × m` reference of §2.

| tied times | weights | ess / n control, treated | log-rank z FH(0,0) | stratified z FH(0,0) | `cox_rhogamma()` SE FH(0,0) | log-rank z FH(0,1) | SE FH(0,1) | SE of the survival difference, median time | band critical value (1,000 draws) | SE of the RMST difference, last time |
|---|---|---|---|---|---|---|---|---|---|---|
| no | ipw | 0.923, 0.732 | 2.6327 → 2.9575 (2.9575) | 2.6009 → 2.9288 (2.9288) | 0.15556 → 0.14095 (0.14095) | 2.5601 → 2.8150 | 0.18277 → 0.16832 | 0.05532 → 0.05299 (0.05299) | 2.989 → 3.040 | 3.3556 → 3.1121 |
| no | stabilized | 0.923, 0.732 | 2.5342 → 2.8569 (2.8569) | 2.5135 → 2.8386 (2.8386) | 0.15628 → 0.14152 (0.14152) | 2.4201 → 2.6709 | 0.18367 → 0.16911 | 0.05532 → 0.05299 (0.05299) | 2.989 → 3.040 | 3.3556 → 3.1121 |
| yes | ipw | 0.923, 0.732 | 2.6149 → 2.9370 (2.9370) | 2.5785 → 2.9043 (2.9043) | 0.15338 → 0.13900 (0.13900) | 2.5708 → 2.8210 | 0.18229 → 0.16822 | 0.05661 → 0.05432 (0.05432) | 2.824 → 2.925 | 3.3308 → 3.0579 |
| yes | stabilized | 0.923, 0.732 | 2.5204 → 2.8399 (2.8399) | 2.4915 → 2.8129 (2.8129) | 0.15414 → 0.13966 (0.13966) | 2.4324 → 2.6782 | 0.18326 → 0.16913 | 0.05661 → 0.05432 (0.05432) | 2.824 → 2.925 | 3.3308 → 3.0579 |

Reading:
- Under `"estimated"` the log-rank z is 11–13% larger in absolute value, the Cox SE 8–9% smaller, the SE of the survival difference 4% smaller and the SE of the RMST difference 7–8% smaller. The band critical value rises by 0.05–0.10 while the standard error it multiplies falls.
- Every FH(0,0) number equals its `coxph()`-route or reference value to the printed precision, with and without tied times, stratified included.
- `lr`, `bhat`, `dhat` and `rmst_time` are identical under the two options (`identical()` is TRUE in all four rows).
- The Kaplan–Meier columns are the same for `"ipw"` and `"stabilized"`; the log-rank and Cox columns are not, under either option, because the two weight forms weight the arms differently in the pooled risk set.
- `check_results()` on a result made with `"estimated"` agrees: 8.746812 three times (no tied times), 8.626088 three times (tied times).

## C. M7

M7 ran (kickoff: `M7: GO`) with the package at `c195953`: 15 cells, 22,500 replicates, plus the bootstrap benchmark (cells 38 and 39, 500 replicates × 200 resamples). Pilot projection on 12 cores: 17.0 min for the simulation and 1.5 min for the bootstrap (limit 120). Wall-clock: 29.9 min and 2.2 min, with the pressure test, the vignette knits and the checks running beside it. Scripts: `stage2e/oc_common.R`, `oc_sim.R`, `oc_boot.R`, `oc_summarise.R` (adapted copies of the `oc/` scripts, which were not edited). Per-replicate results: `stage2e/raw/` (ignored by git). Outputs: `m7_required.txt`, `m7_summary.csv` (1,523 rows, the full time profile), `m7_tables.md`, `m7_stdout.txt`, `m7_pilot.txt`.

Design notes: cell 40 uses `"ipw"` weights (the task leaves the form open); the truth for RMST under proportional hazards is the closed form of the evaluation's integral (equal to its `integrate()` value: 4.013514 at month 36); in the reproduction cells the weight column is the evaluation's own and the model is passed as `ps_weights(model = )`.

### Required

`stage2e/m7_required.txt`, verbatim:

```
FAIL 1 (cells 7, 8, 9, 15, 19, replicates 1-500; package "estimated" against the evaluation's script: sig2_surv0, sig2_surv1, sig2_dhat at the time points, band critical value and limits, RMST variance, critical value and band limits, sig_bhat_asy and |z| for FH(0,0) and FH(0,1)): 3 of 2500 replicates failing; largest difference 0.000418 (limit 1e-8); the failing replicates are in the Cox / log-rank part only, see m7_req1_diagnosis.txt
OK 2 (cells 31, 32; "stabilized" against "ipw" from the same model, both options: sig2_dhat, band critical value, RMST variance): 0 of 4000 replicates failing; largest difference 4.72e-14 (limit 1e-8)
OK 3 (every cell; dhat, rmst_time, bhat and lr identical under "known" and "estimated"): 0 of 22500 replicates failing; 0 comparisons differ; 15 cells
OK 4 (cell 34, tied times; "estimated" FH(0,0) against the coxph residual routes, Breslow): 0 of 2000 replicates failing; largest difference of z^2 1.56e-10 (limit 1e-8); largest difference of the SE (limit 1e-6): 3.5e-10
```

**Requirement 1 is not met as written: 3 of 2,500 replicates (cell 9 replicate 388; cell 15 replicates 104 and 210).** `stage2e/m7_req1_diagnosis.R` → `m7_req1_diagnosis.txt`:

| cell | replicates | the time fix changes a time | alignment confirmation not passed | Cox / log-rank part differs, alignment passed | of those, no time changed | largest Cox / log-rank difference where no time changed | largest Kaplan–Meier, band, RMST difference |
|---|---|---|---|---|---|---|---|
| 7 | 500 | 0 | 0 | 0 | 0 | 3.6e-12 | 2.1e-14 |
| 8 | 500 | 1 | 1 | 0 | 0 | 1.1e-12 | 3.6e-14 |
| 9 | 500 | 2 | 1 | 1 | 0 | 5.8e-14 | 6.9e-14 |
| 15 | 500 | 2 | 0 | 2 | 0 | 4.1e-12 | 7.3e-14 |
| 19 | 500 | 0 | 0 | 0 | 0 | 2.5e-14 | 3.7e-14 |

- The Kaplan–Meier variances, the band and the RMST quantities agree with the script in all 2,500 replicates (largest difference 7.3e-14).
- `sig_bhat_asy` and |z| agree within 4.1e-12 in the 2,495 replicates in which `survival::aeqSurv()` (the time fix `df_counting()` applies) changes no time.
- In 5 replicates the time fix merges two raw times that differ by 7e-08 to 2.5e-07 into a tie. In 2 of them the script's alignment confirmation does not pass and the requirement excludes them. In the other 3 it passes, because it compares the times with `all.equal()`, and the values differ: SE by 3.4e-07 to 2.4e-05 relative, |z| by 2.5e-05 to 4.2e-04.
- In those 3 replicates the script's `cox_robust()` treats the two times as distinct (it has no handling of tied times), the package treats them as tied. The package equals the tie-aware reference of the tests within 8e-15 and, for FH(0,0), the `coxph()` residual routes (Breslow) on the time-fixed data within 7e-13.
- The same package-against-script difference is present under `"known"` in those replicates (for example cell 15 replicate 104, FH(0,0): |z| 0.939176 package, 0.939020 script), so it does not come from the option.

### Table M7

Nominal 0.95; two-sided 5% log-rank test. Rows: cell × option (analyses A, B, C in cell 40; the two bootstrap rows in cells 38 and 39, first 500 replicates). "(500)" marks 500 replicates (the reproduction cells, shown for completeness). ESS ratio = mean of ess / n, control / treated. Length ratio = mean interval length under `"estimated"` over that under `"known"` (in cell 40, B over C). `—` = not evaluated. Monte Carlo error: ±0.005 at 2,000 replicates, ±0.010 at 500.

| k | Cell | Option | ESS ratio c / t | delta(18) | band [3, 42] | RMST(36) | log HR FH(0,0) | log HR FH(0,1) | SE/SD FH(0,0) | reject FH(0,0) | reject FH(0,1) | length ratio delta(18) | length ratio RMST(36) | length ratio log HR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | E-200 | known (500) | 0.84 / 0.70 | 0.966 | 0.962 | 0.962 | 0.970 | 0.948 | 1.120 | — | — | — | — | — |
| 7 | E-200 | estimated (500) | 0.84 / 0.70 | 0.922 | 0.934 | 0.934 | 0.934 | 0.926 | 0.948 | — | — | 0.869 | 0.828 | 0.847 |
| 8 | E-400 | known (500) | 0.82 / 0.69 | 0.972 | 0.978 | 0.976 | 0.980 | 0.962 | 1.189 | — | — | — | — | — |
| 8 | E-400 | estimated (500) | 0.82 / 0.69 | 0.934 | 0.956 | 0.934 | 0.956 | 0.952 | 1.002 | — | — | 0.863 | 0.820 | 0.843 |
| 9 | E-800 | known (500) | 0.83 / 0.69 | 0.962 | 0.988 | 0.972 | 0.980 | 0.960 | 1.125 | — | — | — | — | — |
| 9 | E-800 | estimated (500) | 0.83 / 0.69 | 0.936 | 0.960 | 0.926 | 0.944 | 0.944 | 0.952 | — | — | 0.866 | 0.821 | 0.847 |
| 15 | E-400-null | known (500) | 0.83 / 0.70 | 0.964 | 0.962 | 0.966 | 0.958 | 0.948 | 1.040 | 0.044 | 0.064 | — | — | — |
| 15 | E-400-null | estimated (500) | 0.83 / 0.70 | 0.938 | 0.950 | 0.928 | 0.940 | 0.914 | 0.879 | 0.070 | 0.092 | 0.870 | 0.826 | 0.845 |
| 19 | E-400-delay | known (500) | 0.84 / 0.70 | 0.980 | 0.984 | 0.984 | — | — | 1.153 | — | — | — | — | — |
| 19 | E-400-delay | estimated (500) | 0.84 / 0.70 | 0.944 | 0.966 | 0.954 | — | — | 0.974 | — | — | 0.865 | 0.824 | 0.845 |
| 31 | Es-400 | known | 0.83 / 0.69 | 0.970 | 0.976 | 0.981 | 0.976 | 0.952 | 1.154 | — | — | — | — | — |
| 31 | Es-400 | estimated | 0.83 / 0.69 | 0.941 | 0.949 | 0.946 | 0.944 | 0.936 | 0.976 | — | — | 0.865 | 0.821 | 0.846 |
| 32 | Es-400-null | known | 0.83 / 0.69 | 0.971 | 0.962 | 0.976 | 0.968 | 0.940 | 1.099 | 0.035 | 0.067 | — | — | — |
| 32 | Es-400-null | estimated | 0.83 / 0.69 | 0.943 | 0.939 | 0.945 | 0.943 | 0.925 | 0.930 | 0.060 | 0.085 | 0.870 | 0.826 | 0.847 |
| 33 | Es-400-null-strata | known | 0.83 / 0.69 | — | — | — | — | — | — | 0.057 | 0.076 | — | — | — |
| 33 | Es-400-null-strata | estimated | 0.83 / 0.69 | — | — | — | — | — | — | 0.065 | 0.081 | — | — | — |
| 34 | Et-400-null | known | 0.83 / 0.69 | 0.973 | 0.963 | — | 0.974 | 0.951 | 1.161 | 0.028 | 0.056 | — | — | — |
| 34 | Et-400-null | estimated | 0.83 / 0.69 | 0.948 | 0.935 | — | 0.940 | 0.935 | 0.979 | 0.060 | 0.070 | 0.867 | — | 0.844 |
| 35 | Fe-400 | known | 0.50 / 0.39 | 0.959 | 0.897 | 0.957 | 0.933 | 0.898 | 1.006 | — | — | — | — | — |
| 35 | Fe-400 | estimated | 0.50 / 0.39 | 0.916 | 0.864 | 0.916 | 0.901 | 0.880 | 0.848 | — | — | 0.850 | 0.822 | 0.843 |
| 36 | Fe-800 | known | 0.47 / 0.33 | 0.963 | 0.920 | 0.951 | 0.936 | 0.909 | 0.946 | — | — | — | — | — |
| 36 | Fe-800 | estimated | 0.47 / 0.33 | 0.925 | 0.892 | 0.903 | 0.894 | 0.892 | 0.807 | — | — | 0.855 | 0.831 | 0.854 |
| 37 | Fe-400-null | known | 0.51 / 0.37 | 0.942 | 0.865 | 0.926 | 0.914 | 0.864 | 0.897 | 0.095 | 0.143 | — | — | — |
| 37 | Fe-400-null | estimated | 0.51 / 0.37 | 0.905 | 0.818 | 0.883 | 0.871 | 0.847 | 0.760 | 0.136 | 0.164 | 0.860 | 0.830 | 0.847 |
| 38 | V-1450-null | known | 0.95 / 0.36 | 0.943 | 0.875 | 0.931 | 0.906 | — | 0.925 | 0.103 | — | — | — | — |
| 38 | V-1450-null | estimated | 0.95 / 0.36 | 0.925 | 0.840 | 0.910 | 0.889 | — | 0.861 | 0.117 | — | 0.933 | 0.916 | 0.931 |
| 38 | V-1450-null | boot (500) | — | 0.906 | — | 0.898 | 0.882 | — | 0.846 | — | — | — | — | — |
| 38 | V-1450-null | boot_pct (500) | — | 0.912 | — | 0.906 | 0.876 | — | 0.846 | — | — | — | — | — |
| 39 | V-1450-ph | known | 0.95 / 0.35 | 0.948 | 0.932 | 0.946 | 0.920 | — | 0.962 | — | — | — | — | — |
| 39 | V-1450-ph | estimated | 0.95 / 0.35 | 0.931 | 0.897 | 0.926 | 0.902 | — | 0.893 | — | — | 0.920 | 0.905 | 0.928 |
| 39 | V-1450-ph | boot (500) | — | 0.932 | — | 0.920 | 0.892 | — | 0.889 | — | — | — | — | — |
| 39 | V-1450-ph | boot_pct (500) | — | 0.932 | — | 0.922 | 0.894 | — | 0.889 | — | — | — | — | — |
| 40 | Sub-800-null | A | 0.86 / 0.95 | 0.948 | — | — | 0.942 | 0.942 | 0.985 | 0.059 | 0.065 | — | — | — |
| 40 | Sub-800-null | B | 0.86 / 0.95 | 0.948 | — | — | 0.951 | 0.941 | 0.986 | 0.053 | 0.066 | 0.957 | — | 0.927 |
| 40 | Sub-800-null | C | 0.86 / 0.95 | 0.958 | — | — | 0.966 | 0.947 | 1.064 | 0.036 | 0.059 | — | — | — |

Bootstrap benchmark on the same first 500 replicates (coverage; SE/SD in brackets):

| k | Cell | Quantity | known | estimated | boot | boot_pct |
|---|---|---|---|---|---|---|
| 38 | V-1450-null | delta(18) | 0.930 (0.96) | 0.912 (0.90) | 0.906 (0.91) | 0.912 (0.91) |
| 38 | V-1450-null | R(36) | 0.926 (0.96) | 0.900 (0.88) | 0.898 (0.89) | 0.906 (0.89) |
| 38 | V-1450-null | log HR | 0.898 (0.89) | 0.884 (0.83) | 0.882 (0.85) | 0.876 (0.85) |
| 39 | V-1450-ph | delta(18) | 0.948 (1.01) | 0.932 (0.94) | 0.932 (0.95) | 0.932 (0.95) |
| 39 | V-1450-ph | R(36) | 0.936 (1.01) | 0.920 (0.91) | 0.920 (0.93) | 0.922 (0.93) |
| 39 | V-1450-ph | log HR | 0.914 (0.94) | 0.890 (0.87) | 0.892 (0.89) | 0.894 (0.89) |

Correlation of the log-rank score with the survival difference at month 18 (FH(0,0), `"estimated"`):

| k | Cell | mean estimated correlation (wlr_dhat_estimates) | empirical correlation over replicates |
|---|---|---|---|
| 31 | Es-400 | 0.708 | 0.722 |
| 32 | Es-400-null | 0.737 | 0.725 |

Failures: no `ps_weights()` error and no failed call in any cell; one non-finite interval (cell 7, one replicate, one late time point).

Reading:
- **Moderate confounding, treated-arm ESS ratio 0.69 (cells 31, 32, 34).** `"known"` covers 0.968–0.981 for the difference at 18, RMST at 36 and the log hazard ratio FH(0,0); `"estimated"` covers 0.940–0.948. Band: 0.962–0.976 → 0.935–0.949. Null rejection, FH(0,0): 0.028–0.035 → 0.060. SE/SD: 1.10–1.16 → 0.93–0.98.
- **FH(0,1) and strata.** FH(0,1) coverage 0.940–0.952 → 0.925–0.936 and null rejection 0.056–0.067 → 0.070–0.085. Stratified log-rank (cell 33): 0.057 → 0.065 for FH(0,0), 0.076 → 0.081 for FH(0,1).
- **Tied times (cell 34)** behave as the untied cell 32: 0.973 → 0.948 for the difference at 18, rejection 0.028 → 0.060.
- **Strong confounding, treated-arm ESS ratio 0.33–0.39 (cells 35–37).** Both options are below nominal and `"estimated"` is lower: difference at 18, 0.942–0.963 → 0.905–0.925; log hazard ratio, 0.914–0.936 → 0.871–0.901; band, 0.865–0.920 → 0.818–0.892; null rejection 0.095 → 0.136 (FH(0,1): 0.143 → 0.164). SE/SD 0.90–1.01 → 0.76–0.85.
- **Vignette-shaped cells, treated-arm ESS ratio 0.35 (cells 38, 39).** Difference at 18: 0.943, 0.948 → 0.925, 0.931. Log hazard ratio: 0.906, 0.920 → 0.889, 0.902. Band: 0.875, 0.932 → 0.840, 0.897. Null rejection 0.103 → 0.117.
- **Bootstrap, same 500 replicates.** `"estimated"` and the refitting bootstrap agree within Monte Carlo error on every quantity (coverage 0.884–0.932 against 0.882–0.932; SE/SD 0.83–0.94 against 0.85–0.95); `"known"` is 0.014–0.026 higher.
- **Interval length.** `"estimated"` intervals are 13–18% shorter than `"known"` in cells 31–37 and 7–10% shorter in cells 38, 39.
- **Subgroup (cell 40).** Full-sample weights with `"known"` (A): 0.948 for the difference at 18, 0.942 for the log hazard ratio, rejection 0.059. Refitted weights with `"estimated"` (B): 0.948, 0.951, 0.053. Refitted weights with `"known"` (C): 0.958, 0.966, 0.036.
- **Correlation of the log-rank score with the difference at 18.** Mean estimate 0.708 and 0.737, empirical 0.722 and 0.725.
- **Reproduction cells (500 replicates each).** Coverage under `"estimated"` 0.922–0.966 across the four quantities (Monte Carlo error ±0.010).

## D. Vignette impact inventory (no edits)

Files: `vignettes/weightedsurv_examples.Rmd`, the five files in `vignettes/articles/`, `README.md`. The three files with weighted calls (`weightedsurv_examples.Rmd`, `articles/martingale_resampling_applications.Rmd`, `articles/weightedcox_methodology.Rmd`) were knitted with `knitr::knit` under the START build and under HEAD, each installed into its own temporary library (`stage2e/vignette_impact.R`, a copy of the Stage 2d script that also records `cox_results$score`, `logrank_results` and `check_results()` output; numbers in `stage2e/vignette_numbers_2e.txt`). `articles/weightedkm_methodology.Rmd` and `articles/weightedrmst_methodology.Rmd` have no evaluated chunk other than setup. `articles/weightedcox_wald-vs-logrank_simulations.Rmd` and `README.md` pass no weights.

### D1. Riders

| class of chunk | chunks found | displayed quantities before → after |
|---|---|---|
| (a) calls `plotKM.band_subgroups()` without `draws.band` | none. All 11 calls pass it: `weightedsurv_examples.Rmd` 350, 373, 396, 524, 546, 588, 608; `martingale_resampling_applications.Rmd` 276, 321, 362, and 423 (`eval=FALSE`). The unevaluated snippet at `weightedrmst_methodology.Rmd:283` passes it too | no change |
| (b) shows `cox_results$score`, `logrank_results` or `check_results()` output for weighted data | none. `check_results(dfcount_gbsg)` (`weightedsurv_examples.Rmd:777`) and `res$logrank_results$chisq` (`:794`) are unweighted | no change: both identical before and after |
| (c) analyses a weight column that contains a zero | none. `sw.weights` ranges over 0.50–12.88 (Rotterdam) and 0.22–12.98 (cross-trial) | no change |

- **Numbers.** 45 recorded objects: 37 identical before and after; 8 changed, each a weighted `df_counting()` result, and each only in the two stored outputs Rb changes. No chunk prints those two outputs, and the plots that use these objects draw `zlogrank_text` and `cox_text`, which did not change.

| objects (chunk at line) | `cox_results$score` | `logrank_results$chisq` | `logrank_results$pval` | `z.score^2` (unchanged) |
|---|---|---|---|---|
| Rotterdam, weighted: `weightedsurv_examples.Rmd` 310, 624, 702, 711, 743; `martingale_resampling_applications.Rmd` 196; `weightedcox_methodology.Rmd` 315 | 34.70 → 12.06 | 4.771 → 12.07 | 0.0290 → 0.00051 | 12.07 |
| cross-trial, weighted: `weightedsurv_examples.Rmd` 476 | 13.34 → 4.438 | 15.43 → 4.453 | 0.0000855 → 0.0348 | 4.453 |

- **Figures.** 26 figure files, all byte-identical before and after. No side-by-side sheet was produced.

### D2. Where the option could apply

`stage2e/vignette_option.R` → `vignette_option.txt`, on the data frames as the vignettes build them (saved from the HEAD knit).

| example | code that makes the weight column | `ps_weights()` with the same formula, `type = "stabilized"`: largest absolute difference from the column | calls that use the column, and whether each analyses exactly the model's rows |
|---|---|---|---|
| Rotterdam (n = 1,546; 339 treated) | `weightedsurv_examples.Rmd:299–304`: `glm(treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er, family = "binomial", data = df_rotterdam)`, marginal model `glm(treat ~ 1, …)`, `sw.weights` = P̂(A) / ê for treated, (1 − P̂(A)) / (1 − ê) for controls | 1.2e-11 | `get_dfcounting()` 313, 703: yes. `plotKM.band_subgroups()` 350, 373 (`df = df_rotterdam`, no subgroups): yes. `cumulative_rmst_bands()` 386: yes. `cox_rhogamma(draws = 1000)`, five schemes, 630–682: yes (takes the `df_counting()` result). `get_dfcounting(check.seKM = TRUE, draws = 5000)` 712: yes for the rows, but `draws > 0` is not available with `"estimated"` |
| cross-trial (n = 1,453; 246 treated) | `weightedsurv_examples.Rmd:444–450`: `glm(treat ~ age + meno + size + grade3 + nodes + pgr + er, family = "binomial", data = df_CTC)`, the same form | 4.1e-08 | `get_dfcounting()` 479: yes. `wlr_cumulative()` 506: yes (takes the result). `plotKM.band_subgroups()` 546: yes. `cumulative_rmst_bands()` 558: yes |
| Rotterdam, martingale article | `martingale_resampling_applications.Rmd:174–185`, the same model and form | 1.2e-11 | `df_counting()` 207: yes. `plotKM.band_subgroups()` 276, 362: yes. `cumulative_rmst_bands()` 389: yes |
| Rotterdam, Cox methodology article | `weightedcox_methodology.Rmd:329–335`, the same model and form | 1.2e-11 | `get_dfcounting()` 338: yes. `cox_rhogamma(draws = 1000)` 345, 349: yes |

- Every model has no missing value in its variables and is fitted on exactly the rows of the data frame analysed. No weighted call analyses a subset of the model's rows (the subgroup plot at `weightedsurv_examples.Rmd:588` is on the unweighted GBSG data).
- The cross-trial difference of 4.1e-08 comes from the marginal model: the vignette takes the arm proportion from `glm(treat ~ 1)$fitted`, `ps_weights()` from `mean(treat)`. `"estimated"` accepts that column (the comparison is `all.equal()` with tolerance 1e-8, a mean relative difference).
- The truncated columns (`sw.weights_trunc`: `weightedsurv_examples.Rmd:306`, `452`; `martingale_resampling_applications.Rmd:188`) are not passed to any call and are not what `ps_weights()` returns.

### D3. Effective sample size (`df_counting()$ess`)

| example | arm | n | ess | ess / n |
|---|---|---|---|---|
| Rotterdam (all three files) | treat | 339 | 105.2 | 0.310 |
| | control | 1,207 | 1,162.9 | 0.963 |
| cross-trial | treat | 246 | 73.0 | 0.297 |
| | control | 1,207 | 1,098.8 | 0.910 |

### D4. Known against estimated

FH(0,0). The band and RMST rows use the vignettes' own `plotKM.band_subgroups()` and `cumulative_rmst_bands()` arguments (`tau_add = 42`, `by.risk = 12`, 1,000 draws). "Median time" is the median of the observed times.

| | Rotterdam, known | Rotterdam, estimated | cross-trial, known | cross-trial, estimated |
|---|---|---|---|---|
| log-rank z (`df_counting()`) | 3.474 | 3.906 | 2.110 | 2.302 |
| one-sided p = 1 − Φ(z) | 0.00026 | 0.000047 | 0.0174 | 0.0107 |
| `logrank_results`: chi-square, p | 12.07, 0.00051 | 15.26, 0.000094 | 4.453, 0.0348 | 5.301, 0.0213 |
| hazard ratio and 95% interval (`cox_rhogamma()`) | 0.6319 (0.4868, 0.8201) | 0.6319 (0.4929, 0.8100) | 0.6601 (0.4445, 0.9803) | 0.6601 (0.4563, 0.9547) |
| SE of the log hazard ratio | 0.1330 | 0.1267 | 0.2018 | 0.1883 |
| `cox_results$cox_text`; `$se` | HR = 0.63 (0.49, 0.82); 0.1331 | HR = 0.63 (0.49, 0.81); 0.1268 | HR = 0.66 (0.44, 0.98); 0.2020 | HR = 0.66 (0.46, 0.96); 0.1885 |
| SE of the survival difference at the median time (45.2; 42.0 months) | 0.04898 | 0.04598 | 0.05734 | 0.05365 |
| band critical value | 3.002 | 3.028 | 3.051 | 3.075 |
| mean width of the band | 0.2490 | 0.2378 | 0.3068 | 0.2910 |
| RMST difference and 95% interval at the last time (118.4; 77.9 months) | 16.82 (7.97, 25.67) | 16.82 (8.54, 25.11) | 8.12 (0.96, 15.28) | 8.12 (1.32, 14.91) |
| RMST band critical value | 2.684 | 2.711 | 2.603 | 2.576 |

- The martingale-article and Cox-methodology Rotterdam examples use the same data and weights as the first column pair and give the same numbers. The Cox methodology article also fits FH(0,1): z 2.308 → 2.438, hazard ratio 0.6935 (0.5031, 0.9559) → 0.6935 (0.5060, 0.9504).
- Every estimate is identical under the two options. Under `"estimated"` the standard errors are 5–7% smaller, the log-rank z 9–12% larger and the bands 4–5% narrower.
- Figures: `stage2e/vignette_known_vs_estimated.pdf`, 3 pages (Rotterdam; cross-trial; Rotterdam in the martingale article). Each page shows the weighted `plotKM.band_subgroups()` plot and the `cumulative_rmst_bands()` plot under `"known"` (left) and `"estimated"` (right), the band plots on a common vertical range.

### D5. Text no longer accurate after this task

| file:line | text | why |
|---|---|---|
| `articles/weightedkm_methodology.Rmd:207` | "The weights are treated as fixed and known: there is no adjustment for estimated weights." | true of the default only; `weight_variance = "estimated"` now carries the estimation of a logistic propensity model into the standard error |
| `weightedsurv_examples.Rmd:773–775` (code comment) | "`dfcount_gbsg` contains logrank (chisq) stats … from 3 sources: … (2) survdiff which only provides rho options" | accurate for the unweighted object it is written about; for a weighted result the second source is now the robust score test of `coxph()`, not `survdiff()` |

Checked and not listed:
- `articles/weightedkm_methodology.Rmd:199`, "the additional variance contribution from the estimated propensity score is typically small": a statement about the cited paper, not about the package. In the examples of D4 carrying the estimation through changes the standard errors by 5–7%.
- `articles/weightedkm_methodology.Rmd:205–207` (the weighted Greenwood formula), `articles/martingale_resampling_applications.Rmd:419` and `articles/weightedrmst_methodology.Rmd:219–225` (weighted resampling): listed in earlier stages for the robust variance and the influence-based resampling; this task does not change their status.
- `articles/weightedcox_methodology.Rmd:424`, "`check_results()` | Verify score-test equivalence": holds, and now holds for weighted results too.
- `weightedsurv_examples.Rmd:632`, the code comment recommending asymptotic standard errors with case weights: not a statement about a formula.
- No vignette text mentions the default of `draws.band`, the effective sample size, or zero weights.
- `articles/martingale_resampling_applications.Rmd:426` and `:436` use `"ps_weights"` as a column name in an unevaluated snippet; the name now coincides with the new function. The snippet is unaffected.

## E. Recorded, not changed

**A constant weight other than 1** keeps the equal-weight formulas (`gbsg`, golden dataset):

| weights | log-rank z FH(0,0) | `sig2_lr` | `cox_rhogamma()` SE FH(0,0) | SE of S₁, median grid point | `logrank_results$method`, chisq | `cox_results$score` | ess (control, treated) |
|---|---|---|---|---|---|---|---|
| none | 2.9265646847 | 70.984135 | 0.1250116764 | 0.0304656687 | survdiff, 8.564781 | 8.566089 | 440.0, 246.0 |
| all 1 | 2.9265646847 | 70.984135 | 0.1250116764 | 0.0304656687 | survdiff, 8.564781 | 8.566089 | 440.0, 246.0 |
| all 2 | 4.1424853192 | 141.714922 | 0.0882994992 | 0.0304656687 | survdiff, 8.564781 | 17.132178 | 440.0, 246.0 |

- With all weights 2 the log-rank z is √2 times the unweighted one and the Cox SE 1/√2 of it; `cox_results$score` doubles (the model-based score test); the Kaplan–Meier standard error, `logrank_results` and `ess` are those of the unweighted data.

Also noticed and left alone:
- **Time fix and near-ties.** `df_counting()` and `KM_diff()` apply `survival::aeqSurv()`. In 5 of the 2,500 reproduction replicates of M7 (continuous times) it merges two raw times that differ by about 1e-7 into a tie. The package then treats them as tied, as `coxph()` does. This is the source of the three replicates of M7 requirement 1 (section C).
- **`cumulative_rmst_bands()`** uses the fixed constants 1.96 (pointwise limits) and the 0.95 quantile (band) under both options.
- **`wlr_cumulative()` and `wlr_dhat_estimates()`** on the result of a stratified `df_counting()` call compute the unstratified statistic, as before, under both options.
- **`check.seKM = TRUE` with `"estimated"`.** The comparison plot is drawn against `survfit()`'s standard error, which treats the weights as known; the two curves differ by design and no discrepancy warning is raised.
- **`logrank_results` under `"estimated"`** is always the `coxph()` route (`method = "coxph_robust_score"`), also when the weights happen to be equal, in line with "the test of whether the weights are all equal is not applied".
- **Rb and the scheme.** With unequal positive weights `logrank_results` is the `rho = 0` test for every scheme; before, `survdiff(rho = scheme_params$rho)` was used for `scheme = "fh"`.
- **forestsearch 0.5.0** calls `plotKM.band_subgroups()` with `draws.band` supplied (its own default is 1000), so Rc does not change its plots; its Kaplan–Meier plots use `show.logrank = TRUE`, which draws `zlogrank_text`.

## F. Seed probe, pressure test, API, CRAN check, forestsearch

**Probe** (`dev/audit/km_audit_probe.R --no-write` → `stage2e/probe_2e.txt`). Statuses are unchanged from Stage 2d: OK1–OK3 PASS, I3 the only REPRODUCED, P5 INFO, the other 25 NOT_REPRODUCED. One line changed text: P5, timing and memory only.

**Pressure test** (`stage2e/km_pressure_test_s2e.R`, a copy of `stage2d/km_pressure_test_s2d.R` that differs only in its header and default output directory; 156 s, run while M7 was using 12 cores; `km_pressure_results.csv`, `km_pressure_stdout.txt`). Failing cells by family:

| family | checked cells | failing Stage 2d | failing Stage 2e |
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

- No change to explain: of the 612 rows, none changed `n_checks`, `n_fail` or its `worst` value. T10 fails in the same three cells as in Stages 2 to 2d, for the reasons given there.
- The pressure test has no cell with a zero weight among other weights that reads a variance, no check of `logrank_results` or `cox_results$score` on weighted data, no `plotKM.band_subgroups()` call without `draws.band`, and no call with `"estimated"`; the riders and the option are covered by the new tests, the golden compare and Gate G2 checks 4–8.

**API** (`Rscript dev/audit/stage2e/api.R compare` against `api_s2d.rds`), verbatim:

```
API compare (items done: Rd,Rb,Rc,Ra,W1,W2,W3,W4): 46 exports before, 47 now
DIFFERENCES from api_s2d.rds:
  cumulative_rmst_bands: appended weight_variance = c("known", "estimated"), weight_model = NULL
  df_counting: appended weight_variance = c("known", "estimated"), weight_model = NULL
  KM_diff: appended weight_variance = c("known", "estimated"), weight_model = NULL
  plotKM.band_subgroups: appended weight_variance = c("known", "estimated"), weight_model = NULL
  export added: ps_weights
  plotKM.band_subgroups: default of 'draws.band' changed: 20 -> 1000
API: the listed differences are the only differences from api_s2d.rds
```

`print.ps_weights` is registered as an S3 method (`S3method(print, ps_weights)` in `NAMESPACE`), not exported by name.

**CRAN check.** `R CMD check --as-cran --timings weightedsurv_0.2.0.tar.gz`, run as in Stage 2d (`RSTUDIO_PANDOC` set to the quarto tools directory, `DISPLAY` unset). The tarball came from `R CMD build` of a clean `git clone` at `c195953`, the last commit touching `R/`, `man/` or `NEWS.md`. Log: `stage2e/00check_2e.log`.

```
Status: 2 NOTEs
```

- **NOTEs.** The same two local-tooling NOTEs as in Stage 2d: pandoc not on PATH for README/NEWS (top-level files), and old HTML Tidy (HTML version of the manual).
- **Tests.** `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]` in 17 s. Examples, `--run-donttest` examples and the vignette rebuild are OK. The `ps_weights()` example runs in 0.2 s.

**forestsearch 0.5.0** (CRAN source). `R CMD check --no-manual --ignore-vignettes` with `R_LIBS` set to a temporary library holding the HEAD build (`c195953`). It ran for 77 s. Logs: `stage2e/forestsearch_00check_2e.log`, `forestsearch_testthat_2e.Rout`.

```
Status: OK
[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]
```

As before, no forestsearch test or example calls weightedsurv.

## G. Gate G2 output

`bash dev/audit/stage2e/g2.sh` with `START_SHA=222878d…` and `MASTER_SHA=166919e…` (the block of §8 as extracted from the task document; full output in `stage2e/g2_output.txt`). **Check 9 prints FAIL**: M7 requirement 1 (section C). Every other line prints OK.

```
OK 1a
OK 1b
API compare (items done: Rd,Rb,Rc,Ra,W1,W2,W3,W4): 46 exports before, 47 now
DIFFERENCES from api_s2d.rds:
  cumulative_rmst_bands: appended weight_variance = c("known", "estimated"), weight_model = NULL
  df_counting: appended weight_variance = c("known", "estimated"), weight_model = NULL
  KM_diff: appended weight_variance = c("known", "estimated"), weight_model = NULL
  plotKM.band_subgroups: appended weight_variance = c("known", "estimated"), weight_model = NULL
  export added: ps_weights
  plotKM.band_subgroups: default of 'draws.band' changed: 20 -> 1000
API: the listed differences are the only differences from api_s2d.rds
OK 2a
OK 2b
final whitelist = union over: Rd, Rb, Rc, Ra, W1, W2, W3, W4
golden compare, whitelist 'final' (Rd,Rb,Rc,Ra,W1,W2,W3,W4): 18946 fields in golden, 20068 in new run, 1122 added
ADDED fields (1122), by element:
  ess@names x 106
  weight_variance x 340
LISTS THAT GAINED ELEMENTS (592):
  + ess, weight_variance in 106 lists
  + method in 106 lists
  + se in 106 lists
  + weight_variance in 274 lists
WHITELISTED fields that moved (504):
  [504 lines, one per moved field, omitted here; all are in stage2e/g2_output.txt]
Warnings: no call raises a warning it did not raise in golden_s2d.rds
GOLDEN: identical outside the whitelist
OK 3
OK 4
OK 5
resampled / analytic variance of the difference (median over the grid): 0.993
OK 6
Cox: resampling SE / analytic SE 1.002; analytic SE estimated 0.1409, known 0.1556
OK 7
OK 8
FAIL: M7 requirements
OK 10a
OK 10b
OK 11a
OK 11b
OK 12 (if no FAIL above)
c195953 stage2e docs: estimated weights (W1-W4), zero weights (Rd), weighted survival-based outputs (Rb), draws.band default (Rc), ess (Ra) in NEWS and help pages
4035b67 W4: cox_rhogamma() with weight_variance = "estimated"
0977d68 W3: log-rank with weight_variance = "estimated"
f406510 W2: Kaplan-Meier variances, bands and RMST with weight_variance = "estimated"
704c156 W1: ps_weights(), propensity-score weights with the influence of the coefficients and the derivative of the weights
556e51a Ra: df_counting()$ess, the effective sample size by arm
e08a27d Rc: plotKM.band_subgroups() default draws.band = 1000
d244eb7 Rb: cox_results$score and logrank_results with subject weights
7a68390 Rd: the weights test looks at the positive weights (zero-weight rows are inert)
356703f stage2e: baseline (golden_s2d.rds, api_s2d.rds, golden.R and api.R for Stage 2e)
fd2e815 stage2e: task document
```

Check 13, after the final commit:

```
OK 13
```
