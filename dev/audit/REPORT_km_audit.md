# weightedsurv — Stage-1 audit and Kaplan–Meier pressure test

**Top line.** No gate failed. G0 passed (clean tree; `pkgload::load_all` OK). G1 passed (exit 0; controls OK1–OK3 PASS; 25 REPRODUCED, 1 NOT_REPRODUCED, 1 INFO). G2 is in section G. M1 was run (kickoff: `M1: GO`). Nothing under `R/`, `NAMESPACE`, `DESCRIPTION` or `vignettes/` was changed. Nothing was installed or pushed.

Audit-execution notes (environment, not package findings):
- The first `R CMD check --as-cran` hung at "checking use of S3 registration": `tcltk_init` was trying to open a Tk frame on the XQuartz `DISPLAY`. It was re-run with `DISPLAY` unset and completed; that run is in `00check.log`.
- pandoc is not on `PATH`. `R CMD build` and the check used `RSTUDIO_PANDOC=/Applications/RStudio.app/.../quarto/bin/tools/aarch64`. `urlchecker` used the same directory prepended to `PATH`.
- The probe's P5 memory figure (21,897,706 "MB") is an artefact of the probe script. In R 4.5 `gc()` has a `limit (Mb)` column, so `gc()[, 6]` is the `max used` cell count, not Mb. T9 reports memory from the last column instead.

---

## A. Provenance

| item | value |
|---|---|
| START_SHA (HEAD before intake) | `edd61a2cecf4b08352719af0a24b10fa1e0b041c`, the commit chat audited |
| branch | `master` |
| status at start | clean |
| R / survival | R 4.5.2 (2025-10-31), aarch64-apple-darwin20, macOS 27.2.1 / survival 3.8.11. Chat used R 4.3.3 / survival 3.5.8. |
| `edd61a2` ancestor of HEAD | yes, and HEAD at start *is* `edd61a2`. `git diff --stat edd61a2 HEAD -- R/ NAMESPACE DESCRIPTION vignettes/` is empty. |
| other packages | all Suggests installed (adjustedCurves, data.table, DiagrammeR, dplyr, gt, knitr, pkgdown, pammtools, rmarkdown); also tibble, urlchecker. qpdf and checkbashisms are missing; HTML Tidy is too old. |

Full detail is in `provenance.txt`.

**Diff against CRAN 0.1.0** (`diff_vs_cran_0.1.0.txt`): `R/` and `NAMESPACE` are byte-identical to the CRAN 0.1.0 tarball, and `DESCRIPTION` has the same version (0.1.0). In `DESCRIPTION`, the field layout is reformatted, `URL` gains the pkgdown site, `Roxygen: list(markdown = TRUE)` is added, and `pkgdown` joins Suggests. In `vignettes/`, `references.bib` (186 lines) became `weightedsurv_references.bib` (2,920 lines) and the one `bibliography:` line in `weightedsurv_examples.Rmd` changed to match. `vignettes/articles/` (5 Rmd) was added and is build-ignored; `.DS_Store`, a git-tracked `.Rhistory` and `.gitignore` are also present. Every functional finding below therefore applies equally to what is on CRAN today.

---

## B. Finding register

Status values: REPRODUCED / NOT REPRODUCED / CHANGED for probe-backed items, CONFIRMED / CORRECTED for S-leads. Line numbers refer to the tree at `edd61a2`. Class is chat's suggestion, or the analogous label for new items; it authorises nothing.

| ID | area | finding | chat (edd61a2) | local tree: status and key number | file:line | class |
|---|---|---|---|---|---|---|
| OK1 | control | Unweighted S and Greenwood SE equal `survfit` | 3e-16 | PASS, 2.2e-16. T1: 7,776 survival and 1,296 unweighted-SE checks, worst 9.6e-15 / 1.1e-16 | — | — |
| OK2 | control | Weighted S equals `survfit(weights=)` | 3e-16 | PASS, 5.6e-16. T1: all 15,292 weighted and unweighted survival checks (KM_diff + df_counting) within 1e-10 | — | — |
| OK3 | control | Weighted SE, natural-scale weights, no ties, vs `survfit` robust | 0.92–1.00 | PASS, ratio 0.923–1.005 | — | — |
| K1 | KM estimate | `modify_tau = TRUE` evaluates the estimator on a grid without event times; it is only KM when `at_points` holds every event time | n=60: 0.43, S −0.31 | REPRODUCED: max error 0.4295, min S −0.3079. T2.3: error 0.517 (time in years) and 0.149 (days), 30/30 datasets each. T10 (via `plotKM.band_subgroups`): 0.040 | `km_wlr_calculations_helpers.R:450-455` (grid); `df_counting_helpers.R:233-235`, `km_wlr_calculations_helpers.R:263-265` (lumped dN and Y on the grid) | B |
| K1b | KM estimate | Same cause through a user `at_points` | 0.23 | REPRODUCED: 0.2283. T2.7 (half the event times removed): 30/30 datasets, error up to 1.62. T8 (grid starting after the 3rd event): 5.7e-4 | `km_wlr_calculations_helpers.R:409` (no check that the grid holds event times) | B |
| K2 | KM quantile | S equal to the level on an interval: package returns the left end, survfit the midpoint; the midpoint rule fires only through floating-point residue | 18/30 | REPRODUCED: 23/30 medians. T4 plateau cases fail at every level (0.9: 6, 0.75: 13, 0.5: 23, 0.25: 12, 0.1: 3 of 30). All 172 point-estimate mismatches in T4 are of this kind | `km_wlr_calculations_helpers.R:78-83` (`qjt1` is strictly > `qprob`, so `round(qjt1,12)==qprob` holds only by rounding) | B |
| K2b | KM quantile | Through `df_counting()` defaults the mismatch is silent | 5.0 vs 5.5 | REPRODUCED: 5.00 vs 5.50, 0 warnings | `df_counting.R:378-391` (signal only when `verbose`) | B |
| K3 | KM variance | Weighted SE changes when weights are rescaled | 0.18/0.64/1.60/6.31 | REPRODUCED: identical 0.18/0.64/1.60/6.31. T2.2: SE(c·w)/SE(w) ranges up to 292 (c=0.01) and 853 (c=1/Σw). M1: normalised weights give SE/SD 0.45 → 3.12 and coverage 0.63 → 1.00 | `df_counting_helpers.R:251-256`: `((dN_w − dN_w/Y_w)/Y_w)^2` is not scale-free | M |
| K4 | KM variance | Weighted SE with tied event times too large | 2.06/1.77/1.52/1.25 | REPRODUCED: identical. T1 heavy ties: median per-arm SE ratio to survfit robust 1.39–1.46 (non-normalised weights), worst 8.5. M1, ties at 0.5: SE/SD 2.31/1.97/1.71, coverage ≈1.000 | `df_counting_helpers.R:256` (squares the summed weighted increment instead of Σw²dN) | M |
| K5 | KM variance | All weights = 2 uses the equal-weight (Greenwood) branch; nearly-2 uses the weighted branch | 0.0297 vs 0.0416 | REPRODUCED: 0.0297 vs 0.0416 (w=1: 0.0420). T2.4: duplicated rows and weight 2 give identical S and SE (both take the Greenwood branch) | `df_counting_helpers.R:237-242` | M |
| K6 | KM estimate | Weighted KM ending on an event leaves ±1e-17 residue; negative values stop `check_km_curve()` | 181/400; 41/60 | REPRODUCED: 168/400 negative; plot stops in 38/60 (probe), 42/60 (T10). T3: of 3,415 KM_diff arms ending on an event, 1,206 are not exactly 0 (|residue| ≤ 9.6e-15); 623/7,776 curves contain values below 0. survP / survG show the same pattern | `km_wlr_calculations_helpers.R:263-265`; stop at `kmplotting_helpers.R:110-118` | B |
| K7 | KM estimate | Unsorted `at_points` accepted | −166 to −38 | REPRODUCED: S0 from −166.59 to −38.00. T8: max error 833, no error or warning | `km_wlr_calculations_helpers.R:409`; negative dN at `km_wlr_calculations_helpers.R:263` | B |
| R1 | robustness | Empty arm in `KM_diff()` | `sig2w_multiplier` not found | REPRODUCED. Same error for one observation (T8), a stratum with one arm empty (T7, all 6 cells) and a `plotKM.band_subgroups` subgroup with an empty arm (T10) | `df_counting_helpers.R:237-242,272` (neither branch runs when `length(W)==0`) | B |
| R2 | robustness | Arm with no events in `df_counting()` | subscript out of bounds | REPRODUCED. Occurs only on the default `check.KM = TRUE` path; with `check.KM = FALSE` a result is returned | `df_counting.R:344-348` (`split_*[[2]]` with one survfit stratum) | B |
| R3a | robustness | `NA` time in `KM_diff()` | TRUE/FALSE error | REPRODUCED | `km_wlr_calculations_helpers.R:441` (`is.unsorted(NA)`) | B |
| R3b | robustness | `NA` time in `df_counting()` | TRUE/FALSE error | REPRODUCED | `df_counting.R:222` | B |
| R3c | robustness | `NA` weight in `KM_diff()` → all-NA survival, no warning | silent | REPRODUCED: 0 warnings; T8 classifies it as silent wrong result. In `df_counting` and `plotKM.band_subgroups` the same input gives a cryptic error | `km_wlr_calculations_helpers.R:437` (no weight validation) | B |
| R4a | robustness | Negative weights accepted by `KM_diff()` | silent | REPRODUCED: no error, no warning (T8: silent wrong result) | `km_wlr_calculations_helpers.R:437` | B |
| R4b | robustness | Zero weights in `df_counting()`: KM correct, self-check compares unequal-length vectors | false discrepancy | REPRODUCED: KM correct, 1 false "Discrepancy in KM curve fit" plus a recycling warning. With default `get.cox = TRUE`, `coxph` rejects zero weights outright (N6) | `df_counting.R:339-348,402` | B |
| R5 | robustness | `draws.band < 20` with default `show_resamples` | subscript out of bounds | REPRODUCED: fails for 1, 5, 19; 20 and 21 OK | `km_wlr_calculations_helpers.R:502-505` (`Zdhat_star[, 1:20]`) | B |
| R6 | robustness | Stratified `df_counting()` without `library(survival)` | `could not find function "strata"` | NOT REPRODUCED: works with and without survival attached (probe; T7, 18 cells). Cause is the dependency, not the tree: survival 3.7-3 added an internal `strata()` for `coxph` formulas, and 3.8-1 extended this to `survdiff`/`survfit`/`concordance`. Under survival 3.5.8 (chat) it fails. `DESCRIPTION` imports survival with no version floor | `df_counting.R:158,195` (unimported `strata()` in formulas); `DESCRIPTION` Imports | B |
| R7 | robustness | Event coded 1/2: `df_counting()` accepts, `KM_diff()` rejects | inconsistent | REPRODUCED. T8: `df_counting` gives the correct KM (Surv treats 1/2 as 0/1); `KM_diff` and `plotKM.band_subgroups` give an informative error | `km_wlr_calculations_helpers.R:424` vs `df_counting.R:142` | B |
| I1 | inference | Band critical value hard-coded at the 0.95 quantile; `alpha` ignored | 2.9796 | REPRODUCED: 2.9796 for α 0.05 and 0.20. T5: identical c_alpha for α ∈ {0.01, 0.05, 0.10, 0.20} in 5/5 datasets | `km_wlr_calculations_helpers.R:500` | B |
| I2a | inference | `wlr_dhat_estimates()`: pooled S read one grid point after `tzero` | 0.947 | REPRODUCED: ratio 0.9467 | `km_wlr_calculations_helpers.R:208,214` | B |
| I2b | inference | `tzero` past the last time gives `dhat` of length 0 | silent | REPRODUCED. T8: `tzero` before the first time also gives length 0 (expected 0) | `km_wlr_calculations_helpers.R:208-213` | B |
| I3 | inference | `cumulative_rmst_bands()` uses the trapezoid rule on a step function | 4.2836 vs 4.2688 | REPRODUCED: 4.2836 vs exact 4.2688 | `rmst_calculations.R:53-56,61-64,100-103` | B (moves every RMST value) |
| I4 | inference | Quantile cross-check always uses survfit's 0.95 limits | false at 0.90 | REPRODUCED: 0 warnings at 0.95, 2 at 0.90 | `df_counting.R:353,365` | B |
| I5 | inference | `qprob ≠ 0.5`: `quantile_results` uses "S drops to qprob", `quantile_check` asks survfit for the `qprob`-th percentile | 15.54 vs 4.29 | REPRODUCED: 15.54 vs 4.29. T4 states the convention: `qprob` means the time at which S drops to `qprob` | `km_wlr_calculations_helpers.R:75` vs `df_counting.R:352-353` | B |
| P2 | package | `set.seed()` inside functions overwrites the caller's RNG state | confirmed | REPRODUCED: `runif` after seeds 123 and 999 both give 0.053202. T5: `.Random.seed` changed | `km_wlr_calculations_helpers.R:483`; `df_counting_helpers.R:261`; `weightedCox.R:417` | B |
| P5 | package | O(n²) matrices in `df_counting()` | 4,000: 3.0 s / 540 MB | INFO: probe n=2,000 0.3 s (its MB figure is an artefact, see top line). T9 `df_counting()`: 0.02 s / 0.20 s / 1.05 s / 3.71 s at n = 500 / 2,000 / 4,000 / 8,000. Peak session memory: 227 / 306 / 634 / 1,884 MB, including about 220 MB of baseline | `df_counting_helpers.R:233-234,244-245` (`outer()` n × grid); `df_counting.R:237-248` (grid = all n times) | — |
| S1 | weights | `MB`: `mb_tstar` = last time → all weights 1; beyond → `1/S_left`, infinite where `S_left=0` | lead | CONFIRMED: unique weight = 1 at `mb_tstar == max(tpoints)`; an infinite weight appears when `S_left` reaches 0 and `mb_tstar` is beyond | `weighted_helpers.R:114-117` | B |
| S2 | weights | `requires_cols` typo means the `survG` requirement for Schemper is never enforced | lead | CONFIRMED, with a correction to the effect: the missing column is not silently tolerated. `dfcounting$survG` partially matches `survG_all` (length n), so `data.frame()` fails with the cryptic "arguments imply differing number of rows: 125, 120" | `km_wlr_calculations_helpers.R:180,192,195` and `541,553,556` | B |
| S3 | KM estimate | Default `at_points = sort(df[[tte]])` uses raw times; estimation uses `aeqSurv`-fixed times | lead | CONFIRMED in the source. With a 9.4e-12 near-tie pair both raw values are in the grid, but the S0 values were identical in the probe; no numeric effect observed | `km_wlr_calculations_helpers.R:409,433-435` | B |
| S4 | RMST | 1.96 and the 0.95 quantile hard-coded; `draws.band = 0` fit fails inside `apply()` | lead | CONFIRMED: "dim(X) must have a positive length" (T10). Also: `fit$risk.points` is read, but `KM_diff()` never returns it | `rmst_calculations.R:49,61,70-71,110` | B |
| S5 | plotting / quantile | Log-transform limits not truncated at 1; `NaN` where S = 0 | lead | CONFIRMED: T10 has 3 NaN limit values and the upper limit > 1 at 4 points (max 1.02). Negative K6 residue makes `log()` warn "NaNs produced" inside `km_quantile` on ordinary weighted data (T8) | `kmplotting_helpers.R:18-19`; `km_wlr_calculations_helpers.R:111-112` | B |
| S6 | plotting / RMST | `plotKM.band_subgroups(modify_tau=TRUE)` reuses the grid for subgroups; `cumulative_rmst_bands` reuses `fit$at_points`, so K1 propagates | lead | CONFIRMED: fit_itt error 0.040 vs survfit (T10); subgroup calls take `at_points = fit$at_points` | `kmplotting_helpers.R:462-472,512-516`; `rmst_calculations.R:46,91-96` | B |
| S7 | package | `get_validated_weights()` uses `cat()` when `mb_tstar` is missing | lead | CONFIRMED: prints "Missing mb_tstar argument in scheme_params you have:", then errors from `validate_scheme_params` | `weighted_helpers.R:253` | B |
| S8 | robustness | `df_counting()` does not validate `strata.name` or the event coding | lead | CONFIRMED: "object 'no_such_column' not found"; an event coded 3 gives "missing value where TRUE/FALSE needed" | `df_counting.R:142,194,217` | B |
| S9 | resampling | `calculate_risk_event_counts(draws>0)` reseeds with the same `seedstart` for every arm and stratum | lead | CONFIRMED: output identical whatever the caller's seed; `df_counting` passes one `seedstart` to every call | `df_counting_helpers.R:261`; `df_counting.R:263,284,487,495,507,513` | B |
| S10 | housekeeping | Inner `check_km_curve()` shadows the export; dead branch in `wlr_dhat_estimates`; `ans$strata` assigned twice; `pkgdown` in Suggests | lead | CONFIRMED: inner definition `df_counting.R:392` (signature differs from `kmplotting_helpers.R:110`); `km_wlr_calculations_helpers.R:209` (`at_points[loc] <= tzero` cannot hold when `loc` is the first index > `tzero`); `df_counting.R:234,450`; `pkgdown` used nowhere in R/, vignettes or man. Also unused `dJ` and `n` at `df_counting_helpers.R:243,253` | as listed | — |
| N1 | robustness | Default `df_counting()` (`check.KM=TRUE`) fails when there are no events at all | new | "first argument must be a vector" (T8 all-censored; T6) | `df_counting.R:341-342` (`rle()` on NULL strata) | B |
| N2 | robustness | Unweighted `df_counting()` fails when the log-rank z is NaN (zero variance: n = 1 per arm, all censored) | new | "missing value where TRUE/FALSE needed". Part of the 130 T1 `df_counting` errors, all at n ≤ 2 per arm | `df_counting.R:558-559` (`if (abs(z^2 − chisq) > tol)` with NaN) | B |
| N3 | robustness | `df_counting()` fails in `survdiff` on degenerate data, also with `get.cox = FALSE` and with weights | new | "Lapack routine dgesv: system is exactly singular" (n = 1 per arm, tied times) | `df_counting.R:193-198` (unconditional `survdiff`) | B |
| N4 | KM quantile | Quantile CI limit differs from survfit when the log-scale limit curve is non-monotone (package takes the first crossing) | new | T4: 5/3,000 upper-limit comparisons (0/3,000 lower), all 5 on non-monotone curves; e.g. 27.60 vs 44.24 | `km_wlr_calculations_helpers.R:75,111-116` | B |
| N5 | Cox | `cox_rhogamma()` with no events in one arm returns a divergent root silently | new | bhat = −1243.65, no warning. The `get_Cox$root` crash path for `NA` (`find_cox_root` returns `NA`) was not reached | `weightedCox.R:64-70,287-294` | B |
| N6 | robustness | Default `df_counting()` (`get.cox=TRUE`) rejects zero weights through `coxph` | new | "Invalid weights, must be >0". `KM_diff` accepts zero weights and gives the correct KM | `df_counting.R:161-162` | B |
| N7 | robustness | `plotKM.band_subgroups()` fails when an arm has no events (also empty arm, all censored, one observation) | new | "'to' must be a finite number" (T8) | `kmplotting_helpers.R:440-442,456` (`max()` of empty → −Inf → `seq()`) | B |
| N8 | inference | In weighted `df_counting()`, `logrank_results` (chisq, p, `logrank_text`) are from unweighted `survdiff` | new | Static: `survdiff` is called without weights, as the comment at line 188 says | `df_counting.R:186-207` | — |
| N9 | robustness | Negative times accepted silently by `KM_diff`, `df_counting`, `plotKM.band_subgroups` | new | T8: silent wrong result in all three | no time validation in the three functions | B |

**Reading.**
1. With a full set of event times on the grid, the core estimator is exact: 15,292 survival checks agree with `survfit` within 1e-14, including all ties, weights and censoring levels. Unweighted Greenwood SE is exact too.
2. Every chat finding reproduces on this tree except R6. R6 is a survival-version effect (≥ 3.8-1 resolves `strata()` itself); it returns for anyone on survival < 3.7-3.
3. All failing weighted-SE behaviour comes from one line (`df_counting_helpers.R:256`) plus the branch on `length(unique(W))`. M1 shows the documented formula and `survfit` robust are calibrated (SE/SD 0.95–1.00, coverage 0.93–0.95). The package is too wide with ties (coverage ≈1.00) and wrong in both directions with normalised weights (coverage 0.63–1.00).
4. K1 makes `modify_tau = TRUE`, any user grid and every `plotKM.band_subgroups(modify_tau=TRUE)` fit non-KM. The size of the error depends on the time unit, because `tau.seq = 0.25` is in data units.
5. Many robustness failures sit in the self-check code of `df_counting()` (R2, R4b, N1, N2) rather than in the estimator.

---

## C. Pressure test (`km_pressure_test.R` → `km_pressure_results.csv`, `km_scaling.csv`)

Seeded. Two consecutive runs gave byte-identical `km_pressure_results.csv`; the full run takes about 125 s. "Cells" counts rows in the CSV; report-only rows have `n_fail = NA`.

| family | cells run | pass | fail | report-only | checks (failed) | worst discrepancy |
|---|---|---|---|---|---|---|
| REF (reference lookup sanity) | 1 | 1 | 0 | 0 | 12 (0) | 0 |
| T1 agreement with survfit (3,888 datasets) | 141 | 101 | 10 | 30 | 41,178 (130) | survival 9.6e-15 (pass); unweighted SE 1.1e-16 (pass); 130 `df_counting` errors (N2/N3); weighted-SE ratio to survfit robust from 0.001 to 1,010 |
| T2 invariants | 18 | 11 | 6 | 1 | 513 (93) | K3 SE ratio 853; K1 0.517; T2.7 1.62 |
| T3 curve shape (from T1, plus survP/survG) | 288 | 204 | 84 | 0 | 84,883 (5,205) | residue below 0 to −8.4e-15; non-zero ending ≤ 1.3e-14 (K6); non-increasing, S=1 before first event and carry-forward: 0 failures |
| T4 quantiles (200 datasets + 30 plateaus) | 22 | 7 | 14 | 1 | 7,401 (281) | point estimate 4.02 (K2); upper limit 16.6 (N4); lower limit 0 failures |
| T5 bands and resampling | 11 | 4 | 5 | 2 | 401 (9) | c_alpha constant at 3.0211 (I1); R5 errors |
| T6 df_counting consistency | 24 | 6 | 18 | 0 | 888 (231) | numeric checks ≤ 4.3e-14 (all pass); 231/648 clean datasets raise warnings or errors with `check.KM = TRUE` |
| T7 stratified (with and without survival attached) | 18 | 12 | 6 | 0 | 138 (6) | `*_stratum` exact (≤ 3.3e-16); 6 failures = one arm empty in a stratum (R1) |
| T8 input matrix | 83 | 47 | 26 | 10 n/a | 73 (26) | 17 cryptic errors, 9 silent wrong results |
| T9 scaling | 1 | 1 | 0 | 0 | 12 (0) | n = 8,000: `cox_rhogamma` 10.3 s / 2.3 GB; no cell near 6 GB |
| T10 plot paths | 5 | 0 | 5 | 0 | 64 (46) | 42/60 `plot_weighted_km` return NULL (K6) |

T6 detail: `surv0/1` equal `KM_diff` (0); `riskpoints0/1` equal `survfit` `n.risk`, weighted and unweighted (4.3e-14); `survP` and `survG` equal `survfit` (7.8e-16); `z^2` and `lr^2/sig2_lr` equal `survdiff` chisq (1.7e-14). The 231 `check.KM` failures split into false "outside [0,1]" warnings (K6), false "Discrepancy" warnings (zero weights, R4b), and errors (R2, N1).

T8 (outcome per input × function):

| input | KM_diff | df_counting (defaults) | plotKM.band_subgroups |
|---|---|---|---|
| empty arm | cryptic (R1) | informative (survdiff "only 1 group") | cryptic (N7) |
| arm with no events | correct | cryptic (R2) | cryptic (N7) |
| all censored | correct | cryptic (N1) | cryptic (N7) |
| one observation | cryptic (R1) | cryptic ("'x' must be an array…") | cryptic (N7) |
| NA time | cryptic (R3a) | cryptic (R3b) | cryptic |
| NA event | informative | cryptic | informative |
| NA treatment | informative | informative | informative |
| NA weight | silent wrong (R3c) | cryptic | cryptic |
| NA strata | n/a | cryptic (`vapply` length 0, `df_counting.R:535`) | n/a |
| Inf time | correct | cryptic ("NA/NaN/Inf in foreign function call") | correct |
| zero time | correct | correct | correct |
| negative time | silent wrong (N9) | silent wrong (N9) | silent wrong (N9) |
| zero weights | correct | informative (`coxph`, N6) | correct |
| negative weights | silent wrong (R4a) | informative (`coxph`) | informative |
| treatment 1/2 | informative | informative | informative |
| treatment factor / logical / character | correct | correct | correct |
| event 1/2 | informative (R7) | correct (R7) | informative |
| event logical | correct | correct | correct |
| unsorted `at_points` | silent wrong (K7, error 833) | n/a | n/a |
| duplicated `at_points` | correct | n/a | n/a |
| `at_points` beyond last time | correct | n/a | n/a |
| `at_points` not starting at/before first event | silent wrong (K1b, 5.7e-4) | n/a | n/a |
| missing column name | informative | informative | informative |
| tibble / data.table | correct | correct | correct |
| `tzero` before first / after last (`wlr_dhat_estimates`) | silent wrong (length-0 `dhat`, I2b) | | |

Every `df_counting` "correct" row on the weighted base data also carries a false "KM curve has values outside [0,1]" warning and "NaNs produced" (K6/S5): the base dataset's control arm ends on an event.

Other measured, report-only results:
- T1 weighted SE, per-arm median ratio package / survfit robust, by weights: no ties 0.90–0.96 (normalised 1.84); heavy ties 1.39–1.46 (normalised 2.47). Package / documented formula: no ties 0.97–1.00 (normalised 1.95); heavy ties 1.59–1.63 (normalised 2.85). Extremes: 0.001 (n = 2 per arm, moderate ties); 10.8 (heavy ties, n = 500); about 1,000 (normalised weights, n = 500).
- T2.8: SE(`draws=2000`) / SE(`draws=0`), weighted, quartiles 0.36 / 0.88 / 1.01 (worst 0.27).
- T5: var(`surv0_star`) / `sig2_surv0` along the curve: unweighted 1.00 → 0.92 (10th → 90th percentile of the curve; worst 0.77); weighted 0.99 → 1.05 (worst 1.21).
- T9 (`km_scaling.csv`): `KM_diff(draws.band=1000)` takes 0.05 / 0.42 / 1.22 / 2.85 s and `cox_rhogamma(draws=0)` 0.05 / 0.67 / 2.63 / 10.26 s at n = 500 / 2,000 / 4,000 / 8,000. Peak memory at n = 8,000: 1.5 GB and 2.3 GB.

Failures not already in the register: none. Everything that failed in T1–T10 maps to a register ID, including the new N1–N9.

---

## D. M1 — weighted-SE calibration (`km_m1_montecarlo.R` → `km_m1_results.csv`)

Design: 2,000 replicates per cell, n = 200, one arm, T ~ Exp(0.08), C ~ U(0, 30). Runtime 18 s. Each entry is SE/SD (mean SE over the empirical SD of Ŝ), then coverage of Ŝ ± 1.96·SE for exp(−0.08t).

| ties | weights | estimator | t = 5 | t = 10 | t = 15 |
|---|---|---|---|---|---|
| continuous | none | package | 1.01 / 0.948 | 0.98 / 0.942 | 1.00 / 0.943 |
| continuous | none | survfit robust | 1.01 / 0.948 | 0.98 / 0.942 | 1.00 / 0.943 |
| continuous | none | documented | 1.00 / 0.947 | 0.97 / 0.942 | 0.99 / 0.939 |
| continuous | U(0.3,3) | package | 0.99 / 0.941 | 0.96 / 0.940 | 0.97 / 0.937 |
| continuous | U(0.3,3) | survfit robust | 1.00 / 0.944 | 0.98 / 0.945 | 0.99 / 0.943 |
| continuous | U(0.3,3) | documented | 0.99 / 0.943 | 0.97 / 0.942 | 0.98 / 0.939 |
| continuous | U(0.3,3)/Σ | package | **0.45 / 0.628** | **1.29 / 0.974** | **3.12 / 1.000** |
| continuous | U(0.3,3)/Σ | survfit robust | 1.00 / 0.944 | 0.98 / 0.945 | 0.99 / 0.943 |
| continuous | U(0.3,3)/Σ | documented | 0.99 / 0.943 | 0.97 / 0.942 | 0.98 / 0.939 |
| ceiling 0.5 | none | package | 1.01 / 0.945 | 0.98 / 0.945 | 1.00 / 0.947 |
| ceiling 0.5 | none | survfit robust | 1.01 / 0.945 | 0.98 / 0.945 | 1.00 / 0.947 |
| ceiling 0.5 | none | documented | 0.98 / 0.940 | 0.95 / 0.937 | 0.97 / 0.938 |
| ceiling 0.5 | U(0.3,3) | package | **2.31 / 0.999** | **1.97 / 1.000** | **1.71 / 1.000** |
| ceiling 0.5 | U(0.3,3) | survfit robust | 1.00 / 0.942 | 0.98 / 0.946 | 0.99 / 0.944 |
| ceiling 0.5 | U(0.3,3) | documented | 0.97 / 0.934 | 0.95 / 0.940 | 0.96 / 0.934 |
| ceiling 0.5 | U(0.3,3)/Σ | package | **0.88 / 0.886** | **2.18 / 0.993** | **4.32 / 1.000** |
| ceiling 0.5 | U(0.3,3)/Σ | survfit robust | 1.00 / 0.942 | 0.98 / 0.946 | 0.99 / 0.944 |
| ceiling 0.5 | U(0.3,3)/Σ | documented | 0.97 / 0.934 | 0.95 / 0.940 | 0.96 / 0.934 |

Reading:
- Unweighted, and weighted with natural-scale weights and no ties, the package SE is calibrated (SE/SD 0.96–1.01, coverage 0.937–0.948).
- With tied times and weights, the package SE is 1.7–2.3× the true SD (coverage ≈1.00). With normalised weights it is wrong in both directions: too small early (coverage 0.63), too large late.
- `survfit` robust is calibrated in every cell. The documented formula is invariant to weight scale but runs slightly small with ties (SE/SD 0.95–0.97, coverage 0.934–0.940).

---

## E. CRAN-readiness baseline (changes nothing)

Commands: `R CMD build .` in a scratch directory (vignettes built), with `RSTUDIO_PANDOC` set to the RStudio-bundled pandoc. Then `R CMD check --as-cran --timings weightedsurv_0.1.0.tar.gz` with `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE` and `DISPLAY` unset. No Suggests were missing, so `_R_CHECK_FORCE_SUGGESTS_`, `--no-build-vignettes` and `--ignore-vignettes` were not needed. Check time was 47 s elapsed (build 15 s). The tarball was not committed.

**Status: 1 WARNING, 2 NOTEs** (`00check.log`), verbatim:

```
* checking CRAN incoming feasibility ... [3s/15s] WARNING
Maintainer: ‘Larry Leon <larry.leon.05@post.harvard.edu>’

Insufficient package version (submitted: 0.1.0, existing: 0.1.0)
```
```
* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.
```
```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
```
Also: `checking installed package size ... INFO installed size is 7.8Mb; doc 7.5Mb`. Both NOTEs come from local tooling (pandoc not on `PATH`, old HTML Tidy).

Five slowest examples (elapsed s): `KM_diff` 0.460, `create_baseline_table` 0.041, `weightedsurv-package` 0.030, `df_counting` 0.026, `cox_rhogamma` 0.025. Vignette rebuild took 10 s.

Reverse dependencies: `forestsearch` (0.5.0, Imports `weightedsurv (>= 0.1.0)`). It calls `df_counting` (7 places), `plot_weighted_km` (5) and `plotKM.band_subgroups` (1).

| checklist item | yes/no | where |
|---|---|---|
| `tests/` present | no | — |
| `NEWS.md` present | no | — |
| `cran-comments.md` present | no | — |
| `set.seed()` inside package functions | yes | `km_wlr_calculations_helpers.R:483`, `df_counting_helpers.R:261`, `weightedCox.R:417` |
| `cat()` / `print()` not behind `verbose` | yes | `weighted_helpers.R:253` (unconditional `cat`). `check_results()` prints by default because `verbose = TRUE` is its default (`df_counting_helpers.R:106,134-143`). The others are behind `verbose`: `weightedCox.R:21-22,314`, `df_counting.R:416` |
| estimation functions that draw plots as a side effect | yes | `KM_diff()` draws `matplot` when `draws.band > 0` with default `show_resamples = TRUE` (`km_wlr_calculations_helpers.R:502-505`); `df_counting(check.seKM = TRUE)` (`df_counting.R:441-447`); `cumulative_rmst_bands(plot = TRUE)` by default (`rmst_calculations.R:121-150`) |
| every `par()` change restored with `on.exit` | yes | the only `par()` change is `df_counting.R:442-444`, restored by `on.exit`; `kmplotting_helpers.R:171` only reads `par("usr")` |
| `requireNamespace()` guards on packages in Imports | yes | `survival`: `df_counting.R:155,185`; `ggplot2`: `weighted_helpers.R:325` |
| functions used in formulas but not imported | yes | `strata()` at `df_counting.R:158,195` (R6). Masked by survival ≥ 3.8-1 |
| every Suggests package used, and used conditionally | no | `pkgdown` is used nowhere. `gt` is conditional in R/ (`baseline_table.R:66`). The vignette loads DiagrammeR, dplyr, data.table, gt, adjustedCurves and pammtools with unconditional `library()` (`vignettes/weightedsurv_examples.Rmd:71,152,270,625-626,719-720`) |
| examples over 5 s | no | max 0.46 s |
| `\dontrun` / `\donttest` | `\donttest` only | `man/check_results.Rd:51`, `man/create_baseline_table.Rd:51`; no `\dontrun` |
| URLs valid | yes | `urlchecker::url_check()`: "All URLs are correct!" (11 URLs) |

---

## F. Exports (`exports_inventory.csv`)

46 exports: **16 user-facing, 30 helpers.**

- User-facing: `KM_diff`, `df_counting`, `cox_rhogamma`, `cumulative_rmst_bands`, `plotKM.band_subgroups`, `plot_weighted_km`, `KM_plot_2sample_weighted_counting`, `plot_weight_schemes`, `create_baseline_table`, `wt.rg.S`, `wlr_dhat_estimates`, `check_results`, `plot_km`, `get_dfcounting`, `km_quantile`, `km_quantile_table`.
- Helpers used in a vignette: `KM_estimates`, `cox_score_rhogamma`, `wlr_cumulative`.
- Helpers not used outside `R/` (27): `N_rhogamma`, `add_legends`, `add_median_annotation`, `add_risk_table`, `calculate_risk_event_counts`, `check_km_curve`, `ci_cox`, `count_weighted`, `cox_rhogamma_resample`, `extract_and_calc_weights`, `extract_group_data`, `find_cox_root`, `format_pval`, `get_censoring_and_events`, `get_event_risk_matrices`, `get_riskpoints`, `get_validated_weights`, `get_weights`, `kmq_calculations`, `plot_km_confint_polygon`, `plot_km_curves_counting`, `resampling_survival`, `risk_weighted`, `safe_run`, `score_calculation`, `validate_input`, `validate_scheme_params`.
- Of these, 3 have no caller anywhere in `R/`: `N_rhogamma`, `count_weighted`, `risk_weighted`.
- The exported `check_km_curve` (`kmplotting_helpers.R:110`) is not the function `df_counting` calls; `df_counting` calls its own inner definition (S10).
- The reverse dependency forestsearch uses only `df_counting`, `plot_weighted_km` and `plotKM.band_subgroups`.

---

## G. Post-conditions (Gate G2)
