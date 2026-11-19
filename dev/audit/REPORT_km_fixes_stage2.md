# weightedsurv 0.2.0: Stage-2 fixes from the KM audit

**Top line.** Branch `km-fixes-0.2.0`, START_SHA `166919ec7128f165fb380847f549fa26ca1db985`; final SHA = the commit that adds this report (parent `c5fddce`). Gates G0, G1 and G2 passed (G2 output in section G). All nine authorised packages are DONE: B7, B6, V, B4, B2, B1, B3, B5, B10. None is NOT DONE. B8, B9 and M were held and not touched. M1b was run (kickoff: `M1b: GO`, report `dev/audit/REPORT_km_m1b.md`). Nothing was merged or pushed. Installs went only into temporary libraries under the session scratchpad.

Execution notes (findings, not package results):
- **Follow-up commits.** B4 has two follow-ups, `d165614` and `8696d62`:
  - `R CMD check` flagged `plot.new` as not imported in the new empty-panel branch.
  - The Stage-1 pressure test (T7) showed that an arm missing from a stratum got a stratum curve of 1 rather than `NA`.
  
  B2 has one follow-up, `c5fddce`, which removes a stray `tests/testthat/Rplots.pdf` that a draft test had left and that went in with the B2 commit.
- **Cumulative whitelists.** Each package's golden compare used the whitelists of the packages committed before it (for example B3 was compared with `B1,B3`), since `golden_before.rds` predates every package. `golden.R compare` was extended to accept a comma-separated list. The per-package outputs for the packages that moved fields are `golden_compare_B1.txt`, `_B3.txt` and `_B5.txt`.
- **Golden scope.** The golden snapshot leaves the `strata.name` call out for GBSG `er <= 0`: grade 1 has no treated patient there, so 0.1.0 stops with R1. Every other dataset × call ran on 0.1.0 (5,409 fields).
- **roxygen2 version.** It is 8.1.0 here, against `RoxygenNote: 7.3.3`. `devtools::document()` therefore reformatted the `importFrom` lines of `NAMESPACE` (the `export()` lines are unchanged) and replaced `RoxygenNote` with `Config/roxygen2/version: 8.1.0`. It also rewrote the family links in four Rd files and added an "Authors:" block to `weightedsurv-package.Rd`.
- **testthat floor.** `DESCRIPTION` declares `testthat (>= 3.0.0)` as the task specifies. The tests, however, use `local_mocked_bindings()` (testthat 3.1.7) and `expect_no_error()` / `expect_no_warning()` (3.1.5). They ran under testthat 3.3.2.
- **Stage-1 claim about `n`.** Stage 1 listed `n` in `calculate_risk_event_counts()` as unused. It was used once, by `rnorm(draws * n)`. V inlines it as `length(U)`.

---

## A. Packages

| ID | class | status | commit | files touched (R/, other) | golden fields that moved | arguments appended | tests added |
|---|---|---|---|---|---|---|---|
| B7 | B | DONE | `4cf5640` | `df_counting.R`; NAMESPACE, DESCRIPTION, 5 Rd (roxygen regeneration) | none | none | `test-b7-strata-import.R` (2 blocks), `helper-data.R` |
| B6 | B | DONE | `135b2ff` | `df_counting_helpers.R`, `km_wlr_calculations_helpers.R`, `weightedCox.R`, `weighted_helpers.R`, new `internal_helpers.R` | none | none | `test-b6-side-effects.R` (4) |
| V | moves code | DONE | `cb4c794` | `df_counting.R`, `df_counting_helpers.R`, `weighted_helpers.R` | none | none | `test-v-housekeeping.R` (2) |
| B4 | B | DONE | `32af3f2`, `d165614`, `8696d62` | `df_counting.R`, `df_counting_helpers.R`, `internal_helpers.R`, `km_wlr_calculations_helpers.R`, `kmplotting_helpers.R`, `weightedCox.R`; 3 Rd, NEWS | none | none | `test-b4-validation.R` (12) |
| B2 | B | DONE | `3b85230`, `c5fddce` | `km_wlr_calculations_helpers.R`, `kmplotting_helpers.R`; 2 Rd | none (differences below 1e-12) | none | `test-b2-curve-values.R` (5) |
| B1 | B | DONE | `f382b61` | `km_wlr_calculations_helpers.R`, `df_counting_helpers.R` (docs), `internal_helpers.R`; 3 Rd | every numeric output of the 8 `KM_diff(modify_tau = TRUE)` calls: 105 fields in 7 of 8 datasets (`sim_ties`: its 0.25 grid already holds every event time). `surv0`/`surv1` moved by up to 0.049 and now equal `survfit` at the grid. `cumulative_rmst_bands()` was not fed such a fit and did not move | none | `test-b1-arbitrary-points.R` (6) |
| B3 | B | DONE | `f1b7816` | `km_wlr_calculations_helpers.R`, `df_counting.R`, `internal_helpers.R`; 3 Rd | 28 fields. `quantile_results$quantile` moved by 0.5 in veteran (`default`, `MB`, `strata` calls; a plateau at S = 0.5). `km_quantile_arm1$qhat` moved by 0.5 in veteran. `quantile_check` of the `qprob = 0.25` call moved in all 8 datasets (convention now `probs = 0.75`). Every moved value equals `survfit`'s (contract check in `golden.R`) | none | `test-b3-quantiles.R` (6) |
| B5 | B | DONE | `06ec1dc` | `km_wlr_calculations_helpers.R`, `rmst_calculations.R` (input check only), `weighted_helpers.R`; 3 Rd | `wlr_dhat/sig2_dhat` and `cov_wlr_dhat` in 6 of 8 datasets. Old/new for `sig2_dhat`: veteran 0.970, gbsg_er_le0 0.952, sim_cont 0.980, sim_ties 0.919, sim_weighted 0.987, sim_strata 0.987, gbsg and gbsg_er_gt0 1.000. For `cov_wlr_dhat`: 0.985, 0.976, 0.990, 0.959, 0.993, 0.993, 1.000, 1.000. `cor_wlr_dhat` unchanged everywhere | none | `test-b5-inference.R` (6) |
| B10 | B | DONE | `060da51` | `kmplotting_helpers.R`, `internal_helpers.R`; 2 Rd | none | none | `test-b10-km-rendering.R` (4) |

Infrastructure: `b68df5f` adds `tests/testthat.R` and sets `DESCRIPTION` to `Version: 0.2.0`, with testthat in Suggests, `Config/testthat/edition: 3`, and `pkgdown` removed. It also adds `NEWS.md`. Final test suite: 9 files, 47 `test_that` blocks, 475 expectations, 1.7 s inside `R CMD check`. No test asserts a weighted-SE value.

## B. Finding register after the fixes

Before = Stage-1 status at `166919e` (identical to CRAN 0.1.0). After = final tree. Probe items come from `stage2/probe_after.txt`. S-items come from `stage2/static_leads_after.txt`; that script stops at its N2 block because it expects `df_counting()` to error there, and it no longer does. N-items come from `stage2/n_findings_after.txt`.

| ID | before | after | how verified |
|---|---|---|---|
| OK1 | PASS 2.2e-16 | PASS 2.2e-16 | probe |
| OK2 | PASS 5.6e-16 | PASS 5.6e-16 | probe |
| OK3 | PASS 0.923–1.005 | PASS 0.923–1.005 | probe |
| K1 | REPRODUCED, error 0.43, S down to −0.31 | NOT_REPRODUCED, error 0.0000, min S 0.12 (B1) | probe; T2.3 and T10 pass; B1 tests; golden B1 contract check |
| K1b | REPRODUCED 0.23 | NOT_REPRODUCED 0.0000 (B1) | probe; T2.7, T8 pass; B1 tests |
| K2 | REPRODUCED 23/30 | NOT_REPRODUCED 0/30 (B3) | probe; T4 0/7,401 failures; B3 tests |
| K2b | REPRODUCED 5.00 vs 5.50, silent | NOT_REPRODUCED 5.50 vs 5.50 (B3) | probe |
| K3 | REPRODUCED 0.18/0.64/1.60/6.31 | REPRODUCED, unchanged (M held) | probe; T2 K3 cells |
| K4 | REPRODUCED 2.06/1.77/1.52/1.25 | REPRODUCED, unchanged (M held) | probe; M1b |
| K5 | REPRODUCED 0.0297 vs 0.0416 | REPRODUCED, unchanged (M held) | probe |
| K6 | REPRODUCED 168/400; plot stops 38/60 | NOT_REPRODUCED 0/400; 0/60 (B2) | probe; T3 0 failures; T10 60/60 plots; B2 tests |
| K7 | REPRODUCED, S0 −166 to −38 | NOT_REPRODUCED: error "at_points must be non-decreasing" (B1) | probe; T8 |
| R1 | REPRODUCED (`sig2w_multiplier`) | NOT_REPRODUCED: "No observations in the control arm (treat == 0)." (B4) | probe; T8; B4 tests |
| R2 | REPRODUCED (subscript out of bounds) | NOT_REPRODUCED: result, 0 warnings (B4) | probe; T8 |
| R3a | REPRODUCED | NOT_REPRODUCED: "Column 'time' has 1 missing (NA) value(s)…" (B4) | probe; T8 |
| R3b | REPRODUCED | NOT_REPRODUCED, same message (B4) | probe; T8 |
| R3c | REPRODUCED (silent NA) | NOT_REPRODUCED: error naming column `w` (B4) | probe; T8 |
| R4a | REPRODUCED (silent) | NOT_REPRODUCED: "…(weights) has negative values" (B4) | probe; T8 |
| R4b | REPRODUCED (false discrepancy) | NOT_REPRODUCED: KM correct, 0 discrepancy warnings (B4) | probe; T6 0 failures; B4 tests |
| R5 | REPRODUCED for 1, 5, 19 | NOT_REPRODUCED (B4) | probe; T5; B4 tests |
| R6 | NOT_REPRODUCED (survival 3.8.11) | NOT_REPRODUCED; `strata` now imported (B7) | probe; T7 with/without survival attached; B7 test |
| R7 | REPRODUCED | NOT_REPRODUCED: 1/2 accepted by all three functions (B4) | probe; T8; B4 tests |
| I1 | REPRODUCED 2.9796 at α 0.05 and 0.20 | NOT_REPRODUCED 2.9796 vs 2.3417 (B5) | probe; T5; B5 tests |
| I2a | REPRODUCED ratio 0.9467 | NOT_REPRODUCED ratio 1.0000 (B5) | probe; B5 tests; golden ratios |
| I2b | REPRODUCED length 0 | NOT_REPRODUCED length 1; before first point: 0 (B5) | probe; T8; B5 tests |
| I3 | REPRODUCED 4.2836 vs 4.2688 | REPRODUCED, unchanged (B8 held) | probe; G2 check 2a |
| I4 | REPRODUCED 0 vs 2 warnings | NOT_REPRODUCED 0 vs 0 (B3) | probe; B3 tests |
| I5 | REPRODUCED 15.54 vs 4.29 | NOT_REPRODUCED: `quantile_check` 15.54 = `quantile_results` (B3) | probe; B3 tests |
| P2 | REPRODUCED | NOT_REPRODUCED: 0.287578 / 0.389071 (B6) | probe; T5; B6 tests |
| P5 | INFO, 0.3 s at n = 2,000 | INFO, 0.3 s; T9: `df_counting` 3.61 s at n = 8,000 (3.71 before) | probe; T9 |
| S1 | CONFIRMED | Fixed. At `mb_tstar == max`, weights 1/max(S(t−), S(last)); beyond: finite (B5) | static_leads_after; B5 tests |
| S2 | CONFIRMED (cryptic `data.frame` error) | Fixed: "…missing required elements: survG (scheme = 'schemper' needs the censoring curve…)"; no partial matching (B5) | static_leads_after; B5 tests |
| S3 | CONFIRMED, no numeric effect | Unchanged default grid. B1 makes any grid exact, so raw versus fixed times cannot matter | static_leads_after; B1 tests |
| S4 | CONFIRMED ("dim(X)…") | Fixed: "needs a KM_diff() fit made with draws.band > 0". Levels and integration unchanged (B8 held) | T10; B5 tests |
| S5 | CONFIRMED (NaN, upper > 1) | Fixed for `km_quantile` (NA where S = 0, upper ≤ 1), `plot_km_confint_polygon` (capped, closes) and the drawn band (B2, B10) | T8/T10; Step P band max 1.000; B2, B10 tests |
| S6 | CONFIRMED (0.040) | Fixed through B1: `fit_itt` vs `survfit` within 1e-6 | T10 cell passes |
| S7 | CONFIRMED (`cat`) | Fixed: stdout empty; the names are in the error (B6) | static_leads_after; B6 tests |
| S8 | CONFIRMED | Fixed: "Missing required column(s) in df: no_such_column"; event 3 is an informative error (B4) | static_leads_after; B4 tests |
| S9 | CONFIRMED (recorded only) | Unchanged by design: one `seedstart` for every arm and stratum | static_leads_after |
| S10 | CONFIRMED | Fixed. Inner function renamed `check_km_fit`; duplicate `ans$strata` and unused `dJ` removed; `n` inlined (V); dead branch removed (B5); `pkgdown` removed from Suggests | V tests; diff |
| N1 | cryptic ("first argument must be a vector") | Result plus one warning listing the NA quantities (B4) | n_findings_after; T8 |
| N2 | cryptic (TRUE/FALSE) | Result: `z.score = NA`, one warning (B4) | n_findings_after; T1 0 errors |
| N3 | cryptic (dgesv singular) | Result: log-rank `chisq = NA`, one warning with the `survdiff` message (B4) | n_findings_after; B4 tests |
| N4 | 5/3,000 upper limits differ | 0/3,000 (B3) | T4 |
| N5 | silent divergent root | Same value (bhat −1243.65) plus warning "…not finite / not identified" (B4) | n_findings_after; B4 tests |
| N6 | `coxph` error on zero weights | Result: Cox on positive-weight rows, HR 1.104 (B4) | n_findings_after; B4 tests |
| N7 | cryptic ("'to' must be a finite number") | "max_tau is undefined… Supply Maxtau…"; with `Maxtau` a result (B4) | n_findings_after; T8; B4 tests |
| N8 | static: log-rank is unweighted `survdiff` | Unchanged (not in any authorised package) | n_findings_after |
| N9 | silent wrong result in all three | Error "…has negative values; times must be >= 0." in all three (B4) | n_findings_after; T8 |

## C. Pressure test before / after

`dev/audit/stage2/km_pressure_test_s2.R` is a copy of the Stage-1 script. It differs only in the `--out=` output directory (default `dev/audit/stage2/`) and an `--only=` family filter used during development; the tests themselves are unchanged. Full run: 115 s, exit 0, written to `stage2/km_pressure_results.csv` and `km_scaling.csv`. "Before" is the Stage-1 `dev/audit/km_pressure_results.csv`.

| family | cells | pass before → after | fail before → after | report-only before → after | checks (failed) before → after |
|---|---|---|---|---|---|
| REF | 1 | 1 → 1 | 0 → 0 | 0 → 0 | 12 (0) → 12 (0) |
| T1 | 141 | 101 → 111 | 10 → 0 | 30 → 30 | 41178 (130) → 41438 (0) |
| T2 | 18 | 11 → 14 | 6 → 3 | 1 → 1 | 513 (93) → 513 (3) |
| T3 | 288 | 204 → 288 | 84 → 0 | 0 → 0 | 84883 (5205) → 86376 (0) |
| T4 | 22 | 7 → 21 | 14 → 0 | 1 → 1 | 7401 (281) → 7401 (0) |
| T5 | 11 | 4 → 9 | 5 → 0 | 2 → 2 | 401 (9) → 401 (0) |
| T6 | 24 | 6 → 24 | 18 → 0 | 0 → 0 | 888 (231) → 888 (0) |
| T7 | 18 | 12 → 18 | 6 → 0 | 0 → 0 | 138 (6) → 192 (0) |
| T8 | 83 | 47 → 73 | 26 → 0 | 10 → 10 | 73 (26) → 73 (0) |
| T9 | 1 | 1 → 1 | 0 → 0 | 0 → 0 | 12 (0) → 12 (0) |
| T10 | 5 | 0 → 2 | 5 → 3 | 0 → 0 | 64 (46) → 64 (3) |

The check counts in T1, T3 and T7 rise because datasets that used to stop `df_counting()` (N2/N3, R1 in a stratum) now return results that are checked. T8 now has no cryptic error and no silent wrong result. All 73 classified cells are correct results or informative errors. The `Inf`-time cells of `KM_diff()` and `plotKM.band_subgroups()` went from correct result to informative error, per the B4 contract.

Remaining failures (6 cells):
- **T2, three cells:** "weights w → c·w, c = 0.01 / 10 / 1/Σw: SE ratio (K3)", worst 50.5 / 1.41 / 168. Attributed to **M** (held).
- **T10, "KM_plot_2sample_weighted_counting(conf.int = TRUE) when a curve reaches 0 (S5)":** the cell recomputes the log-transform limits itself from `surv0`/`sig2_surv0` with the old formula, and counts its own NaN and values above 1 (3 NaN, max 1.02). It does not read what is drawn. The plot drew without error (`worst = drew`, `p$error` NULL). The band actually drawn comes from `km_step_vertices()`, which caps at 1 and has no NaN (Step P: band max 1.000 in all 74 arms; B2/B10 tests). The NaN in the cell is NOT a package failure.
- **T10, "plotKM.band_subgroups() with a subgroup whose treatment arm is empty (R1)":** this now gives the informative error "Subgroup 'sg == 'B'': No observations in the treatment arm (treat == 1)." The B4 contract requires an error that names the subgroup; the cell counts any error as a failure.
- **T10, "cumulative_rmst_bands() given KM_diff(draws.band = 0) (S4)":** this now gives the informative error "cumulative_rmst_bands() needs a KM_diff() fit made with draws.band > 0 (fit$dhat_star is missing)." The B5/S4 contract requires it; the cell counts any error as a failure.

Report-only rows: the T1 weighted-SE ratios are unchanged from Stage 1, because the variance formula is untouched (M). For example, normalised weights with no ties give a median per-arm ratio to `survfit` robust of 1.84, range 0.205–1,010. T9 at n = 8,000: `df_counting` 3.61 s / 1.6 GB, `KM_diff(draws.band = 1000)` 3.16 s, `cox_rhogamma` 11.44 s / 2.0 GB.

## D. Plot comparison (GBSG, `km_plot_compare_before.csv` → `km_plot_compare_after.csv`)

Columns:
- "Marks before it": censored observations before the arm's first event.
- "First drawn vertex": the first point the package puts on the device.
- "Max vertex diff": the maximum absolute difference between the drawn step function and `survfit`'s (`c(0, time)`, `c(1, surv)`, and `lower` / `upper` with `conf.type = "log"`). It is evaluated at every jump and every midpoint in [0, last time], wherever the package draws something.

In "before", the band is evaluated as drawn: a polygon, linear between vertices. All 37 subgroups ran in both modes (`status = ok`).

| subgroup | arm | n | events | first event | marks before it | first drawn vertex before → after | band max before → after | max vertex diff, curve: before → after | max vertex diff, band: before → after |
|---|---|---|---|---|---|---|---|---|---|
| er <= 0 | 0 | 56 | 29 | 6.41 | 0 | (6.41, 0.982) → (0.00, 1.000) | 1.017 → 1.000 | 0 → 0 | 2.6e-02 → 1.1e-16 |
| er <= 0 | 1 | 26 | 16 | 5.55 | 2 | (5.55, 0.958) → (0.00, 1.000) | 1.042 → 1.000 | 5.6e-17 → 5.6e-17 | 5.4e-02 → 1.1e-16 |
| er > 0 | 0 | 384 | 176 | 2.37 | 10 | (2.37, 0.997) → (0.00, 1.000) | 1.003 → 1.000 | 1.1e-16 → 1.1e-16 | 2.3e-01 → 3.3e-16 |
| er > 0 | 1 | 220 | 78 | 5.82 | 4 | (5.82, 0.991) → (0.00, 1.000) | 1.004 → 1.000 | 1.1e-16 → 1.1e-16 | 4.7e-02 → 2.2e-16 |
| er <= 5 | 0 | 102 | 59 | 3.71 | 2 | (3.71, 0.990) → (0.00, 1.000) | 1.010 → 1.000 | 1.1e-16 → 1.1e-16 | 3.0e-02 → 2.2e-16 |
| er <= 5 | 1 | 52 | 28 | 5.55 | 2 | (5.55, 0.980) → (0.00, 1.000) | 1.020 → 1.000 | 0 → 0 | 4.7e-02 → 1.1e-16 |
| er > 5 | 0 | 338 | 146 | 2.37 | 8 | (2.37, 0.997) → (0.00, 1.000) | 1.003 → 1.000 | 1.1e-16 → 1.1e-16 | 2.5e-01 → 3.3e-16 |
| er > 5 | 1 | 194 | 66 | 5.82 | 4 | (5.82, 0.995) → (0.00, 1.000) | 1.005 → 1.000 | 1.1e-16 → 1.1e-16 | 4.9e-02 → 2.2e-16 |
| er <= 20 | 0 | 181 | 97 | 3.71 | 3 | (3.71, 0.994) → (0.00, 1.000) | 1.005 → 1.000 | 1.1e-16 → 1.1e-16 | 1.7e-01 → 2.2e-16 |
| er <= 20 | 1 | 91 | 40 | 5.55 | 4 | (5.55, 0.989) → (0.00, 1.000) | 1.011 → 1.000 | 0 → 0 | 3.0e-02 → 1.1e-16 |
| er > 20 | 0 | 259 | 108 | 2.37 | 7 | (2.37, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 2.8e-16 → 2.8e-16 | 7.6e-02 → 4.4e-16 |
| er > 20 | 1 | 155 | 54 | 5.82 | 2 | (5.82, 0.993) → (0.00, 1.000) | 1.006 → 1.000 | 0 → 0 | 4.5e-02 → 1.1e-16 |
| er <= 50 | 0 | 264 | 129 | 2.37 | 8 | (2.37, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 1.1e-16 → 1.1e-16 | 1.8e-01 → 2.2e-16 |
| er <= 50 | 1 | 128 | 51 | 5.55 | 4 | (5.55, 0.992) → (0.00, 1.000) | 1.008 → 1.000 | 0 → 0 | 2.3e-02 → 1.1e-16 |
| er > 50 | 0 | 176 | 76 | 6.74 | 2 | (6.74, 0.994) → (0.00, 1.000) | 1.006 → 1.000 | 2.2e-16 → 2.2e-16 | 8.3e-02 → 3.3e-16 |
| er > 50 | 1 | 118 | 43 | 5.82 | 2 | (5.82, 0.991) → (0.00, 1.000) | 1.008 → 1.000 | 1.7e-16 → 1.7e-16 | 6.5e-02 → 2.2e-16 |
| er > 100 | 0 | 102 | 40 | 9.23 | 1 | (9.23, 0.990) → (0.00, 1.000) | 1.010 → 1.000 | 1.1e-16 → 1.1e-16 | 9.3e-02 → 2.2e-16 |
| er > 100 | 1 | 84 | 25 | 5.82 | 1 | (5.82, 0.988) → (0.00, 1.000) | 1.012 → 1.000 | 1.1e-16 → 1.1e-16 | 1.1e-01 → 2.2e-16 |
| er <= 0 & meno == 0 | 0 | 33 | 16 | 9.36 | 0 | (9.36, 0.970) → (0.00, 1.000) | 1.030 → 1.000 | 2.2e-16 → 2.2e-16 | 4.9e-02 → 3.3e-16 |
| er <= 0 & meno == 0 | 1 | 7 | 4 | 5.55 | 1 | (5.55, 0.833) → (0.00, 1.000) | 1.221 → 1.000 | 1.1e-16 → 1.1e-16 | 2.2e-01 → 5.6e-17 |
| er <= 0 & meno == 1 | 0 | 23 | 13 | 6.41 | 0 | (6.41, 0.957) → (0.00, 1.000) | 1.044 → 1.000 | 1.1e-16 → 1.1e-16 | 5.5e-02 → 2.2e-16 |
| er <= 0 & meno == 1 | 1 | 19 | 12 | 5.91 | 1 | (5.91, 0.944) → (0.00, 1.000) | 1.056 → 1.000 | 5.6e-17 → 5.6e-17 | 5.6e-02 → 2.2e-16 |
| er <= 0 & grade == 3 | 0 | 29 | 15 | 6.41 | 0 | (6.41, 0.966) → (0.00, 1.000) | 1.034 → 1.000 | 1.7e-16 → 1.7e-16 | 5.7e-02 → 3.3e-16 |
| er <= 0 & grade == 3 | 1 | 8 | 6 | 5.91 | 0 | (5.91, 0.875) → (0.00, 1.000) | 1.137 → 1.000 | 5.6e-17 → 5.6e-17 | 1.4e-01 → 3.3e-16 |
| er <= 0 & nodes > 3 | 0 | 28 | 19 | 6.41 | 0 | (6.41, 0.964) → (0.00, 1.000) | 1.036 → 1.000 | 1.1e-16 → 1.1e-16 | 4.5e-02 → 2.2e-16 |
| er <= 0 & nodes > 3 | 1 | 16 | 12 | 5.55 | 1 | (5.55, 0.933) → (0.00, 1.000) | 1.069 → 1.000 | 5.6e-17 → 5.6e-17 | 6.9e-02 → 2.2e-16 |
| er <= 0 & pgr <= 0 | 0 | 38 | 21 | 9.63 | 0 | (9.63, 0.974) → (0.00, 1.000) | 1.026 → 1.000 | 1.1e-16 → 1.1e-16 | 3.5e-02 → 2.2e-16 |
| er <= 0 & pgr <= 0 | 1 | 22 | 13 | 5.91 | 2 | (5.91, 0.950) → (0.00, 1.000) | 1.050 → 1.000 | 1.1e-16 → 1.1e-16 | 6.8e-02 → 2.2e-16 |
| er <= 0 & age <= 50 | 0 | 33 | 16 | 9.36 | 0 | (9.36, 0.970) → (0.00, 1.000) | 1.030 → 1.000 | 0 → 0 | 5.4e-02 → 1.1e-16 |
| er <= 0 & age <= 50 | 1 | 8 | 5 | 5.55 | 1 | (5.55, 0.857) → (0.00, 1.000) | 1.168 → 1.000 | 1.1e-16 → 1.1e-16 | 2.2e-01 → 1.1e-16 |
| er > 0 & pgr > 0 | 0 | 364 | 163 | 2.37 | 10 | (2.37, 0.997) → (0.00, 1.000) | 1.003 → 1.000 | 1.7e-16 → 1.7e-16 | 2.4e-01 → 2.2e-16 |
| er > 0 & pgr > 0 | 1 | 212 | 72 | 5.82 | 3 | (5.82, 0.990) → (0.00, 1.000) | 1.004 → 1.000 | 1.1e-16 → 1.1e-16 | 5.2e-02 → 2.2e-16 |
| er > 0 & grade == 3 | 0 | 82 | 41 | 3.22 | 2 | (3.22, 0.988) → (0.00, 1.000) | 1.012 → 1.000 | 0 → 0 | 6.2e-02 → 1.1e-16 |
| er > 0 & grade == 3 | 1 | 42 | 17 | 5.82 | 1 | (5.82, 0.951) → (0.00, 1.000) | 1.019 → 1.000 | 1.1e-16 → 1.1e-16 | 1.4e-01 → 2.2e-16 |
| er > 0 & age <= 40 | 0 | 50 | 31 | 5.68 | 1 | (5.68, 0.980) → (0.00, 1.000) | 1.020 → 1.000 | 5.6e-17 → 5.6e-17 | 2.8e-02 → 1.1e-16 |
| er > 0 & age <= 40 | 1 | 12 | 4 | 9.40 | 0 | (9.40, 0.917) → (0.00, 1.000) | 1.087 → 1.000 | 0 → 0 | 1.5e-01 → 0 |
| pgr <= 0 | 0 | 58 | 34 | 8.18 | 1 | (8.18, 0.982) → (0.00, 1.000) | 1.017 → 1.000 | 1.1e-16 → 1.1e-16 | 2.6e-02 → 2.8e-16 |
| pgr <= 0 | 1 | 30 | 19 | 5.91 | 3 | (5.91, 0.963) → (0.00, 1.000) | 1.037 → 1.000 | 2.2e-16 → 2.2e-16 | 5.8e-02 → 3.3e-16 |
| pgr > 0 | 0 | 382 | 171 | 2.37 | 10 | (2.37, 0.997) → (0.00, 1.000) | 1.003 → 1.000 | 1.1e-16 → 1.1e-16 | 1.6e-01 → 2.2e-16 |
| pgr > 0 | 1 | 216 | 75 | 5.55 | 3 | (5.55, 0.995) → (0.00, 1.000) | 1.005 → 1.000 | 0 → 0 | 5.2e-02 → 1.1e-16 |
| pgr <= 20 | 0 | 178 | 99 | 2.37 | 4 | (2.37, 0.994) → (0.00, 1.000) | 1.006 → 1.000 | 1.1e-16 → 1.1e-16 | 5.5e-02 → 2.2e-16 |
| pgr <= 20 | 1 | 99 | 55 | 5.55 | 5 | (5.55, 0.989) → (0.00, 1.000) | 1.010 → 1.000 | 1.1e-16 → 1.1e-16 | 6.6e-02 → 2.2e-16 |
| pgr > 20 | 0 | 262 | 106 | 3.22 | 6 | (3.22, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 1.1e-16 → 1.1e-16 | 8.1e-02 → 2.2e-16 |
| pgr > 20 | 1 | 147 | 39 | 5.82 | 1 | (5.82, 0.993) → (0.00, 1.000) | 1.007 → 1.000 | 0 → 0 | 3.5e-02 → 1.1e-16 |
| meno == 0 | 0 | 231 | 97 | 3.94 | 7 | (3.94, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 5.6e-17 → 5.6e-17 | 3.4e-02 → 1.1e-16 |
| meno == 0 | 1 | 59 | 22 | 5.55 | 1 | (5.55, 0.983) → (0.00, 1.000) | 1.017 → 1.000 | 0 → 0 | 4.6e-02 → 1.1e-16 |
| meno == 1 | 0 | 209 | 108 | 2.37 | 3 | (2.37, 0.995) → (0.00, 1.000) | 1.005 → 1.000 | 1.1e-16 → 1.1e-16 | 1.7e-01 → 2.2e-16 |
| meno == 1 | 1 | 187 | 72 | 5.82 | 5 | (5.82, 0.995) → (0.00, 1.000) | 1.005 → 1.000 | 1.1e-16 → 1.1e-16 | 5.0e-02 → 2.2e-16 |
| grade == 1 | 0 | 48 | 12 | 15.64 | 3 | (15.64, 0.978) → (0.00, 1.000) | 1.022 → 1.000 | 1.1e-16 → 1.1e-16 | 1.5e-01 → 1.1e-16 |
| grade == 1 | 1 | 33 | 6 | 18.37 | 5 | (18.37, 0.964) → (0.00, 1.000) | 1.038 → 1.000 | 1.1e-16 → 1.1e-16 | 2.5e-01 → 1.1e-16 |
| grade == 2 | 0 | 281 | 137 | 2.37 | 7 | (2.37, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 2.2e-16 → 2.2e-16 | 1.1e-01 → 3.3e-16 |
| grade == 2 | 1 | 163 | 65 | 5.55 | 4 | (5.55, 0.994) → (0.00, 1.000) | 1.006 → 1.000 | 0 → 0 | 5.6e-02 → 1.1e-16 |
| grade == 3 | 0 | 111 | 56 | 3.22 | 2 | (3.22, 0.991) → (0.00, 1.000) | 1.009 → 1.000 | 5.6e-17 → 5.6e-17 | 6.5e-02 → 1.1e-16 |
| grade == 3 | 1 | 50 | 23 | 5.82 | 1 | (5.82, 0.959) → (0.00, 1.000) | 1.016 → 1.000 | 0 → 0 | 1.2e-01 → 5.6e-17 |
| nodes <= 3 | 0 | 248 | 88 | 6.28 | 6 | (6.28, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 2.2e-16 → 2.2e-16 | 2.2e-02 → 3.3e-16 |
| nodes <= 3 | 1 | 128 | 31 | 12.12 | 8 | (12.12, 0.992) → (0.00, 1.000) | 1.008 → 1.000 | 1.1e-16 → 1.1e-16 | 8.5e-02 → 2.2e-16 |
| nodes > 3 | 0 | 192 | 117 | 2.37 | 6 | (2.37, 0.995) → (0.00, 1.000) | 1.005 → 1.000 | 1.1e-16 → 1.1e-16 | 1.1e-01 → 2.2e-16 |
| nodes > 3 | 1 | 118 | 63 | 5.55 | 3 | (5.55, 0.991) → (0.00, 1.000) | 1.008 → 1.000 | 0 → 0 | 3.8e-02 → 1.1e-16 |
| nodes > 9 | 0 | 67 | 48 | 3.22 | 1 | (3.22, 0.985) → (0.00, 1.000) | 1.015 → 1.000 | 1.1e-16 → 1.1e-16 | 3.1e-02 → 2.2e-16 |
| nodes > 9 | 1 | 36 | 25 | 5.82 | 1 | (5.82, 0.971) → (0.00, 1.000) | 1.028 → 1.000 | 0 → 0 | 9.6e-02 → 1.1e-16 |
| size <= 20 | 0 | 113 | 47 | 2.37 | 4 | (2.37, 0.991) → (0.00, 1.000) | 1.009 → 1.000 | 1.1e-16 → 1.1e-16 | 3.2e-02 → 2.2e-16 |
| size <= 20 | 1 | 67 | 18 | 5.55 | 4 | (5.55, 0.984) → (0.00, 1.000) | 1.015 → 1.000 | 1.1e-16 → 1.1e-16 | 5.4e-02 → 2.2e-16 |
| size > 20 | 0 | 327 | 158 | 3.22 | 6 | (3.22, 0.997) → (0.00, 1.000) | 1.003 → 1.000 | 1.1e-16 → 1.1e-16 | 1.2e-01 → 2.2e-16 |
| size > 20 | 1 | 179 | 76 | 5.82 | 2 | (5.82, 0.989) → (0.00, 1.000) | 1.004 → 1.000 | 0 → 0 | 5.2e-02 → 1.1e-16 |
| size > 50 | 0 | 39 | 23 | 3.71 | 1 | (3.71, 0.974) → (0.00, 1.000) | 1.026 → 1.000 | 0 → 0 | 5.6e-02 → 1.1e-16 |
| size > 50 | 1 | 14 | 7 | 8.94 | 1 | (8.94, 0.923) → (0.00, 1.000) | 1.161 → 1.000 | 1.1e-16 → 1.1e-16 | 2.6e-01 → 2.2e-16 |
| age <= 40 | 0 | 58 | 35 | 5.68 | 1 | (5.68, 0.982) → (0.00, 1.000) | 1.017 → 1.000 | 1.1e-16 → 1.1e-16 | 2.5e-02 → 1.1e-16 |
| age <= 40 | 1 | 15 | 5 | 9.40 | 1 | (9.40, 0.929) → (0.00, 1.000) | 1.074 → 1.000 | 1.1e-16 → 1.1e-16 | 1.1e-01 → 2.2e-16 |
| age > 40 & age <= 60 | 0 | 290 | 121 | 2.37 | 7 | (2.37, 0.996) → (0.00, 1.000) | 1.003 → 1.000 | 2.2e-16 → 2.2e-16 | 2.4e-02 → 3.3e-16 |
| age > 40 & age <= 60 | 1 | 135 | 54 | 5.55 | 3 | (5.55, 0.992) → (0.00, 1.000) | 1.007 → 1.000 | 1.7e-16 → 1.7e-16 | 2.5e-02 → 2.2e-16 |
| age > 60 | 0 | 92 | 49 | 6.41 | 3 | (6.41, 0.989) → (0.00, 1.000) | 1.011 → 1.000 | 5.6e-17 → 5.6e-17 | 1.0e-01 → 1.1e-16 |
| age > 60 | 1 | 96 | 35 | 5.82 | 2 | (5.82, 0.989) → (0.00, 1.000) | 1.010 → 1.000 | 1.1e-16 → 1.1e-16 | 1.2e-01 → 2.2e-16 |
| (all patients) | 0 | 440 | 205 | 2.37 | 10 | (2.37, 0.998) → (0.00, 1.000) | 1.002 → 1.000 | 1.1e-16 → 1.1e-16 | 1.6e-01 → 2.2e-16 |
| (all patients) | 1 | 246 | 94 | 5.55 | 6 | (5.55, 0.996) → (0.00, 1.000) | 1.004 → 1.000 | 1.7e-16 → 1.7e-16 | 4.6e-02 → 2.2e-16 |

Reading:
1. **Before**, in all 74 arms, the curve starts at the arm's first event time and not at (0, 1). The first drawn x ranges from 2.37 to 18.37 months (grade 1, treated: (18.37, 0.964)). In 65 of 74 arms, 236 censoring marks sit at S = 1 with no line under them. **After**, every curve starts at (0, 1.000).
2. Where both draw, the curve values are identical to `survfit` before and after (max 2.8e-16). The drawing differences are the missing first segment and the band, not the estimates.
3. **Before**, the band exceeds 1 in all 74 arms (max 1.221, `er <= 0 & meno == 0`, treated). Because the vertices are joined by straight lines, it differs from `survfit`'s step limits by 0.022 to 0.257 (worst: `size > 50`, treated). **After**, the band maximum is 1.000 everywhere and its vertices equal `survfit`'s within 4.4e-16.
4. **After**, curve and band are step functions from 0 to the arm's last observed time. The censoring marks, risk table, axes and legends are unchanged (B10 test; pages 1–37 of `km_plots_after.pdf`).
5. The difference-curve sheets (`plotKM.band_subgroups(draws.band = 0)`, 10 sheets of 4 subgroups) agree with `S1 − S0` from `survfit` at `xpoints` in all 37 subgroups, before and after (max 2.8e-16). B10 did not touch this function, and nothing in those sheets departs from `survfit`.

## E. CRAN check and forestsearch

`R CMD check --as-cran --timings weightedsurv_0.2.0.tar.gz`. The tarball came from `R CMD build` of a clean `git clone` of `c5fddce`; the R code there is identical to the final tree. The check ran with `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`, `DISPLAY` unset and RStudio's pandoc in `RSTUDIO_PANDOC`. Log: `stage2/00check_after.log`.

```
* checking CRAN incoming feasibility ... [3s/21s] OK
* checking for future file timestamps ... NOTE
unable to verify current time
* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
Status: 3 NOTEs
```

- **NOTEs.** 0 ERROR, 0 WARNING. Two NOTEs are the local-tooling notes from Stage 1 (pandoc not on `PATH`, old HTML Tidy). The "future file timestamps" NOTE is also local: the time servers R queries (worldtimeapi.org, worldclockapi.com) did not answer from this machine. The Stage-1 WARNING ("Insufficient package version") is gone because the version is now 0.2.0.
- **NOTEs in earlier runs.** A first check on a build from the working directory showed two more NOTEs:
  - a hidden `.claude` directory, a local untracked tool directory that a clean clone does not contain;
  - `plot.new` not imported. This was a real defect, fixed in `d165614`.
- **Tests and timings.** `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 475 ]`, 1.7 s elapsed. The slowest example is `KM_diff` at 0.48 s.

forestsearch 0.5.0, the CRAN source:
- **How it was run.** `R CMD check --no-manual --ignore-vignettes forestsearch_0.5.0.tar.gz`, with `R_LIBS` set to a temporary library holding the new weightedsurv 0.2.0 (built from `8696d62`, whose R code equals the final tree). All forestsearch dependencies were already installed in the user library, so nothing else was installed. Runtime 46 s, well inside the 30-minute box. Logs: `stage2/forestsearch_00check.log`, `stage2/forestsearch_testthat.Rout`.
- **Result.** `Status: OK`. Tests: `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]`. No new failure, so the comparison run with CRAN 0.1.0 was not needed.
- **Coverage.** No forestsearch test file or example calls `df_counting`, `plot_weighted_km` or `plotKM.band_subgroups`, so this check does not exercise forestsearch's use of weightedsurv.

## F. Policy choices made inside B4 that a user could notice

1. Non-finite weights (`Inf`) are an error, as negative weights are.
2. Event 1/2 coding is recognised exactly as `survival::Surv()` does: all values in {1, 2} with a maximum of 2. A column of all 1s is read as 0/1 (all events), as `Surv()` reads it. `df_counting()` and `plotKM.band_subgroups()` recode the column to 0/1 before Cox, log-rank and subgroup subsetting, so a subgroup with only 1s cannot be misread.
3. An arm whose observations all have weight 0 is treated as empty: "No observations with positive weight in the … arm".
4. `Inf` times are now an error in `KM_diff()` and `plotKM.band_subgroups()`, which used to return a correct result for them.
5. In `df_counting()`, when an arm has no events among the positive-weight rows, `coxph()` is not called. `cox_results` keeps its six element names with `cox_fit = NA`, `hr`, `hr_ci`, `pval` and `score` `NA`, and `cox_text = "HR = NA"`. A `coxph()` error is caught the same way.
6. In `df_counting()`, when there are no events at all, `survdiff()` is not called, so `chisq` and `pval` are `NA` and `logrank_text` is "Logrank p = NA". A `survdiff()` error (for example N3) is caught the same way.
7. The weighted log-rank `z.score` is `NA`, not `NaN`, when its variance is 0. `zlogrank_text` then shows `p = NA`.
8. All reasons are collected into one warning: "df_counting(): some quantities could not be computed and are NA: …".
9. In stratified `df_counting()`, an arm absent from a stratum gets `NA` for its `*_stratum` curve and variance. Its counts are 0, and the reason goes into the same warning.
10. In `df_counting()`, the survfit-based self-checks (KM curve, quantiles) use the positive-weight rows. The KM curve is compared at `survfit`'s event times.
11. `plotKM.band_subgroups()` raises the `max_tau` error when any arm has no events and `Maxtau` is `NULL`. With `Maxtau` supplied it runs.
12. `cox_rhogamma()` also returns its documented `status = "fail"` list, with a warning, when root finding returns `NA`. It used to stop with "$ operator is invalid for atomic vectors".
13. The `KM_diff(show_resamples = TRUE)` axis label states how many resamples are drawn ("1st 5" when `draws.band = 5`).
14. A missing `weight.name` column is now reported up front in `KM_diff()`.

## G. Gate G2 output

Run from the repo root with `START_SHA=166919ec7128f165fb380847f549fa26ca1db985` on HEAD `c5fddce`, before the commit that adds this report, `probe_after.txt` and `tests_after.txt`. Check 9 runs after that commit; its output is printed at the end of the session.

```
OK 1a
OK 1b
OK 2a
OK 2b
OK 3a
OK 3b
API compare: 46 exports before, 46 now
API: unchanged
OK 4a
final whitelist = union over: B7, B6, V, B4, B2, B1, B3, B5, B10
golden compare, whitelist 'final': 5409 fields in golden, 5409 in new run, 0 added
WHITELISTED fields that moved (145):
  veteran/df_counting_default/quantile_results/quantile: max abs diff 0.5
  veteran/df_counting_MB/quantile_results/quantile: max abs diff 0.5
  veteran/df_counting_q25/quantile_check/time: max abs diff 135
  veteran/df_counting_q25/quantile_check/lower: max abs diff 123
  veteran/df_counting_q25/quantile_check/upper: max abs diff 297
  veteran/KM_diff_band_modtau/surv0: max abs diff 0.0177
  veteran/KM_diff_band_modtau/sig2_surv0: max abs diff 0.000306
  veteran/KM_diff_band_modtau/surv1: max abs diff 0.0192
  veteran/KM_diff_band_modtau/sig2_surv1: max abs diff 0.000292
  veteran/KM_diff_band_modtau/dhat: max abs diff 0.00806
  veteran/KM_diff_band_modtau/sig2_dhat: max abs diff 0.000597
  veteran/KM_diff_band_modtau/lower: max abs diff 0.0187
  veteran/KM_diff_band_modtau/upper: max abs diff 0.0203
  veteran/KM_diff_band_modtau/dhat_star: max abs diff 0.0494
  veteran/KM_diff_band_modtau/Zdhat_star: max abs diff 1.22
  veteran/KM_diff_band_modtau/surv0_star: max abs diff 0.0426
  veteran/KM_diff_band_modtau/surv1_star: max abs diff 0.015
  veteran/KM_diff_band_modtau/c_alpha_band: max abs diff 0.00159
  veteran/KM_diff_band_modtau/sb_lower: max abs diff 0.0289
  veteran/KM_diff_band_modtau/sb_upper: max abs diff 0.0305
  veteran/wlr_dhat/cov_wlr_dhat: max abs diff 0.0051
  veteran/wlr_dhat/sig2_dhat: max abs diff 0.000213
  veteran/km_quantile_arm1/qhat: max abs diff 0.5
  veteran/df_counting_strata/quantile_results/quantile: max abs diff 0.5
  gbsg/df_counting_q25/quantile_check/time: NA pattern differs
  gbsg/df_counting_q25/quantile_check/lower: NA pattern differs
  gbsg/df_counting_q25/quantile_check/upper: NA pattern differs
  gbsg/KM_diff_band_modtau/surv0: max abs diff 0.00503
  gbsg/KM_diff_band_modtau/sig2_surv0: max abs diff 1.25e-05
  gbsg/KM_diff_band_modtau/surv1: max abs diff 0.00617
  gbsg/KM_diff_band_modtau/sig2_surv1: max abs diff 6.44e-05
  gbsg/KM_diff_band_modtau/dhat: max abs diff 0.00151
  gbsg/KM_diff_band_modtau/sig2_dhat: max abs diff 7.68e-05
  gbsg/KM_diff_band_modtau/lower: max abs diff 0.00253
  gbsg/KM_diff_band_modtau/upper: max abs diff 0.00206
  gbsg/KM_diff_band_modtau/dhat_star: max abs diff 0.00465
  gbsg/KM_diff_band_modtau/Zdhat_star: max abs diff 0.0778
  gbsg/KM_diff_band_modtau/surv0_star: max abs diff 0.00148
  gbsg/KM_diff_band_modtau/surv1_star: max abs diff 0.00438
  gbsg/KM_diff_band_modtau/c_alpha_band: max abs diff 0.00285
  gbsg/KM_diff_band_modtau/sb_lower: max abs diff 0.00331
  gbsg/KM_diff_band_modtau/sb_upper: max abs diff 0.00244
  gbsg_er_le0/df_counting_q25/quantile_check/time: NA pattern differs
  gbsg_er_le0/df_counting_q25/quantile_check/lower: NA pattern differs
  gbsg_er_le0/df_counting_q25/quantile_check/upper: NA pattern differs
  gbsg_er_le0/KM_diff_band_modtau/surv0: max abs diff 0.0132
  gbsg_er_le0/KM_diff_band_modtau/sig2_surv0: max abs diff 0.000167
  gbsg_er_le0/KM_diff_band_modtau/surv1: max abs diff 0.0494
  gbsg_er_le0/KM_diff_band_modtau/sig2_surv1: max abs diff 0.00103
  gbsg_er_le0/KM_diff_band_modtau/dhat: max abs diff 0.0388
  gbsg_er_le0/KM_diff_band_modtau/sig2_dhat: max abs diff 0.00114
  gbsg_er_le0/KM_diff_band_modtau/lower: max abs diff 0.0413
  gbsg_er_le0/KM_diff_band_modtau/upper: max abs diff 0.0437
  gbsg_er_le0/KM_diff_band_modtau/dhat_star: max abs diff 0.0225
  gbsg_er_le0/KM_diff_band_modtau/Zdhat_star: max abs diff 0.18
  gbsg_er_le0/KM_diff_band_modtau/surv0_star: max abs diff 0.00602
  gbsg_er_le0/KM_diff_band_modtau/surv1_star: max abs diff 0.0203
  gbsg_er_le0/KM_diff_band_modtau/c_alpha_band: max abs diff 0.00479
  gbsg_er_le0/KM_diff_band_modtau/sb_lower: max abs diff 0.0421
  gbsg_er_le0/KM_diff_band_modtau/sb_upper: max abs diff 0.0457
  gbsg_er_le0/wlr_dhat/cov_wlr_dhat: max abs diff 0.00748
  gbsg_er_le0/wlr_dhat/sig2_dhat: max abs diff 0.000836
  gbsg_er_gt0/df_counting_q25/quantile_check/time: NA pattern differs
  gbsg_er_gt0/df_counting_q25/quantile_check/lower: NA pattern differs
  gbsg_er_gt0/df_counting_q25/quantile_check/upper: NA pattern differs
  gbsg_er_gt0/KM_diff_band_modtau/surv0: max abs diff 0.00549
  gbsg_er_gt0/KM_diff_band_modtau/sig2_surv0: max abs diff 1.62e-05
  gbsg_er_gt0/KM_diff_band_modtau/surv1: max abs diff 0.00703
  gbsg_er_gt0/KM_diff_band_modtau/sig2_surv1: max abs diff 9.1e-05
  gbsg_er_gt0/KM_diff_band_modtau/dhat: max abs diff 0.00211
  gbsg_er_gt0/KM_diff_band_modtau/sig2_dhat: max abs diff 0.000107
  gbsg_er_gt0/KM_diff_band_modtau/lower: max abs diff 0.00329
  gbsg_er_gt0/KM_diff_band_modtau/upper: max abs diff 0.00278
  gbsg_er_gt0/KM_diff_band_modtau/dhat_star: max abs diff 0.00497
  gbsg_er_gt0/KM_diff_band_modtau/Zdhat_star: max abs diff 0.0833
  gbsg_er_gt0/KM_diff_band_modtau/surv0_star: max abs diff 0.00132
  gbsg_er_gt0/KM_diff_band_modtau/surv1_star: max abs diff 0.00466
  gbsg_er_gt0/KM_diff_band_modtau/c_alpha_band: max abs diff 0.00106
  gbsg_er_gt0/KM_diff_band_modtau/sb_lower: max abs diff 0.00418
  gbsg_er_gt0/KM_diff_band_modtau/sb_upper: max abs diff 0.00309
  sim_cont/df_counting_q25/quantile_check/time: max abs diff 17.5
  sim_cont/df_counting_q25/quantile_check/lower: max abs diff 13.3
  sim_cont/df_counting_q25/quantile_check/upper: NA pattern differs
  sim_cont/KM_diff_band_modtau/surv0: max abs diff 0.0284
  sim_cont/KM_diff_band_modtau/sig2_surv0: max abs diff 0.000349
  sim_cont/KM_diff_band_modtau/surv1: max abs diff 0.0214
  sim_cont/KM_diff_band_modtau/sig2_surv1: max abs diff 0.000143
  sim_cont/KM_diff_band_modtau/dhat: max abs diff 0.0164
  sim_cont/KM_diff_band_modtau/sig2_dhat: max abs diff 0.000435
  sim_cont/KM_diff_band_modtau/lower: max abs diff 0.0137
  sim_cont/KM_diff_band_modtau/upper: max abs diff 0.0209
  sim_cont/KM_diff_band_modtau/dhat_star: max abs diff 0.0181
  sim_cont/KM_diff_band_modtau/Zdhat_star: max abs diff 0.207
  sim_cont/KM_diff_band_modtau/surv0_star: max abs diff 0.0179
  sim_cont/KM_diff_band_modtau/surv1_star: max abs diff 0.00829
  sim_cont/KM_diff_band_modtau/c_alpha_band: max abs diff 0.0333
  sim_cont/KM_diff_band_modtau/sb_lower: max abs diff 0.0186
  sim_cont/KM_diff_band_modtau/sb_upper: max abs diff 0.0206
  sim_cont/wlr_dhat/cov_wlr_dhat: max abs diff 0.00275
  sim_cont/wlr_dhat/sig2_dhat: max abs diff 0.000102
  sim_ties/df_counting_q25/quantile_check/time: max abs diff 21.5
  sim_ties/df_counting_q25/quantile_check/lower: max abs diff 12.5
  sim_ties/df_counting_q25/quantile_check/upper: NA pattern differs
  sim_ties/wlr_dhat/cov_wlr_dhat: max abs diff 0.0108
  sim_ties/wlr_dhat/sig2_dhat: max abs diff 0.000401
  sim_weighted/df_counting_q25/quantile_check/time: max abs diff 15
  sim_weighted/df_counting_q25/quantile_check/lower: max abs diff 13
  sim_weighted/df_counting_q25/quantile_check/upper: NA pattern differs
  sim_weighted/KM_diff_band_modtau/surv0: max abs diff 0.0255
  sim_weighted/KM_diff_band_modtau/sig2_surv0: max abs diff 9.38e-05
  sim_weighted/KM_diff_band_modtau/surv1: max abs diff 0.0211
  sim_weighted/KM_diff_band_modtau/sig2_surv1: max abs diff 9.82e-05
  sim_weighted/KM_diff_band_modtau/dhat: max abs diff 0.00676
  sim_weighted/KM_diff_band_modtau/sig2_dhat: max abs diff 0.000155
  sim_weighted/KM_diff_band_modtau/lower: max abs diff 0.00572
  sim_weighted/KM_diff_band_modtau/upper: max abs diff 0.00899
  sim_weighted/KM_diff_band_modtau/dhat_star: max abs diff 0.0177
  sim_weighted/KM_diff_band_modtau/Zdhat_star: max abs diff 0.269
  sim_weighted/KM_diff_band_modtau/surv0_star: max abs diff 0.0168
  sim_weighted/KM_diff_band_modtau/surv1_star: max abs diff 0.00812
  sim_weighted/KM_diff_band_modtau/c_alpha_band: max abs diff 0.0425
  sim_weighted/KM_diff_band_modtau/sb_lower: max abs diff 0.00914
  sim_weighted/KM_diff_band_modtau/sb_upper: max abs diff 0.00787
  sim_weighted/wlr_dhat/cov_wlr_dhat: max abs diff 0.00199
  sim_weighted/wlr_dhat/sig2_dhat: max abs diff 4.75e-05
  sim_strata/df_counting_q25/quantile_check/time: max abs diff 20.1
  sim_strata/df_counting_q25/quantile_check/lower: max abs diff 15.4
  sim_strata/df_counting_q25/quantile_check/upper: NA pattern differs
  sim_strata/KM_diff_band_modtau/surv0: max abs diff 0.0278
  sim_strata/KM_diff_band_modtau/sig2_surv0: max abs diff 8.54e-05
  sim_strata/KM_diff_band_modtau/surv1: max abs diff 0.0151
  sim_strata/KM_diff_band_modtau/sig2_surv1: max abs diff 6.54e-05
  sim_strata/KM_diff_band_modtau/dhat: max abs diff 0.0127
  sim_strata/KM_diff_band_modtau/sig2_dhat: max abs diff 0.000144
  sim_strata/KM_diff_band_modtau/lower: max abs diff 0.0136
  sim_strata/KM_diff_band_modtau/upper: max abs diff 0.0118
  sim_strata/KM_diff_band_modtau/dhat_star: max abs diff 0.0133
  sim_strata/KM_diff_band_modtau/Zdhat_star: max abs diff 0.191
  sim_strata/KM_diff_band_modtau/surv0_star: max abs diff 0.0144
  sim_strata/KM_diff_band_modtau/surv1_star: max abs diff 0.00373
  sim_strata/KM_diff_band_modtau/c_alpha_band: max abs diff 0.0407
  sim_strata/KM_diff_band_modtau/sb_lower: max abs diff 0.0167
  sim_strata/KM_diff_band_modtau/sb_upper: max abs diff 0.00869
  sim_strata/wlr_dhat/cov_wlr_dhat: max abs diff 0.00193
  sim_strata/wlr_dhat/sig2_dhat: max abs diff 4.63e-05
B5 old/new ratios:
  veteran       sig2_dhat     old 0.00696039 new 0.00717291 old/new 0.970372
  veteran       cov_wlr_dhat  old 0.336652 new 0.341753 old/new 0.985075
  veteran       cor_wlr_dhat  old 0.731734 new 0.731734 old/new 1.000000
  gbsg          sig2_dhat     old 0.00156313 new 0.00156313 old/new 1.000000
  gbsg          cov_wlr_dhat  old 0.28045 new 0.28045 old/new 1.000000
  gbsg          cor_wlr_dhat  old 0.841932 new 0.841932 old/new 1.000000
  gbsg_er_le0   sig2_dhat     old 0.0165187 new 0.017355 old/new 0.951814
  gbsg_er_le0   cov_wlr_dhat  old 0.299372 new 0.306856 old/new 0.975610
  gbsg_er_le0   cor_wlr_dhat  old 0.840606 new 0.840606 old/new 1.000000
  gbsg_er_gt0   sig2_dhat     old 0.0017199 new 0.0017199 old/new 1.000000
  gbsg_er_gt0   cov_wlr_dhat  old 0.273051 new 0.273051 old/new 1.000000
  gbsg_er_gt0   cor_wlr_dhat  old 0.842115 new 0.842115 old/new 1.000000
  sim_cont      sig2_dhat     old 0.00504629 new 0.00514875 old/new 0.980100
  sim_cont      cov_wlr_dhat  old 0.272497 new 0.275249 old/new 0.990000
  sim_cont      cor_wlr_dhat  old 0.711725 new 0.711725 old/new 1.000000
  sim_ties      sig2_dhat     old 0.00456888 new 0.00497035 old/new 0.919226
  sim_ties      cov_wlr_dhat  old 0.250373 new 0.261142 old/new 0.958763
  sim_ties      cor_wlr_dhat  old 0.708103 new 0.708103 old/new 1.000000
  sim_weighted  sig2_dhat     old 0.00352589 new 0.00357337 old/new 0.986711
  sim_weighted  cov_wlr_dhat  old 0.29581 new 0.297796 old/new 0.993333
  sim_weighted  cor_wlr_dhat  old 0.733596 new 0.733596 old/new 1.000000
  sim_strata    sig2_dhat     old 0.00343809 new 0.00348439 old/new 0.986711
  sim_strata    cov_wlr_dhat  old 0.287845 new 0.289777 old/new 0.993333
  sim_strata    cor_wlr_dhat  old 0.745232 new 0.745232 old/new 1.000000
GOLDEN: identical outside the whitelist
OK 4b
OK 5a
OK 5b
OK 6a
OK 6b
OK 7 (if no FAIL above)
c5fddce B2: remove tests/testthat/Rplots.pdf committed by mistake with the B2 tests
7c3a87a stage2: verification artifacts (pressure test, plot comparison after, R CMD check, forestsearch, N/S after-checks)
8696d62 B4: NA stratum curve for an arm absent from a stratum
d165614 B4: call graphics::plot.new() in df_counting()'s SE-comparison plot
b68df5f stage2: package infrastructure (tests, DESCRIPTION 0.2.0, NEWS.md)
060da51 B10: KM curves drawn as step functions from time zero, with a step band
06ec1dc B5: inference details
f1b7816 B3: quantiles agree with survfit
f382b61 B1: Kaplan-Meier at arbitrary time points in KM_diff()
3b85230 B2: exact curve values and CI limits
32af3f2 B4: input validation and degenerate data
cb4c794 V: housekeeping (moves existing code)
135b2ff B6: restore the caller's RNG state; no cat() in get_validated_weights()
4cf5640 B7: import survival::strata
9f2a76b stage2: baseline capture (golden, API, plot comparison before) -- G1 passed
921bb9d M1b: weighted Greenwood added to the weighted-SE Monte Carlo
c0b5d35 stage2: add task doc
```
