# REPORT — propensity-weight simulation results by time

Task: `dev/tasks/TASK_ps_weight_results_by_time.md` (2026-11-20). Branch `km-fixes-0.2.0`, started at `37a0ac0`.
Work folder: `dev/audit/ps_by_time/`. Findings only.

## 1. Gates

| # | Result | Evidence |
|---|---|---|
| G1 | PASS | `git diff --name-only 37a0ac0..HEAD`: only `dev/tasks/TASK_ps_weight_results_by_time.md`, `dev/audit/ps_by_time/*`, `dev/audit/REPORT_ps_weight_results_by_time.md`, `vignettes/articles/estimated_weights_methodology.Rmd`, `vignettes/articles/sim_results/ps_weight_results_by_time.csv` (new), and the pointer sentences of `vignettes/weightedsurv_examples.Rmd`, `vignettes/articles/weightedkm_methodology.Rmd`, `vignettes/articles/weightedcox_methodology.Rmd`. `git diff 37a0ac0..HEAD -- R/ man/ tests/ NAMESPACE DESCRIPTION dev/audit/stage2e/ dev/audit/oc/`: empty (0 lines). `vignettes/articles/sim_results/wald_vs_logrank_full.rds`: no diff. `cran-comments.md` unchanged. `git ls-files docs`: 0. Rechecked after the final commit (`g1_paths.txt`). |
| G2 | PASS | `ps9_gates.txt`. Mean realised ESS ratio over 2,000 replicates: Es9-400 treated 0.9034, control 0.9550; Es9-400-null treated 0.9032, control 0.9547 (limits [0.88, 0.94] and ≥ 0.95; the control ratio is 0.0050 and 0.0047 above its limit). Realised P(treated), data regenerated without package calls: cell 31 0.3985, cell 41 0.3973, cell 42 0.3969 (largest difference 0.0016; limit 0.02). Population: 0.3977 against 0.3978 (`calibration.txt`). |
| G3 | PASS | Pilot (`pilot.txt`): 25 replicates per cell (50 in all) on 12 workers, 1.7 s and 1.5 s; projected full run 4.3 min (limit 15). Full run once (`ps9_stdout.txt`): 2 cells × 2,000 replicates, cell 41 84.0 s, cell 42 78.5 s, wall-clock 162.8 s (2.7 min), 12 workers. Failed replicates (any failed call, or `ps_weights()` error): 0 and 0. |
| G4 | PASS | `ps9_gates.txt`: 92 true values of Es9-400 (S0, S1, difference, RMST at 11 months, log HR; both options) compared with cell 31: 92 identical, largest difference 0. Es9-400-null: 48 true values of the difference, RMST difference and log HR, all exactly 0. |
| G5 | PASS | `ps9_gates.txt`: the stored CSV read back; 385 copied rows (M7 286, OC 99; 9 with bootstrap columns), 5,576 values compared with `identical()` against `m7_summary.csv` / `oc_summary.csv`: 0 differ. Cross-check against the inventory CSV: 385 rows matched, 6,826 values, 0 differ beyond 1e-12 relative; 54 inventory values not carried (true-weight bootstrap columns of OC cells 5 and 10, outside the benchmark cells 8, 38, 39 named in the task). |
| G6 | PASS | `g6_render.txt`: article rendered with `options(warn = 2)` in 1.6 s (main vignette 7.8 s), weightedsurv 0.2.0 built from HEAD in a temporary library. `g6_g7_checks.txt`: 17 code chunks, none contains `dev/`; the only read call is `read.csv("sim_results/ps_weight_results_by_time.csv")`; Tables 1–5 and Figure 1 present in the HTML; the month-18 table and its "Source" line are gone. |
| G7 | PASS | `g6_g7_checks.txt` lists the 99 numerals of the new prose outside inline R. Every number taken from results is inline R from the CSV. The listed numerals are: design constants (N(0, 1); 0.6, 0.8; hazard 0.04; −0.4 and HR 0.670; uniform(0, 72); propensity coefficients −0.44, 0.51, −0.5, 1, 2, −2.3, 1.3; N = 400, 1,450; calibration target 0.91; 500 bootstrap replicates, 200 resamples; nominal 0.95; 2.5% / 97.5% quantiles; band range [3, 42]); stored grid months (1, 3, 6, 24, 30, 36, 48); percentile labels 10th / 90th; regime names as stored in the CSV column `regime` (Near 0.9, About 0.69, 14% treated); cell names (Es-400, Fe-400, …); table, figure and section numbers; FH(0,0); package version 0.2.0 and commit ids (f55001f, c195953, d9b7558). Formula values (MCSE, event-time percentiles, % treated of the 14%-treated design) are computed inline. |
| G8 | PASS | `grep -rn "dev/" vignettes/ --include="*.Rmd"` returns one line: `vignettes/articles/estimated_weights_methodology.Rmd:346`, the provenance sentence of section 7.1. `vignettes/weightedsurv_examples.Rmd`: no `dev/`. |
| G9 | PASS | `g9_00check.log`, `g9_testthat.Rout`, `g9_tarball.txt`. `R CMD check --as-cran` on the tarball built by `R CMD build` from a clean `git clone` at `a19e189` (`DISPLAY` unset; pandoc from `/Applications/quarto/bin/tools` on the `PATH`): `Status: 2 NOTEs`, 0 errors, 0 warnings; the two NOTEs are those of Stage 3 (README URL 404; HTML Tidy). Tests `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`. Tarball: 112 files, no `sim_results`, no `vignettes/articles`, no `dev/` (the only `.rds` is `build/vignette.rds`; no `.csv`). |
| G10 | PASS | `g10_check_pkgdown.txt`: "No problems found." `g10_build_site_log.txt`: `build_site()` completed in 1 min 32 s. `docs/articles/estimated_weights_methodology.html` has "Simulation results by time", Tables 1–5 and Figure 1. The simulations article read its stored file (page states "completed: 30,000; failed: 0"; the whole build took 1.5 min). `docs/` is ignored; `git status` shows nothing from it. Not deployed. |

## 2. Calibration

`calibrate.R` → `calibration.txt`. One draw of X ~ N(0, 1), n = 10^6, seed 20261120; true propensity; population ESS ratio of arm a computed as 1 / (P · E[1/e]) (treated) and 1 / ((1 − P) · E[1/(1 − e)]) (control); `a` solved so that P(treated) equals that of cells 31/32 (logit e = −0.5 + X), `b` so that the treated ratio is 0.91.

| | a | b | P(treated) | Population ESS ratio treated / control | Realised mean ESS ratio treated / control (2,000 replicates, fitted stabilized weights) |
|---|---|---|---|---|---|
| Solution | −0.438943 | 0.505434 | | | |
| Fixed for the run | −0.44 | 0.51 | 0.3977 (cells 31/32: 0.3978) | 0.9084 / 0.9579 | Es9-400 0.9034 / 0.9550; Es9-400-null 0.9032 / 0.9547 |
| Cells 31/32, for reference | −0.5 | 1 | 0.3978 | 0.6766 / 0.8307 | 0.69 / 0.83 |

Realised ESS ratio over replicates: treated 0.553–0.995, control 0.816–0.997.

Run: `ps9_sim.R` (sources `dev/audit/stage2e/oc_common.R`; evaluates the 22 function definitions of `dev/audit/stage2e/oc_sim.R` without its top-level code; neither file edited). Package: `pkgload::load_all()` of the working tree at `f55001f` (`R/` unchanged since `ac76f82`: `git diff ac76f82..HEAD -- R/` empty). Seeds: replicate r of cell k seeded `1e6 * k + r`, k = 41 (Es9-400), 42 (Es9-400-null), r = 1..2,000. Data generator: `gen()` of `oc_common.R` with e = plogis(−0.44 + 0.51 X), same draws in the same order. Per replicate: the M7 content of cells 31/32 (S0, S1, difference, RMST with pointwise intervals at the 11 months; difference and RMST bands on [3, 42]; log HR FH(0,0), FH(0,1); log-rank z; ESS ratio by arm) under both options, without the requirement-2 / `wlr_dhat_estimates()` extras of cells 31/32. Requirement 3 of M7 (estimates identical under the two options): 0 of 4,000 replicates failing. Summary: `ps9_summarise.R` evaluates the definitions and per-cell loop of `dev/audit/stage2e/oc_summarise.R` on `dev/audit/ps_by_time/raw/` → `ps9_summary.csv` (216 rows). Per-replicate files: 8.5 MB in all (over 5 MB), kept in `dev/audit/ps_by_time/raw/m7_cell41.rds`, `m7_cell42.rds` (ignored by `dev/audit/ps_by_time/.gitignore`, `raw*/`).

## 3. The new cells

Coverage of nominal 0.95 intervals with MCSE in brackets; 2,000 replicates (1,997 evaluable at month 48 in Es9-400-null, 2,000 elsewhere). Bias and empirical SD do not depend on the option. Bands: coverage and mean width. (`report_tables.R` → `report_tables.md`.)

| Cell | Estimand | True | Bias | Emp SD | SE known | SE est | Cov known (MCSE) | Cov est (MCSE) | Length known | Length est |
|---|---|---|---|---|---|---|---|---|---|---|
| Es9-400 | difference, month 1 | 0.013 | 0.0007 | 0.0182 | 0.0183 | 0.0176 | 0.945 (0.005) | 0.939 (0.005) | 0.072 | 0.069 |
| Es9-400 | difference, month 6 | 0.065 | 0.0001 | 0.0358 | 0.0399 | 0.0361 | 0.975 (0.004) | 0.953 (0.005) | 0.156 | 0.141 |
| Es9-400 | difference, month 12 | 0.106 | -0.0003 | 0.0452 | 0.0500 | 0.0443 | 0.962 (0.004) | 0.943 (0.005) | 0.196 | 0.174 |
| Es9-400 | difference, month 24 | 0.143 | 0.0002 | 0.0536 | 0.0578 | 0.0518 | 0.962 (0.004) | 0.934 (0.006) | 0.227 | 0.203 |
| Es9-400 | difference, month 36 | 0.144 | 0.0009 | 0.0563 | 0.0601 | 0.0554 | 0.965 (0.004) | 0.947 (0.005) | 0.236 | 0.217 |
| Es9-400 | difference, month 48 | 0.129 | 0.0008 | 0.0599 | 0.0617 | 0.0585 | 0.951 (0.005) | 0.935 (0.005) | 0.242 | 0.229 |
| Es9-400 | RMST, month 1 | 0.006 | -0.0002 | 0.0103 | 0.0101 | 0.0098 | 0.943 (0.005) | 0.936 (0.005) | 0.040 | 0.038 |
| Es9-400 | RMST, month 6 | 0.208 | -0.0008 | 0.1303 | 0.1434 | 0.1295 | 0.973 (0.004) | 0.954 (0.005) | 0.562 | 0.508 |
| Es9-400 | RMST, month 12 | 0.730 | -0.0072 | 0.3282 | 0.3718 | 0.3229 | 0.975 (0.003) | 0.952 (0.005) | 1.458 | 1.266 |
| Es9-400 | RMST, month 24 | 2.271 | -0.0243 | 0.7824 | 0.8989 | 0.7580 | 0.975 (0.003) | 0.941 (0.005) | 3.524 | 2.972 |
| Es9-400 | RMST, month 36 | 4.014 | -0.0351 | 1.2361 | 1.4343 | 1.2074 | 0.976 (0.003) | 0.944 (0.005) | 5.622 | 4.733 |
| Es9-400 | RMST, month 48 | 5.664 | -0.0456 | 1.6848 | 1.9480 | 1.6559 | 0.976 (0.003) | 0.945 (0.005) | 7.636 | 6.491 |
| Es9-400 | log HR FH(0,0) | -0.400 | -0.0030 | 0.1203 | 0.1387 | 0.1188 | 0.976 (0.003) | 0.953 (0.005) | 0.544 | 0.466 |
| Es9-400 | log HR FH(0,1) | -0.400 | -0.0056 | 0.1598 | 0.1659 | 0.1561 | 0.954 (0.005) | 0.945 (0.005) | 0.650 | 0.612 |
| Es9-400 | band [3, 42], difference (width) |  |  |  |  |  | 0.978 (0.003) | 0.955 (0.005) | 0.313 | 0.286 |
| Es9-400 | band [3, 42], RMST (width) |  |  |  |  |  | 0.987 (0.003) | 0.971 (0.004) | 3.234 | 2.787 |
| Es9-400-null | difference, month 1 | 0.000 | -0.0009 | 0.0195 | 0.0198 | 0.0190 | 0.952 (0.005) | 0.945 (0.005) | 0.078 | 0.074 |
| Es9-400-null | difference, month 6 | 0.000 | -0.0001 | 0.0386 | 0.0434 | 0.0387 | 0.977 (0.003) | 0.952 (0.005) | 0.170 | 0.152 |
| Es9-400-null | difference, month 12 | 0.000 | 0.0003 | 0.0477 | 0.0535 | 0.0469 | 0.975 (0.004) | 0.949 (0.005) | 0.210 | 0.184 |
| Es9-400-null | difference, month 24 | 0.000 | 0.0003 | 0.0543 | 0.0586 | 0.0526 | 0.968 (0.004) | 0.951 (0.005) | 0.230 | 0.206 |
| Es9-400-null | difference, month 36 | 0.000 | -0.0011 | 0.0542 | 0.0571 | 0.0532 | 0.959 (0.004) | 0.943 (0.005) | 0.224 | 0.208 |
| Es9-400-null | difference, month 48 | 0.000 | -0.0021 | 0.0539 | 0.0543 | 0.0520 | 0.944 (0.005) | 0.929 (0.006) | 0.213 | 0.204 |
| Es9-400-null | RMST, month 1 | 0.000 | -0.0005 | 0.0110 | 0.0110 | 0.0106 | 0.949 (0.005) | 0.943 (0.005) | 0.043 | 0.042 |
| Es9-400-null | RMST, month 6 | 0.000 | -0.0028 | 0.1398 | 0.1557 | 0.1387 | 0.970 (0.004) | 0.948 (0.005) | 0.611 | 0.544 |
| Es9-400-null | RMST, month 12 | 0.000 | -0.0027 | 0.3460 | 0.4008 | 0.3424 | 0.975 (0.004) | 0.950 (0.005) | 1.571 | 1.342 |
| Es9-400-null | RMST, month 24 | 0.000 | 0.0108 | 0.7984 | 0.9436 | 0.7860 | 0.984 (0.003) | 0.952 (0.005) | 3.699 | 3.081 |
| Es9-400-null | RMST, month 36 | 0.000 | 0.0025 | 1.2455 | 1.4562 | 1.2189 | 0.978 (0.003) | 0.945 (0.005) | 5.708 | 4.778 |
| Es9-400-null | RMST, month 48 | 0.000 | -0.0216 | 1.6598 | 1.9018 | 1.6144 | 0.976 (0.003) | 0.941 (0.005) | 7.455 | 6.329 |
| Es9-400-null | log HR FH(0,0) | 0.000 | 0.0037 | 0.1167 | 0.1332 | 0.1132 | 0.977 (0.003) | 0.943 (0.005) | 0.522 | 0.444 |
| Es9-400-null | log HR FH(0,1) | 0.000 | 0.0089 | 0.1569 | 0.1599 | 0.1505 | 0.951 (0.005) | 0.939 (0.005) | 0.627 | 0.590 |
| Es9-400-null | band [3, 42], difference (width) |  |  |  |  |  | 0.983 (0.003) | 0.956 (0.005) | 0.324 | 0.294 |
| Es9-400-null | band [3, 42], RMST (width) |  |  |  |  |  | 0.987 (0.003) | 0.967 (0.004) | 3.223 | 2.750 |

Reading:

- With treated ESS ratio about 0.90, `"estimated"` covers the survival difference at 0.929–0.956 across the 11 months of both cells; `"known"` at 0.944–0.977. At the six months above, `"estimated"` is within 2 MCSE of 0.95 (0.940–0.960) except at months 1 and 24 in Es9-400 (0.939, 0.934) and month 48 in Es9-400-null (0.929).
- The `"estimated"` standard error is 0.96–1.01 times the empirical SD of the difference; the `"known"` one 1.01–1.12 times. `"estimated"` intervals are 0.88–0.96 of the `"known"` length.
- The bias of the difference is within ±0.0023 at every month in both cells. The RMST bias in Es9-400 grows to −0.046 at month 48 (1.2 MCSE from zero; true value 5.66).
- RMST: `"estimated"` 0.936–0.954, `"known"` 0.943–0.984; log HR FH(0,0): 0.953 / 0.943 against 0.976 / 0.977; FH(0,1): 0.945 / 0.939 against 0.954 / 0.951. Bands: `"estimated"` 0.955–0.971, `"known"` 0.978–0.987.
- Against cells 31/32 (treated ESS ratio 0.69, same outcome model and treated proportion), the late decline is smaller: at month 48 `"estimated"` covers the difference at 0.935 / 0.929 against 0.928 / 0.900, and the bias at month 48 is 0.0008 / −0.0021 against −0.0058 / −0.0042.

## 4. Changes to the article and pointer fixes

`vignettes/articles/estimated_weights_methodology.Rmd`:

- **Section 7 (Conditions and properties).** The paragraphs on the correct model, overlap and effective sample size are unchanged. Removed: the paragraph "Evaluation evidence" (copied from the Stage 2e report), the month-18 table, its "Source" line and its five reading bullets, including the "late in follow-up" bullet that pointed to `dev/audit/stage2e/m7_summary.csv`. Added: subsection **7.1 Simulation results by time**:
  1. Design paragraph (data-generating model, five regimes with their propensity models, N, weight type and ESS ratios inline; the truth in closed form; evaluable replicates; MCSE 0.005 at 2,000 and 0.010 at 500, computed; percentile of each reported month among observed event times from the formula, control and treated: month 1 at 6 / 5, 6 at 30 / 26, 12 at 52 / 45, 24 at 79 / 72, 36 at 92 / 88, 48 at 97 / 96; 10th and 90th percentiles 1.8 and 33.7 / 2.2 and 38.6, as in the inventory). Provenance sentence (version 0.2.0, commits, records under `dev/audit/`, no URL).
  2. Table 1 (survival difference, months 1, 6, 12, 24, 36, 48; 60 rows), with a note on evaluable replicates.
  3. Figure 1 (coverage of the difference at the 11 months, 10 panels, `"known"` and `"estimated"`, 0.95 ± 2 MCSE band).
  4. Table 2 (RMST difference, same months and columns), Table 3 (log HR FH(0,0) and FH(0,1), one row per cell and statistic), Table 4 (bands on [3, 42], both options), Table 5 (refitting bootstrap on the first 500 replicates of cells 8, 38, 39 beside `"known"` and `"estimated"`), with a paragraph on how the bootstrap and cell 8 were computed.
  5. Reading by regime, true against estimated weights, log HR and bands, bootstrap; every number inline R from the CSV; qualitative statements guarded by `stopifnot()` in the chunk `simres-read`.
- **Section 6.** "The Stage 2e record (`dev/tasks/TASK_km_fixes_stage2e.md`) lists `draws = 0` …" → "The development record of this version lists `draws = 0` …".
- Sections 8 and 9 unchanged.

Pointer fixes (`grep -rn "dev/" vignettes/ --include="*.Rmd"` before: 7 lines; after: 1):

| File | Before | After |
|---|---|---|
| `vignettes/weightedsurv_examples.Rmd:527` | "(`dev/audit/REPORT_km_fixes_stage2e.md` in the package repository) … both fell short of nominal, mostly late in follow-up, …" | "(section "Simulation results by time" of the article `estimated_weights_methodology`) … both fell short of nominal, at the earliest or the latest months depending on the design (Table 1 and Figure 1 of that section), …" |
| `vignettes/articles/weightedkm_methodology.Rmd:201` | "the results are in `dev/audit/REPORT_km_fixes_stage2e.md` in the package repository." | "the results are in the section "Simulation results by time" of the article `estimated_weights_methodology`." |
| `vignettes/articles/weightedcox_methodology.Rmd:153` | "The evaluation is in `dev/audit/REPORT_km_fixes_stage2e.md` in the package repository." | "The evaluation is in the section "Simulation results by time" of the article `estimated_weights_methodology`." |
| `vignettes/articles/estimated_weights_methodology.Rmd` | section 6 `dev/tasks/` sentence; section 7 lines 291, 314, 320 | section 6 as above; section 7 lines removed with the table; one new provenance sentence (line 346) |

Stored file: `vignettes/articles/sim_results/ps_weight_results_by_time.csv`, 481 rows, the 43 columns of the inventory CSV (`regime` holds the task's regime names): 96 new rows (cells 41, 42), 286 M7 rows (cells 31, 32, 35, 37, 38, 39), 96 OC rows (cells 5, 10; `"known"` = `proposed`), 3 OC rows for cell 8 (bootstrap benchmark only). Estimands: S0, S1, difference, RMST at 11 months; log HR FH(0,0), FH(0,1) (not in 38, 39); difference and RMST bands. ESS ratios of OC cells from `dev/audit/oc/raw/` (not in `oc_summary.csv`). `wald_vs_logrank_full.rds` unchanged.

## 5. Check result, site build, tarball

- `R CMD check --as-cran`, tarball from a clean `git clone` at `a19e189` (the last commit touching package files; later commits touch `dev/` only, which `.Rbuildignore` excludes): `Status: 2 NOTEs`, 0 errors, 0 warnings; tests `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`; examples, `--run-donttest` examples and vignettes OK. Same result as Stage 3; `cran-comments.md` not changed. R 4.5.2, macOS, `DISPLAY` unset, pandoc from quarto on the `PATH`.
- `pkgdown::check_pkgdown()`: no problems. `pkgdown::build_site()`: completed locally in 1 min 32 s; not deployed; `docs/` not committed.
- Tarball: `weightedsurv_0.2.0.tar.gz`, 4,329,283 bytes, 112 files, sha256 `dc25652b58aade3ec9cf6648bb2c96bbc6fc2de717501d42b6feeb49553261fd`. Copied to `~/Downloads/weightedsurv_0.2.0.tar.gz` (replaces the Stage 3 tarball).

## 6. State of the branch

At the time of writing (before this report's commit): HEAD `a19e189`; commits since `37a0ac0`:

```
a19e189 vignettes: estimated-weights article section 7.1 'Simulation results by time' ...; dev/ pointers ...; G6, G7 outputs
a13a2ea ps_by_time: gap run Es9-400 / Es9-400-null (...), summary, gates G2-G5; .../ps_weight_results_by_time.csv
f55001f ps_by_time: inventory record (...) copied as given
056f12b task: TASK_ps_weight_results_by_time.md (...)
```

This report, the G9/G10 outputs, `report_tables.*` and `g1_paths.txt` are committed in one further commit (`release`-style message beginning "ps_by_time: REPORT"). Upstream `origin/km-fixes-0.2.0`; nothing pushed; after that commit the branch is 5 commits ahead of it. `git status` after the commit: clean (the per-replicate files in `dev/audit/ps_by_time/raw/` and `docs/` are ignored). No merge, tag, submission or deployment.
