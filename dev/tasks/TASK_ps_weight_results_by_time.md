# TASK: propensity-weight simulation results by time — gap run, article section, final check

Written 2026-11-20. Executor: Claude Code (CC). Approved by Larry ("comprehensive option"): the gap run, the results section in the estimated-weights article, the pointer fix, and the final rebuild and check.

## 0. Rules

- Branch `km-fixes-0.2.0`; HEAD must contain `37a0ac0`. If not, stop.
- First action: copy this file to `dev/tasks/TASK_ps_weight_results_by_time.md` and commit it.
- **No changes under `R/`, `man/`, `tests/`, `NAMESPACE`, `DESCRIPTION`.**
- Allowed paths: `vignettes/articles/estimated_weights_methodology.Rmd`; new folder `vignettes/articles/sim_results/` (one summary CSV added; the existing `wald_vs_logrank_full.rds` unchanged); `vignettes/weightedsurv_examples.Rmd` and the other `.Rmd` files under `vignettes/` **for `dev/` pointer sentences only** (section 4); `cran-comments.md` (only if the check result changes); `dev/`.
- Never push. Never merge. No tag. No submission. No site deployment.
- Do not edit the committed Stage 2e or OC scripts or results under `dev/audit/stage2e/` and `dev/audit/oc/`. New work goes in `dev/audit/ps_by_time/`.
- Do not rerun the forestsearch check (no `R/` change). Do not rerun any stored simulation.
- Gates (section 6) are stop-on-failure: if one fails, stop, write the report, do not ask.
- Report findings only: no recommendation, no follow-up task proposals.

**Approved compute (once):** 2 cells × 2,000 replicates, N = 400, both `weight_variance` options in each replicate. Expected 2–4 minutes at 12 workers (Stage 2e timings: 0.42–0.68 core-seconds per replicate at N = 400). Run a 50-replicate pilot first; if it projects more than 15 minutes for the full run, stop (G3). No bootstrap.

## 1. Sources

- `~/Downloads/REPORT_ps_weight_simulation_inventory.md` and `~/Downloads/ps_weight_simulation_table.csv` (read-only inventory of stored results, 2026-11-20). Copy both into `dev/audit/ps_by_time/` and commit them as the inventory record.
- `dev/audit/stage2e/oc_common.R`, `oc_sim.R`, `oc_summarise.R`, `m7_summary.csv`; `dev/audit/oc/oc_summary.csv`; `dev/audit/REPORT_km_fixes_stage2e.md`.

## 2. Gap run: treated ESS ratio near 0.9

The regime of the restricted Rotterdam example (treated ESS ratio 0.93, control 0.98) has no stored full two-arm cell.

**Design.** Identical to M7 cell 31 (`Es-400`) and cell 32 (`Es-400-null`) in every respect except the propensity model: control hazard 0.04 per month, censoring uniform(0, 72), X ~ N(0, 1) in both assignment and outcome as in `oc_common.R`, N = 400, stabilized weights from `ps_weights()`, both options on the same data and weights, the same per-replicate content as M7 (S0, S1, difference and RMST difference at months 1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48; difference and RMST bands on [3, 42]; log HR FH(0,0) and FH(0,1); log-rank rejection; ESS ratio by arm).

| New cell | Scenario |
|---|---|
| `Es9-400` | PH, HR = exp(−0.4) = 0.670 |
| `Es9-400-null` | null, HR = 1 |

**Propensity model.** `logit e = a + b X`. Choose `b` so that the treated ESS ratio is about 0.91 (target range 0.88–0.94) and the control ESS ratio is at least 0.95, and `a` so that P(treated) matches cells 31/32. Calibrate from the weights alone on one large draw of X (for example 10^6), using the true propensity; no survival fitting. Record `a`, `b`, P(treated), and the population ESS ratios. Then fix them for the run.

**Code.** A new script in `dev/audit/ps_by_time/` that sources `dev/audit/stage2e/oc_common.R` and reuses its functions with the new cell parameters; do not edit the stage 2e scripts. Package: the build of HEAD (`R/` unchanged since `ac76f82`; record the commit). Seeds recorded per cell. Summarise with the same rules as `oc_summarise.R` (truth, evaluable replicates, coverage, MCSE, length).

**Storage.** Commit the summary. Commit per-replicate files too if their total is under 5 MB; otherwise keep them in an ignored `raw/` folder and record their size and location.

## 3. Results section in the estimated-weights article

Replace the month-18 evaluation table of section 7 of `vignettes/articles/estimated_weights_methodology.Rmd` with a results section built from one stored summary file.

**Stored file.** `vignettes/articles/sim_results/ps_weight_results_by_time.csv`: one row per cell × estimand × time point, columns as the inventory CSV (true value, mean estimate, bias, bias MCSE, empirical SD, mean SE known / estimated, coverage known / estimated with MCSE, mean length known / estimated, evaluable replicates, ESS ratios). Rows for these cells only:

| Regime (treated ESS ratio) | Cells | Weights |
|---|---|---|
| Near 0.9 | new `Es9-400`, `Es9-400-null` | estimated, stabilized |
| About 0.69 | M7 31 `Es-400`, M7 32 `Es-400-null` | estimated, stabilized |
| Strong confounding, 0.33–0.39 | M7 35 `Fe-400`, M7 37 `Fe-400-null` | estimated, ipw |
| 14% treated, 0.35 | M7 39 `V-1450-ph`, M7 38 `V-1450-null` | estimated, stabilized |
| Reference: true weights | OC 5 `D-400`, OC 10 `F-400` | true, ipw (`"known"` only) |

Existing rows are copied from `m7_summary.csv` / `oc_summary.csv` unchanged. Add the bootstrap benchmark rows (cells 8, 38, 39; first 500 replicates; difference at month 18, RMST at month 36, log HR) as stored.

The article reads only this CSV (no read from `dev/`). It is not in the package tarball (`vignettes/articles/` is excluded).

**Content of the section.**

1. **Design paragraph.** The data-generating model, the regimes, N, replicates, the propensity models, how the truth was obtained, Monte Carlo standard error of a coverage (0.005 at 2,000 replicates). For each reported month, its percentile in the distribution of observed event times under the design (from the formula; the inventory gives the 10th and 90th percentiles). Provenance in plain text: run for weightedsurv 0.2.0, package commit, scripts and records kept in the package repository under `dev/audit/` (no URL).
2. **Table 1: survival difference.** Months 1, 6, 12, 24, 36, 48. Rows regime × scenario × month; columns: true, bias, empirical SD, SE known, SE estimated, coverage known (MCSE), coverage estimated (MCSE), length known, length estimated. Reference rows show `"known"` only.
3. **Figure 1.** Coverage of the survival difference at all 11 stored months, one panel per regime and scenario, two lines (`"known"`, `"estimated"`), horizontal line at 0.95 with ±2 MCSE. Same data as the CSV.
4. **Table 2: RMST difference.** Same months and columns.
5. **Table 3: log hazard ratio,** FH(0,0) and FH(0,1): one row per cell, same columns without a time index.
6. **Table 4: simultaneous bands on [3, 42].** Coverage and mean width under both options.
7. **Table 5: refitting bootstrap** on the same first 500 replicates, beside `"known"` and `"estimated"` on those replicates.
8. **Reading** (findings only, short, after the tables): by regime and by time; where `"estimated"` is near nominal and where both options fall short; the bias of the estimate and how it grows with time under concentrated weights; that `"estimated"` on fitted weights tracks `"known"` on true weights in the same design; agreement with the bootstrap. Do not present either option as generally better on the strength of one regime. Every number in the reading is inline R from the CSV.

Section 7's existing text on conditions (correct model, overlap, effective sample size) stays. Remove "late in follow-up" claims that the tables do not show, or point them to the table that shows them.

## 4. Pointer fix

`vignettes/weightedsurv_examples.Rmd` (in the package tarball) points to `dev/audit/REPORT_km_fixes_stage2e.md` for the evaluation. Replace with a pointer to the results section of the estimated-weights article. Then `grep -rn "dev/" vignettes/ --include="*.Rmd"` and do the same for every other prose pointer to a `dev/` path; the only `dev/` mention allowed afterwards is the provenance sentence of section 3.1. List each change.

## 5. Final rebuild and check

1. Render the estimated-weights article and the main vignette with `options(warn = 2)`.
2. Build the tarball from a clean `git clone` at the final commit; `R CMD check --as-cran` (`DISPLAY` unset, pandoc on the path as in Stage 3).
3. `pkgdown::check_pkgdown()`; `pkgdown::build_site()` locally (the simulations article reads its stored file; it must not run the full simulation). Do not commit `docs/`. Do not deploy.
4. Update `cran-comments.md` only if the check result differs from Stage 3.

## 6. Gates (stop-on-failure; record PASS/FAIL with evidence)

| # | Check |
|---|---|
| G1 | `git diff --name-only 37a0ac0..HEAD` lists only allowed paths. `git diff 37a0ac0..HEAD -- R/ man/ tests/ NAMESPACE DESCRIPTION dev/audit/stage2e/ dev/audit/oc/` is empty. `vignettes/articles/sim_results/wald_vs_logrank_full.rds` unchanged. |
| G2 | Calibration: in both new cells the mean realised treated ESS ratio is in [0.88, 0.94] and the control ratio ≥ 0.95. P(treated) within 0.02 of cells 31/32. |
| G3 | Compute: pilot projection ≤ 15 minutes; run once; 2 cells × 2,000 replicates; failed replicates ≤ 1% per cell (count recorded); wall-clock recorded. |
| G4 | Truth: the true values of `Es9-400` equal those of cell 31 at every month and estimand (same outcome model; `identical()` or difference < 1e-12); `Es9-400-null` true difference, RMST difference and log HR are 0. |
| G5 | Copied rows: every row of the stored CSV taken from `m7_summary.csv` or `oc_summary.csv` equals its source exactly (script check, count of rows compared). |
| G6 | Article: renders with `options(warn = 2)` in under 60 seconds; no code chunk reads a path under `dev/`; Tables 1–5 and Figure 1 present; the earlier month-18-only table is gone. |
| G7 | Numerals in the new prose: each is inline R from the CSV, a design constant, or a formula constant (list, as in earlier tasks). |
| G8 | `grep -rn "dev/" vignettes/ --include="*.Rmd"` returns only the provenance sentence(s) of section 3.1, listed. `vignettes/weightedsurv_examples.Rmd` contains no `dev/` path. |
| G9 | `R CMD check --as-cran` on the clean-clone tarball: 0 errors, 0 warnings, NOTEs no more than the two of Stage 3; tests `FAIL 0`, `PASS 1691`. Neither file in `vignettes/articles/sim_results/` is in the tarball. |
| G10 | `check_pkgdown()` no problems; `build_site()` completes in under 5 minutes; the estimated-weights article page has the new section; the simulations article did not run the full simulation. `git status` shows nothing from `docs/`. |

## 7. Report

`dev/audit/REPORT_ps_weight_results_by_time.md`:

1. Gate table.
2. Calibration: `a`, `b`, P(treated), population and realised ESS ratios.
3. **One comparative table for the new cells** (survival difference and RMST difference at months 1, 6, 12, 24, 36, 48; log HR; bands): true, bias, empirical SD, SE and coverage and length under both options. Then a short plain-language reading, findings only.
4. Changes to the article (section by section) and the pointer fixes.
5. Check result, site build, tarball path, size and sha256.
6. State of the branch: HEAD, commits since `37a0ac0`, `git status`, commits ahead of the remote.

After the final commit, copy to `~/Downloads/` and print the paths: the report, the rendered HTML of the estimated-weights article, `ps_weight_results_by_time.csv`, and the new `weightedsurv_0.2.0.tar.gz` (it replaces the Stage 3 tarball in `~/Downloads/`). Copy `cran-comments.md` too if it changed. Do not push.
