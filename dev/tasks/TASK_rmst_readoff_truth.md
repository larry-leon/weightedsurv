# TASK: one truth convention for the RMST rows of the simulation results — weightedsurv 0.2.0

Written 2026-11-20. Executor: Claude Code (CC). Approved by Larry: use the read-off truth R(t_j) for every RMST row (recommendation of chat, 2026-11-20).

## 0. Rules

- Branch `release-0.2.0-robust`; HEAD must contain `8107811` and the report commit after it. If not, stop.
- First action: copy this file to `dev/tasks/TASK_rmst_readoff_truth.md` and commit it.
- **Allowed paths:** `vignettes/articles/estimated_weights_methodology.Rmd`, `vignettes/articles/sim_results/ps_weight_results_by_time.csv`, `dev/`. Nothing else. Both vignette paths are excluded from the package tarball, so the tarball in `~/Downloads` (sha256 `a7f38c2e…`) and its checks stand; do not rebuild it.
- Do not edit anything under `dev/audit/stage2e/`, `dev/audit/oc/`, or the committed scripts of `dev/audit/ps_by_time/`. New scripts go in `dev/audit/rmst_readoff/`.
- **No simulation, no bootstrap, no package analysis of any replicate.** Allowed: reading stored per-replicate files; evaluating the closed-form truth function of `dev/audit/stage2e/oc_common.R`; regenerating replicate *data* from recorded seeds with the data generator only (no package call), checked against the stored data checksums.
- Never push. Never merge. No tag. No submission. No site deployment.
- Gates (section 4) are stop-on-failure: if one fails, stop, write the report, do not ask.
- Report findings only: no recommendation, no follow-up task proposals.

## 1. The issue (from `dev/audit/REPORT_release_robustness_0.2.0.md` section 5)

Each replicate's RMST "at month t" is the estimate at the last grid time t_j ≤ t (the default `KM_diff()` grid is the event times). The OC summaries score it against R(t_j) (`dev/audit/oc/oc_sim.R:210-211`, `oc_summarise.R:107`); the M7 summaries, and the new cells 41/42 through `ps9_summarise.R`, score it against the nominal R(t). The fixed-truth bias therefore includes the read-off gap R(t_j) − R(t) < 0. The stored CSV mixes the two conventions in Table 2 of the article, and its `mean_estimate` column is `truth + bias`, which is not the mean of the estimates for per-replicate truth rows. The previous task moved the cell 8 month 36 row to the fixed truth; that change is superseded here.

The survival-difference, S0 and S1 rows are not affected (a Kaplan–Meier curve is constant between event times, so its value at t equals its value at t_j).

## 2. Work

1. **Recover t_j per replicate and month**, for every cell with RMST rows in the CSV (M7 31, 32, 35, 37, 38, 39; new 41, 42; OC 5, 8, 10), in this order of preference:
   a. a stored grid time, grid index or R(t_j) in the per-replicate file;
   b. regenerate the replicate's data from its recorded seed with the generator of `oc_common.R` (no package call), confirm the stored data checksum (`$chk` or equivalent) matches, and apply the grid and read-off rule exactly as coded in the simulation script that produced the cell (read the script; record the lines).
   Record which route each cell used.
2. **Recompute every RMST row** (estimand `R`, all cells, all 11 months) from the per-replicate estimates and limits:
   - `true_value`: the nominal R(t), unchanged (label).
   - new column `mean_truth_readoff`: mean over evaluable replicates of R(t_j).
   - `bias` = mean(R̂(t_j) − R(t_j)); `bias_mcse` its Monte Carlo SE.
   - `mean_estimate` = mean of R̂(t_j) (the actual mean, for every RMST row).
   - `coverage_known`, `coverage_estimated` (and MCSEs) = coverage of R(t_j).
   - `empirical_sd`, mean SEs, mean lengths: recomputed; they should equal the stored values.
   - Bootstrap columns of the benchmark rows (cells 8, 38, 39, RMST month 36): recompute coverage of R(t_j) from `boot_cell*.rds` on the same replicates where stored; record if a column cannot be recomputed and leave it unchanged.
   - Add `mean_truth_readoff` as an empty column for non-RMST rows. Every non-RMST row is otherwise byte-identical.
3. **Band rows** (`R_band_3_42`): read the scripts and record against which truth the RMST band coverage was computed (grid-point truth or not). Change nothing for bands unless their truth is the nominal-month truth applied at a different time; then record it and leave unchanged.
4. **Fallback (pre-approved).** If t_j cannot be recovered for some cell by either route, use the fixed truth R(t) for **all** RMST rows instead: recompute OC rows 5, 8, 10 against R(t) from their per-replicate estimates, leave M7 and new rows as stored (except `mean_estimate`, set to the actual mean), and state the read-off gap in the article (section 3). Record which cell forced the fallback.
5. **Article** (`estimated_weights_methodology.Rmd`, section 7.1):
   - Design paragraph: the RMST at month t is read at the last grid time t_j ≤ t; bias and coverage are against R(t_j) (or, under the fallback, against R(t), with one sentence that the bias then includes the read-off gap R(t_j) − R(t), of order δ(t)·(t − t_j)).
   - Table 2: a note line stating the convention. Keep the displayed columns; the code must handle the new CSV column.
   - Reading: rewrite every RMST sentence from the recomputed values (inline R). Remove any statement that attributes the read-off gap to the estimator (for example "bias accumulates with t" if the recomputed bias does not show it). List each RMST sentence before and after.
   - Nothing else in the article changes.

## 3. Checks

- Render the article with `options(warn = 2)`; `pkgdown::check_pkgdown()`; `pkgdown::build_article("articles/estimated_weights_methodology")` into a scratch destination. Do not commit `docs/`.

## 4. Gates

| # | Check |
|---|---|
| G1 | `git diff --name-only <start HEAD>..HEAD` lists only the allowed paths. `git diff <start HEAD>..HEAD -- . ':(exclude)vignettes/articles' ':(exclude)dev'` is empty (package tarball content unchanged). |
| G2 | Route recorded for every cell. If route b was used: the regenerated data checksum equals the stored one in every replicate used (count). |
| G3 | Reproduction: for OC cells 5 and 10 (already on R(t_j)), the recomputed bias, coverage and `mean_truth_readoff`-based values equal the stored OC summary values within 1e-12 (all months). |
| G4 | Same data, same row: after recomputation, OC cell 8 and M7 cell 31 RMST month 36 agree in `mean_estimate`, `bias`, `bias_mcse`, `coverage_known`, `coverage_estimated` within 1e-12. |
| G5 | Every RMST row: `mean_estimate` equals the mean of the stored per-replicate estimates within 1e-12; `empirical_sd`, mean SEs and mean lengths equal their stored values within 1e-12. |
| G6 | Every non-RMST row of the CSV is byte-identical to before (apart from the added empty column), checked by script with a row count. |
| G7 | Article renders in under 60 s with `options(warn = 2)`; `check_pkgdown()` no problems; `build_article()` completes; Table 2 note and design paragraph state the convention used. |
| G8 | Numerals in changed prose: each inline R from the CSV, a design constant or a formula constant. |

## 5. Report

`dev/audit/REPORT_rmst_readoff_truth.md`:

1. Gate table.
2. Route per cell (and the script lines that define the grid and read-off rule).
3. **One comparative table:** RMST rows, every cell, months 1, 6, 12, 24, 36, 48: bias against R(t) (old) and R(t_j) (new); coverage known and estimated, old and new; mean of R(t) − R(t_j). Then a short plain-language reading.
4. Band truth finding.
5. Article changes, sentence by sentence.
6. Branch state: HEAD, commits since start, `git status`.

After the final commit, copy to `~/Downloads/` and print the paths: the report and the rendered HTML of the article. Do not push.
