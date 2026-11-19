# TASK — weightedsurv: codebase audit and Kaplan–Meier pressure test

**Stage 1 of the CRAN version update** (CRAN currently has 0.1.0, published 2025-12-22).
Written in chat 2026-11-18. Executed by Claude Code in the local `weightedsurv` repo.

| | |
|---|---|
| Touches `R/` | **No.** Nothing that ships is modified. |
| Writes | Only under `dev/tasks/` and `dev/audit/` (`dev/` is in `.Rbuildignore`). |
| Fixes | None in this task. Findings are recorded; fixes are a later task that Larry dispositions by class. |
| Compute | Deterministic test matrix T1–T10: about 4,000 small seeded datasets, expected under 15 min in total. One Monte Carlo step (M1): 6 cells × 2,000 replicates, n = 200, ≤ 5 min wall-clock — run **only if the kickoff says `M1: GO`**. |
| Installs | None. If a suggested package is missing, work without it and record which. |
| Gates | Stop-on-failure, not stop-to-ask. On a failed gate: write what happened at the top of the report, commit, stop. |
| Git | Commit locally. **Never push.** No CRAN, win-builder or R-hub submissions. |

---

## 1. Why this task exists

Larry wants the existing codebase audited, and the Kaplan–Meier estimation functions pressure-tested, before a version update goes to CRAN.

Chat read and ran GitHub HEAD `edd61a2` (2026-02-12) under R 4.3.3 with survival 3.5.8. The local tree may be ahead of that commit. **Treat every chat result below as a lead about your tree, not a fact about it.** Read the current files; re-establish everything by running it.

What chat found on `edd61a2`, in short:

- The core estimator is right. Unweighted survival and Greenwood standard errors equal `survival::survfit` to 1e-16, including with heavy ties. Weighted survival equals `survfit(weights=)` to 1e-16.
- 26 reproducible defects sit around that core: Kaplan–Meier evaluated on grids that omit event times, the weighted variance, quantile ties, input handling, and a few inference details. They are listed in §4 with IDs.
- A reduced `R CMD check` (no vignettes, no manual) was clean. The problems are functional, not check-level. There is no `tests/` directory.

---

## 2. Steps

### Step 0 — Intake (Gate G0)

1. The kickoff has you note the current HEAD as `START_SHA`, then copy this file to `dev/tasks/TASK_km_audit.md` and the seed probe to `dev/audit/km_audit_probe.R`, and commit both.
2. Record in `dev/audit/provenance.txt`: `git rev-parse HEAD`, branch, `git status --porcelain`, R version, `packageVersion("survival")`, `sessionInfo()`.
3. Record whether `edd61a2` is an ancestor of HEAD. If it is, save `git diff --stat edd61a2 HEAD -- R/ NAMESPACE DESCRIPTION vignettes/`.

**G0 fails if** there are uncommitted changes under `R/`, `NAMESPACE` or `DESCRIPTION` at the start, or `pkgload::load_all(".")` errors.

### Step 1 — Seed probes (Gate G1)

```
Rscript dev/audit/km_audit_probe.R > dev/audit/probe_stdout.txt 2>&1
```

The script loads the local tree with `pkgload::load_all`, never the installed CRAN copy. It deliberately does not attach `survival`. It writes `dev/audit/probe_results.csv`.

**G1 fails if** the exit code is non-zero (a `PROBE_ERROR` row), or any control `OK1`–`OK3` is `FAIL`. The controls passed on `edd61a2`; a control failure means the local estimator is not the one chat audited, and the plan needs re-scoping.

`NOT_REPRODUCED` is **not** a failure. It means the local tree differs from `edd61a2` on that point. Record it and find out why (which commit, which lines).

### Step 2 — Source read of the local tree

Read every file in `R/` in full. Then:

1. Build `dev/audit/exports_inventory.csv`, one row per export in `NAMESPACE`: `export, file, line, called_from_R, used_in_vignettes, used_in_examples, role`, where `role` is `user-facing` or `helper`. This feeds a later decision about the exported surface (46 exports on `edd61a2`, several of them small internal helpers).
2. For each static-read lead `S1`–`S10` in §4, confirm or correct it against the current file, with file and line.
3. For each probe-backed finding, record the file and line where the behaviour originates.

### Step 3 — Kaplan–Meier pressure test

Write `dev/audit/km_pressure_test.R` implementing the matrix in §3. Seeded and deterministic. It writes `dev/audit/km_pressure_results.csv`, one row per cell: `family, cell, n_checks, n_fail, worst, note`.

The script exits 0 when tests fail — failures are data. It exits non-zero only when the script itself breaks. Give it a `--no-write` flag that prints the summary without rewriting any file (Gate G2 uses it). Write the T9 timings to a separate `dev/audit/km_scaling.csv`, since they differ run to run.

### Step 4 — M1 Monte Carlo (only if the kickoff says `M1: GO`)

Spec in §3, family M1. If the kickoff says `NO-GO` or says nothing, skip it and write "M1 not run" in the report.

### Step 5 — CRAN-readiness baseline (changes nothing)

1. In a temp directory: `R CMD build .` then `R CMD check --as-cran` on the tarball. Copy `00check.log` to `dev/audit/`. Record the status line, every NOTE and WARNING verbatim, the five slowest examples, and total check time. Do not commit the tarball. If suggested packages are missing, do not install them: set `_R_CHECK_FORCE_SUGGESTS_=false`, and if the vignette cannot build add `--no-build-vignettes` / `--ignore-vignettes`; record exactly which flags you used and which packages were missing.
2. Diff against CRAN 0.1.0: `download.packages("weightedsurv", destdir = <tmp>, type = "source")`, unpack, and `diff -r` its `R/`, `NAMESPACE`, `DESCRIPTION` and `vignettes/` against the tree. Save a summary as `dev/audit/diff_vs_cran_0.1.0.txt`. If offline, record that and continue.
3. Reverse dependencies: `tools::package_dependencies("weightedsurv", reverse = TRUE, which = "all")`. Record the result, or that you were offline.
4. Fill the checklist below, each with yes/no and `file:line` where relevant:

   - `tests/` present; `NEWS.md` present; `cran-comments.md` present
   - `set.seed()` called inside package functions
   - `cat()` / `print()` output not behind a `verbose` argument
   - estimation functions that draw plots as a side effect
   - every `par()` change restored with `on.exit`
   - `requireNamespace()` guards on packages that are in `Imports`
   - functions used in formulas but not imported (see `R6`)
   - every package in `Suggests` actually used, and used conditionally
   - examples over 5 s; use of `\dontrun` / `\donttest`
   - URLs valid (`urlchecker::url_check()` if installed, otherwise skip and say so)

### Step 6 — Report

Write `dev/audit/REPORT_km_audit.md` in the format of §5. Commit.

### Step 7 — Post-conditions (Gate G2)

Run the block in §6 and paste its output as the last section of the report.

---

## 3. Pressure-test matrix

### Comparison rules (read first)

- **Reference for survival and unweighted SE:** `survival::survfit`, evaluated with `summary(fit, times = ..., extend = TRUE)`.
- **References for weighted SE** (report both; neither is pass/fail):
  - `survfit(..., weights = w, robust = TRUE)`;
  - the formula documented in `vignettes/articles/weightedkm_methodology.Rmd`, implemented **in the audit script only**:
    `SE(t) = S(t) * sqrt( sum over s <= t of [ sum_i w_i^2 dN_i(s) ] / Y_w(s)^2 )`.
- **Near-tied times:** the package applies `survival::aeqSurv` before estimating. Build the reference from the same time-fixed vector, or compare at interval midpoints. Evaluating a per-arm `survfit` at raw times that differ from a fixed time by less than 1e-8 produces spurious differences the size of one KM jump. Chat hit this; it is an artefact of the comparison, not a defect.
- **Tolerance:** 1e-10 for quantities that should be identical. Everything else is reported as a ratio or a difference with no pass/fail.

### T1 — Agreement with `survfit` over a design grid

Full factorial, 3 seeds per cell:

| Factor | Levels |
|---|---|
| n per arm | 1, 2, 5, 20, 100, 500 |
| censoring | 0%, ~30%, ~70%, ~95% |
| ties | none (continuous); moderate (rounded to 0.1); heavy (integers 1–8) |
| weights | none; integer 1–3; U(0.3, 3); 1/p with p ~ Beta(2, 8); U(0.3, 3) normalised to sum 1; U(0.3, 3) with 10% set to 0 |
| last observation in each arm | as generated; forced event; forced censored |

For each dataset call `KM_diff()` with default `at_points`, and `df_counting()` with `check.KM = FALSE, get.cox = FALSE`. Check:

- survival, both arms, vs `survfit` — pass/fail at 1e-10;
- unweighted SE vs `survfit` — pass/fail at 1e-10;
- weighted SE — ratio to each of the two references, report only;
- `dhat = surv1 - surv0`, and `sig2_dhat = sig2_surv0 + sig2_surv1` — pass/fail;
- cells where the package errors, or `survfit` errors: record the message.

### T2 — Invariants (no external reference needed)

Each is pass/fail unless marked.

1. Row permutation of `df` changes nothing: survival, SE, quantiles, z-score.
2. Weights `w` vs `c·w` for c in {0.01, 10, 1/sum(w)}: survival identical; SE identical (expected to fail — `K3`; report the ratio).
3. Time rescaling `t → c·t` for c in {1/12, 30.4375}: survival identical at rescaled points; medians scale by c; `KM_diff(modify_tau = TRUE)` on the rescaled data (expected to fail — `K1`).
4. Duplicating every row equals integer weight 2: survival identical.
5. Swapping arm labels: `dhat → −dhat`, z-score → −z-score, quantile rows swap.
6. Adding non-event times to `at_points` leaves survival unchanged at the original points.
7. Removing event times from `at_points` leaves survival unchanged at the retained points (expected to fail — `K1`).
8. `KM_diff(draws = 0)` vs `KM_diff(draws = 2000)` on weighted data: SE ratio within Monte Carlo error; report the ratio.

### T3 — Shape of the curve

For every curve produced in T1: values in [0, 1]; non-increasing; equal to 1 before the first event; carried forward after the last observation; exactly 0 (not ±1e-17) when the last observation is an event (`K6`). Check the same for `survP`, `survG` and the `*_stratum` vectors from `df_counting()`.

### T4 — Quantiles

On 200 seeded datasets spanning the T1 ties and censoring levels, plus the exact-plateau cases (no censoring, n = 2, 4, …, 60):

- `km_quantile()` and `df_counting()$quantile_results` vs `quantile(survfit)` for survival levels 0.9, 0.75, 0.5, 0.25, 0.1;
- state explicitly which convention `qprob` follows — "time at which S drops to `qprob`" or "`qprob`-th percentile of the failure distribution" (`I5`);
- lower and upper limits for `conf_level` in {0.80, 0.90, 0.95}, against `survfit(conf.int = conf_level, conf.type = "log")`;
- `NA` handling when the curve never reaches the level.

### T5 — `KM_diff()` bands and resampling

- `c_alpha_band` is non-increasing in `alpha` over {0.01, 0.05, 0.10, 0.20} (expected to fail — `I1`).
- The simultaneous band contains the pointwise interval at every time.
- With a fixed `seedstart`, two calls return identical output; the caller's `.Random.seed` is unchanged afterwards (expected to fail — `P2`).
- `draws.band` in {1, 5, 19, 20, 21} with default `show_resamples` (`R5`).
- Variance of the columns of `surv0_star` / `surv1_star` across 5,000 draws vs the analytic `sig2_surv0` / `sig2_surv1`, unweighted and weighted: report the ratio along the curve.

### T6 — Internal consistency of `df_counting()`

- `surv0`, `surv1` equal `KM_diff()`'s at shared time points.
- `riskpoints0/1` equal `summary(survfit, times = risk.points)$n.risk`, unweighted and weighted.
- `survP` equals the pooled `survfit`; `survG` equals `survfit` on `1 − status`.
- `z.score^2` equals `survdiff(rho = 0)$chisq` unweighted; `lr` and `sig2_lr` reproduce it.
- `check.KM = TRUE` raises no warning on any clean T1 dataset.

### T7 — Stratified output

With 2, 3 and 6 strata, including a stratum with one arm empty and a stratum with no events: every `*_stratum` vector equals the per-stratum `survfit` at each observation's time. Run once with `survival` attached and once in a clean session without it (`R6`).

### T8 — Input and error-path matrix

For `KM_diff()`, `df_counting()` and `plotKM.band_subgroups()`, feed each input below and classify the outcome as **correct result**, **informative error**, **cryptic error**, or **silent wrong result**. Report as one table: input × function.

empty arm · arm with no events · all censored · one observation · `NA` in time / event / treatment / weight / strata · `Inf` time · zero and negative times · zero / negative / `NA` weights · treatment coded 1/2, as factor, logical, character · event coded 1/2 and logical · unsorted `at_points` · duplicated `at_points` · `at_points` beyond the last time · `at_points` not starting at or before the first event · missing column names · tibble and `data.table` input · `tzero` before the first and after the last time in `wlr_dhat_estimates()`.

### T9 — Scaling

`df_counting()` (defaults), `KM_diff(draws.band = 1000, show_resamples = FALSE)` and `cox_rhogamma(draws = 0)` at n = 500, 2,000, 4,000, 8,000: elapsed time and peak memory (`gc(reset = TRUE)` before, `sum(gc()[, 6])` after). Abort a cell if peak memory passes 6 GB and record that. Chat measured 0.1 s / 150 MB, 0.8 s / 227 MB and 3.0 s / 540 MB at n = 500, 2,000 and 4,000 for `df_counting()`.

### T10 — Plot paths (smoke tests under `pdf(NULL)`)

- `plot_weighted_km()` on 60 weighted datasets whose last observation in each arm is an event (`K6`): count how many draw and how many return `NULL`.
- `KM_plot_2sample_weighted_counting(conf.int = TRUE)` when a curve reaches 0.
- `plotKM.band_subgroups()` with a subgroup that has one arm empty (`R1`), and with `modify_tau = TRUE` (`K1`).
- `cumulative_rmst_bands()` given a fit from `KM_diff(draws.band = 0)` (`S4`).

### M1 — Calibration of the weighted standard error (only if `M1: GO`)

Budget: **6 cells × 2,000 replicates, n = 200, ≤ 5 min.** If it would exceed 10 min, stop and record.

- Data: one arm, T ~ Exp(0.08), C ~ U(0, 30). Truth S(t) = exp(−0.08 t).
- Cells: ties {continuous; times rounded up to 0.5} × weights {none; w ~ U(0.3, 3) independent of T; the same w normalised to sum 1}.
- Evaluate at t = 5, 10, 15.
- Estimators of SE: the package (`calculate_risk_event_counts` + `KM_estimates`); `survfit` robust; the documented formula from the comparison rules.
- Report, per cell × estimator × time: mean(SE) / empirical SD of the estimate, and coverage of `S ± 1.96·SE`.

---

## 4. Finding register from chat (re-establish each on the local tree)

`chat` = result on `edd61a2`. **Class** is the kind of `R/` change a fix would be: **B** changes behaviour (output changes only where it is currently wrong or crashing), **M** changes the method, **—** record only. The class column is chat's suggestion for Larry's later disposition. It authorises nothing.

### Probe-backed

| ID | Finding | chat (`edd61a2`) | Class |
|---|---|---|---|
| OK1 | Unweighted S and Greenwood SE equal `survfit` | max diff 3e-16 | — |
| OK2 | Weighted S equals `survfit(weights=)` | max diff 3e-16 | — |
| OK3 | Weighted SE, natural-scale weights, no ties, vs `survfit` robust | ratio 0.92–1.00 where S > 0.1 | — |
| K1 | `KM_diff(modify_tau = TRUE)` computes KM on a grid without the event times. `calculate_risk_event_counts()` / `KM_estimates()` are only Kaplan–Meier when `at_points` contains every event time | n = 60: error up to 0.43, survival −0.31. n = 200: 0.023. n = 1,000: 0.011. n = 200 with time in years: 0.17 | B |
| K1b | Same cause through a user-supplied `at_points` | `seq(0, 30, 5)`: error 0.23 | B |
| K2 | Median when S equals 0.5 on an interval: left end returned, `survfit` returns the midpoint; which one you get depends on floating-point rounding | 18 of 30 cases differ | B |
| K2b | Through `df_counting()` defaults the mismatch is silent (the quantile cross-check only signals when `verbose = TRUE`) | 5.0 vs 5.5, no warning | B |
| K3 | Weighted SE changes when weights are rescaled | SE(w/Σw) / SE(w) = 0.18, 0.64, 1.60, 6.31 along the curve | M |
| K4 | Weighted SE with tied event times is too large | 2.06×, 1.77×, 1.52×, 1.25× `survfit` robust | M |
| K5 | All weights exactly 2 uses the equal-weight branch; nearly-2 uses the weighted branch | SE 0.0297 vs 0.0416 | M |
| K6 | Weighted KM ending on an event leaves ±1e-17 instead of 0; negative values make `check_km_curve()` stop | negative in 181 of 400; plot function stops in 41 of 60 | B |
| K7 | Unsorted `at_points` accepted | S from −166 to −38, no error | B |
| R1 | Empty arm in `KM_diff()` | `object 'sig2w_multiplier' not found` | B |
| R2 | Arm with no events in `df_counting()` | `subscript out of bounds` | B |
| R3a | `NA` time in `KM_diff()`: warns, then fails | `missing value where TRUE/FALSE needed` | B |
| R3b | `NA` time in `df_counting()` | `missing value where TRUE/FALSE needed` | B |
| R3c | `NA` weight in `KM_diff()` | survival all `NA`, no warning | B |
| R4a | Negative weights in `KM_diff()` | accepted silently | B |
| R4b | Zero weights in `df_counting()`: KM is right, the self-check compares vectors of different length | false "Discrepancy in KM curve fit" | B |
| R5 | `KM_diff(draws.band < 20)` with default `show_resamples = TRUE` | `subscript out of bounds` | B |
| R6 | Stratified `df_counting()` without `library(survival)`: `strata()` is used in a formula but not imported | `could not find function "strata"` | B |
| R7 | Event coded 1/2: `df_counting()` accepts, `KM_diff()` rejects | inconsistent | B |
| I1 | Band critical value hard-coded at the 0.95 quantile; `alpha` ignored | 2.9796 for α = 0.05 and 0.20 | B |
| I2a | `wlr_dhat_estimates()`: pooled S read one grid point after `tzero` | variance ratio 0.947 | B |
| I2b | `wlr_dhat_estimates()` with `tzero` past the last time | `dhat` of length 0, no error | B |
| I3 | `cumulative_rmst_bands()` integrates the KM step function by the trapezoid rule; the RMST vignette documents rectangles | 4.2836 vs exact 4.2688 | B — every RMST value moves; needs Larry's explicit call |
| I4 | Quantile cross-check always uses `survfit`'s 0.95 limits | false discrepancy at `conf_level = 0.90` | B |
| I5 | `qprob ≠ 0.5`: `quantile_results` is the time S drops to `qprob`; `quantile_check` asks `survfit` for the `qprob`-th percentile | 15.54 vs 4.29 at `qprob = 0.25` | B |
| P2 | `set.seed()` inside functions overwrites the caller's RNG state | confirmed | B |
| P5 | O(n²) matrices in `df_counting()` | n = 4,000: 3.0 s, 540 MB | — |

### Static-read leads (no probe yet — confirm or correct against the current file)

| ID | Lead |
|---|---|
| S1 | `wt.rg.S(scheme = "MB")`: when `mb_tstar` equals the last time, `which.max()` on an all-`FALSE` vector returns 1 and every weight becomes 1; when `mb_tstar` is beyond it, weights are `1/S_left`, infinite where `S_left = 0` |
| S2 | `wlr_dhat_estimates()` and `wlr_cumulative()`: `requires_cols` is a typo for `required_cols`, so the `survG` requirement for Schemper weights is never enforced |
| S3 | `KM_diff()` default `at_points = sort(df[[tte.name]])` uses raw times while estimation uses `aeqSurv`-fixed times |
| S4 | `cumulative_rmst_bands()`: 1.96 and the 0.95 quantile are hard-coded; a fit from `KM_diff(draws.band = 0)` has `dhat_star = NULL` and the function fails inside `apply()` |
| S5 | `plot_km_confint_polygon()` and `km_quantile()`: log-transform limits are not truncated at 1 and are `NaN` where S = 0 |
| S6 | `plotKM.band_subgroups(modify_tau = TRUE)` reuses the grid for every subgroup `KM_diff()` call, and `cumulative_rmst_bands()` reuses `fit$at_points`, so `K1` propagates to both |
| S7 | `get_validated_weights()` writes to the console with `cat()` when `mb_tstar` is missing |
| S8 | `df_counting()` does not validate `strata.name` (a missing column gives `object '…' not found`) and does not validate the event coding |
| S9 | `calculate_risk_event_counts(draws > 0)` reseeds with the same `seedstart` for both arms and every stratum |
| S10 | Housekeeping: an inner `check_km_curve()` in `df_counting()` shadows the exported one with a different signature; a branch in `wlr_dhat_estimates()` whose condition can never be true; `ans$strata` assigned twice; `pkgdown` in `Suggests` |

---

## 5. Report format — `dev/audit/REPORT_km_audit.md`

Tables first, then a short plain reading. No dense narrative.

- **Top line.** If a gate failed, say which and why, here.
- **A. Provenance.** SHA, branch, R and survival versions, relation to `edd61a2`, one-paragraph summary of the diff against CRAN 0.1.0.
- **B. Finding register.** One table, one row per ID in §4 plus any new ones (`N1`, `N2`, …):
  `| ID | area | finding | chat (edd61a2) | local tree: status and key number | file:line | class |`
  Status is one of `REPRODUCED`, `NOT REPRODUCED`, `CHANGED`, `CONFIRMED` / `CORRECTED` (for S-leads). Then at most ten lines of plain-language reading.
- **C. Pressure test.** One table, rows T1–T10: cells run, pass, fail, worst discrepancy. Then a list of failures that are not already in the register.
- **D. M1.** One table, cells × estimators: SE/SD ratio and coverage at t = 5, 10, 15. Three lines of reading. Or "M1 not run".
- **E. CRAN baseline.** Check status line, NOTES/WARNINGS verbatim, the checklist from Step 5.
- **F. Exports.** Counts by role, and the helpers not used outside `R/`.
- **G. Post-conditions.** Output of the §6 block.

Record findings only. **No recommendations section and no follow-up task proposals.** If something blocked the audit itself, say so in the top line.

---

## 6. Gate G2 — post-conditions

Run from the repo root with `START_SHA` set. Every check must print `OK`.

```bash
# 1. nothing that ships was touched
git diff --name-only "$START_SHA" HEAD | grep -Ev '^dev/(tasks|audit)/' && echo "FAIL: paths outside dev/tasks, dev/audit" || echo "OK 1"

# 2. everything committed
[ -z "$(git status --porcelain)" ] && echo "OK 2" || echo "FAIL: uncommitted changes"

# 3. deliverables exist
for f in tasks/TASK_km_audit.md audit/km_audit_probe.R audit/probe_results.csv audit/probe_stdout.txt \
         audit/provenance.txt audit/exports_inventory.csv audit/km_pressure_test.R \
         audit/km_pressure_results.csv audit/km_scaling.csv audit/00check.log audit/REPORT_km_audit.md; do
  [ -f "dev/$f" ] || echo "FAIL: missing dev/$f"; done; echo "OK 3 (if no FAIL above)"

# 4. every register ID appears in the report
for id in OK1 OK2 OK3 K1 K1b K2 K2b K3 K4 K5 K6 K7 R1 R2 R3a R3b R3c R4a R4b R5 R6 R7 \
          I1 I2a I2b I3 I4 I5 P2 P5 S1 S2 S3 S4 S5 S6 S7 S8 S9 S10; do
  grep -q "| *$id *|" dev/audit/REPORT_km_audit.md || echo "FAIL: $id not in report"; done; echo "OK 4 (if no FAIL above)"

# 5. scripts re-run cleanly, without rewriting their outputs
Rscript dev/audit/km_audit_probe.R   --no-write > /dev/null 2>&1 && echo "OK 5a" || echo "FAIL: probe"
Rscript dev/audit/km_pressure_test.R --no-write > /dev/null 2>&1 && echo "OK 5b" || echo "FAIL: pressure test"
[ -z "$(git status --porcelain)" ] && echo "OK 5c" || echo "FAIL: re-running left the tree dirty"

# 6. nothing pushed (informational)
git log --oneline "$START_SHA"..HEAD
```

Apart from the timing rows (probe `P5`, and `km_scaling.csv`), two runs of either script must give identical results. If they do not, fix the seeding in the audit script — not in `R/` — re-run and re-commit.
