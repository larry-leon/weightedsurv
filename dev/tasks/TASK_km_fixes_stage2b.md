# TASK — weightedsurv 0.2.0: Stage 2b on branch `km-fixes-0.2.0`

Written in chat 2026-11-19 from `dev/audit/REPORT_km_fixes_stage2.md` and `dev/audit/REPORT_km_m1b.md`.
Executed by Claude Code in the local `weightedsurv` repo.

| | |
|---|---|
| Touches `R/` | **Yes.** Five items, listed by class in §3. |
| Branch | `km-fixes-0.2.0` (exists). Every commit goes there. No merge. **Never push.** |
| Authorised — **changes the method** | **M**: weighted Kaplan–Meier variance becomes weighted Greenwood. Approved by Larry after M1b. |
| Authorised — changes behaviour | **B8** (`mean_rule` option for RMST), **B9** (`show_resamples` default), **B11** (semi-transparent bands) |
| Authorised — housekeeping | **T0**: testthat floor in `DESCRIPTION` |
| **Not authorised** | Anything else. In particular: the resampling variance path (`draws > 0`), the RMST default rule, the RMST draws / variance / sup statistic, the hard-coded RMST levels, and `plotKM.band_subgroups()`. |
| Compute | **M1c**: 6 cells × 2,000 replicates, about 20 s — only if the kickoff says `M1c: GO`. Test and regression suites: minutes. forestsearch check: about 1 min. |
| Installs | Temporary library only. |
| Gates | Stop-on-failure, not stop-to-ask. Same package-gate rule as Stage 2: a package that cannot meet its contract is restored and recorded as NOT DONE; the task continues. |

---

## 1. Ground rules

1. **Read the current file on the branch before editing it.** This document was written without sight of the Stage-2 code; where it names a function introduced in Stage 2, confirm it exists.
2. One commit per item, message prefixed with its ID. Tests go in the same commit.
3. **Exported API: exactly two changes are allowed**, both specified below — `cumulative_rmst_bands()` gains a last argument `mean_rule`; `KM_diff()`'s `show_resamples` default becomes `FALSE`. No other export, argument, default or returned-element name changes. Returned lists may gain elements.
4. **Results do not move** except where an item says so. The baseline is now the branch as it stands at the start of this task, not 0.1.0.
5. `vignettes/`: only the two article edits named in B8 and M. Nothing else. List, without editing, any other sentence in `vignettes/` or `README.md` that these changes make inaccurate.
6. Stage-1 and Stage-2 files under `dev/audit/` are a record. New material goes in `dev/audit/stage2b/`. Copy a script before changing it.
7. Error and warning messages you add must not contain `TRUE/FALSE`, `sig2w_multiplier` or `subscript out of bounds`.
8. If an existing test encodes behaviour that an item here changes by design, update that test in the item's commit and list it in the report. Do not weaken any other test.

---

## 2. Steps

### Step 0 — Intake (Gate G0)

1. Confirm the branch is `km-fixes-0.2.0` and the tree is clean. Note HEAD as `START_SHA` and `git rev-parse master` as `MASTER_SHA`.
2. `mkdir -p dev/audit/stage2b`; copy this file to `dev/tasks/TASK_km_fixes_stage2b.md`; commit.

**G0 fails (stop) if** the branch is wrong, the tree is not clean, `dev/audit/REPORT_km_fixes_stage2.md` is missing, or `pkgload::load_all(".")` errors.

### Step 1 — Baseline, before any `R/` edit (Gate G1)

Using copies of the Stage-2 scripts in `dev/audit/stage2b/`:

- `golden.R make` → `golden_s2.rds`: the same datasets and calls as Stage 2, on the current branch code. Add these calls so the new items are covered:
  - `df_counting()` and `KM_diff()` with non-constant weights and `draws = 0` on every dataset (use a seeded U(0.3, 3) weight column where a dataset has none);
  - the same weighted `KM_diff()` with `draws.band = 200, show_resamples = FALSE`;
  - `cumulative_rmst_bands(plot = FALSE)` on an unweighted and a weighted `draws.band = 200` fit.
- Record, for every golden call, the set of warnings it raises.
- `api.R make` → `api_s2.rds`.

**G1 fails (stop) if** `git diff --quiet "$START_SHA" -- R/ NAMESPACE` is false at this point, or a snapshot call errors.

### Step 2 — Items, in this order

`T0 → B9 → B11 → B8 → M`

For each: implement to the contract in §3, add its tests, run the tests, run `golden.R compare` with the cumulative whitelist, then commit.

### Step 3 — Documentation and NEWS (one commit)

- `NEWS.md`, under 0.2.0: one bullet each for M, B8, B9, B11, in user language. M goes under "Changes in behaviour" and says plainly that every weighted standard error changes and unweighted results do not.
- Help pages as each contract specifies.
- The two article edits named in B8 and M.

### Step 4 — Verification (Gate G2, stop on failure)

Run the block in §5. In addition:

1. **Pressure test.** Re-run the Stage-2 copy of the pressure test into `dev/audit/stage2b/`. Tabulate Stage 2 → Stage 2b by family. Expected: the three T2 "weights w → c·w" cells now pass. Expected new failure, by design: the SE part of T2.4 (duplicated rows versus weight 2), because weights are no longer frequency weights. Report the T1 weighted-SE ratios to `survfit` robust, before and after, by ties × weights. Every other change must be explained.
2. **M1c** (only if `M1c: GO`). Run a copy of the M1b script that differs only in its output path, so that its "package" rows are produced by the new package code. Write `dev/audit/stage2b/km_m1c_results.csv`. Required: the "package" rows equal the "weighted Greenwood" rows in every cell (SE/SD and coverage to the printed precision).
3. **`R CMD check --as-cran`** on a tarball built from a clean clone: 0 ERROR, 0 WARNING; quote every NOTE. Save as `dev/audit/stage2b/00check_2b.log`.
4. **forestsearch**: the same check as in Stage 2, against the new build.
5. **Plot sheets**: run the Stage-2 plot comparison in a third mode, writing `dev/audit/stage2b/km_plots_after2.pdf` and `km_plot_compare_after2.csv`. The numbers must equal the Stage-2 "after" file (geometry is unchanged by B11).

### Step 5 — Report

`dev/audit/REPORT_km_fixes_stage2b.md`, format in §4. Commit. Then run the reporting line from the kickoff.

---

## 3. Items — contracts

### T0 — testthat floor. Housekeeping
`DESCRIPTION`: `testthat (>= 3.1.7)`. Nothing else.

### B9 — `show_resamples` default. Class: changes behaviour
- `KM_diff(show_resamples = FALSE)` becomes the default. Nothing else in `KM_diff()` changes.
- If the roxygen example of `KM_diff()` or `README.md` relies on the old default to show the resample plot, add `show_resamples = TRUE` to that call so what is displayed stays the same.
- Golden whitelist: none.

### B11 — semi-transparent confidence bands. Class: changes behaviour (drawing only)
- In the KM-curve plot (`plot_km_curves_counting()` and its callers), the two bands are filled with semi-transparent colour so each stays visible where they overlap. Keep the present hues (grey for control, light blue for treated); use `grDevices::adjustcolor(<hue>, alpha.f = 0.5)`. Curves and censoring marks are drawn on top, as now.
- No new argument. Geometry is unchanged: the vertex helper from B10 returns what it returned before.
- `grDevices` may be added to `Imports` and `NAMESPACE` (base R).
- `plotKM.band_subgroups()` and the RMST plot are not touched.
- Test: the fill colours passed to `polygon()` have alpha below 1; the vertices are identical to before.
- Golden whitelist: none.

### B8 — `mean_rule` option for RMST. Class: changes behaviour (adds an option; nothing existing moves)
- `cumulative_rmst_bands(..., mean_rule = c("trapezoid", "rectangle"))`: a new **last** argument, resolved with `match.arg()`. Default `"trapezoid"`.
- With the default, every element of the result is identical to today's.
- `"rectangle"` is the left-rectangle cumulative sum: `R(t_j) = sum over k <= j of d(t_{k-1}) * (t_k - t_{k-1})`, with `R = 0` at the first grid point.
- The chosen rule is applied through **one internal helper** to the estimate and to both sets of resampled curves. Nothing else changes: not the draws, the seeds, the variance from the first draw set, the standardisation, the sup statistic, the 0.95 quantile or the 1.96.
- The returned list gains `mean_rule`.
- Required, with `mean_rule = "rectangle"` on a fit whose grid holds every event time (the default `KM_diff()` grid, and `fit_itt` from `plotKM.band_subgroups()`): `rmst_time` at the last grid point, and at three interior grid points, equals the difference of `survfit`'s restricted means (`summary(fit, rmean = tau)`, treated minus control) within 1e-10. Test this on veteran, GBSG and one simulated dataset.
- Help page: describe both rules; say that `"rectangle"` equals `survfit` when the time grid holds every event time, which is the default, and not when the fit was made with `modify_tau = TRUE`.
- Article: in `vignettes/articles/weightedrmst_methodology.Rmd`, correct the description of the integration rule to match (trapezoid by default, rectangle optional). Edit only that passage.
- Golden whitelist: none. The new element `mean_rule` is an addition.

### M — weighted Kaplan–Meier variance. Class: **changes the method**

**Approved estimator** (weighted Greenwood). For a group with weights `w_i`, weighted at-risk total `Y_w(s)`, weighted events `dN_w(s) = sum_i w_i dN_i(s)`:

`Var[S(t)] = S(t)^2 * sum over s <= t of  [ sum_i w_i^2 dN_i(s) ] / [ Y_w(s) * (Y_w(s) - dN_w(s)) ]`

A term is 0 where `Y_w(s) - dN_w(s) <= 1e-10 * Y_w(s)` (the B2 tolerance) or `Y_w(s) = 0`.

**Where it goes.** The analytic path of `calculate_risk_event_counts()` — the `sig2w_multiplier` it returns when no resampling is asked for — and nothing else.

- One formula for all weights. The separate equal-weights formula goes away: at unit weights the expression above *is* Greenwood, so unweighted results do not move.
- The resampling path is untouched. It is used exactly when it is used today: `draws > 0` **and** the weights are not all equal. In every other case the formula above is used. Keep that switch as it is.
- `KM_estimates()` is unchanged, including its default when no multiplier is supplied.
- Signatures and returned names are unchanged.

**What must hold afterwards**

| Property | Requirement |
|---|---|
| Unit weights | Identical to today (Greenwood; equals `survfit`'s `std.err`) within 1e-12. |
| Scale | SE for `w`, `10·w`, `0.01·w` and `w / sum(w)` identical within 1e-9 relative. |
| Constant weights | All weights equal to any `c > 0` give the unit-weight SE. |
| Formula | Equals the expression above computed independently in the test file, with and without tied event times, within 1e-9 relative. |
| Survival estimate | Unchanged everywhere. |
| Resampled curves and resampling variance | Unchanged everywhere. |
| Warnings | No golden call raises a warning it did not raise in `golden_s2.rds`. |

**Golden whitelist.** A field may move only in a call that passes non-constant weights with `draws = 0`, and only if it is a variance or is computed from one: `sig2_*`, standard errors, pointwise limits, simultaneous-band limits, standardised resamples, band critical values, and weighted quantile confidence limits. Must not move anywhere: survival estimates, `dhat`, resampled curves, counts, risk tables, quantile point estimates, log-rank and Cox output, RMST output, and every field of every unweighted call. For each `sig2_surv0` / `sig2_surv1` that moves, `golden.R` checks the new value against the formula above.

**Documentation**
- On the help pages of `calculate_risk_event_counts()`, `KM_diff()` and `df_counting()`, state: the formula; that weights are treated as fixed and known, with no adjustment for estimated weights; that the SE does not depend on the scale of the weights; and that weights are therefore not frequency weights — a weight of 2 is not two copies of the row.
- Article: in `vignettes/articles/weightedkm_methodology.Rmd`, correct the variance formula to the expression above. Edit only that passage.

**Tests.** The properties in the table. No test may compare the weighted SE with `survfit`.

---

## 4. Report format — `dev/audit/REPORT_km_fixes_stage2b.md`

Tables first, then a short plain reading.

- **Top line.** Branch, `START_SHA`, final SHA, gates, any item NOT DONE with the reason.
- **A. Items.** `| ID | class | status | commit | files touched | golden fields that moved | tests added |`
- **B. M in numbers.** One table, rows = golden datasets with weights: weighted SE at the 25th / 50th / 75th percentile of the time grid, before → after, and the ratio to `survfit` robust before → after. Then the M1c table if run (package rows beside weighted Greenwood rows). Three lines of reading.
- **C. B8.** One table, rows = veteran, GBSG, simulated: RMST at the last time under `"trapezoid"`, under `"rectangle"`, and from `survfit`; and the band critical value under each rule.
- **D. Seed probe and pressure test**, Stage 2 → Stage 2b, with each remaining failure attributed.
- **E. API.** The printed list of differences from `api_s2.rds`.
- **F. CRAN check and forestsearch.** Status lines verbatim; test count and runtime.
- **G. Sentences in `vignettes/` or `README.md`** that the changes make inaccurate and that this task did not edit, with file and line.
- **H. Gate G2 output.**

Findings only. No follow-up task proposals.

---

## 5. Gate G2 — post-conditions

Run from the repo root with `START_SHA` and `MASTER_SHA` set, after the last item commit. Every line must print `OK`. Checks 5 and 6 write two files: paste the output into the report, commit the report and those files, then run check 9.

```bash
S=dev/audit/stage2b

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. exports unchanged; only the two allowed API changes
diff <(git show "$START_SHA":NAMESPACE | grep '^export' | sort) <(grep '^export' NAMESPACE | sort) > /dev/null && echo "OK 2a" || echo "FAIL: exports changed"
Rscript -e 'pkgload::load_all(".", quiet = TRUE); f <- formals(cumulative_rmst_bands)
  stopifnot(identical(tail(names(f), 1), "mean_rule"), identical(eval(f$mean_rule), c("trapezoid", "rectangle")),
            isFALSE(formals(KM_diff)$show_resamples)); cat("OK 2b\n")'
Rscript $S/api.R compare && echo "OK 2c" || echo "FAIL: API differs beyond the two allowed changes"

# 3. golden: identical outside the Stage-2b whitelist; contract checks pass
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. M: scale-free, constant weights = unit weights, equals the approved formula (ties included)
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(1); n <- 150; U <- ceiling(rexp(n, .1) * 2) / 2; D <- rbinom(n, 1, .7); w <- runif(n, .3, 3); at <- sort(unique(U))
  f <- function(w) calculate_risk_event_counts(U, D, w, at)$sig2w_multiplier
  stopifnot(isTRUE(all.equal(f(w), f(w / sum(w)), tolerance = 1e-9)), isTRUE(all.equal(f(w), f(10 * w), tolerance = 1e-9)),
            isTRUE(all.equal(f(rep(2, n)), f(rep(1, n)), tolerance = 1e-12)))
  Yw <- colSums(outer(U, at, ">=") * w); N1 <- colSums(outer(U[D == 1], at, "<=") * w[D == 1]); N2 <- colSums(outer(U[D == 1], at, "<=") * w[D == 1]^2)
  d1 <- diff(c(0, N1)); d2 <- diff(c(0, N2)); ref <- ifelse(Yw - d1 > 1e-10 * Yw, d2 / (Yw * (Yw - d1)), 0)
  stopifnot(isTRUE(all.equal(f(w), ref, tolerance = 1e-9)))
  Y1 <- colSums(outer(U, at, ">=")); dn <- diff(c(0, colSums(outer(U[D == 1], at, "<="))))
  stopifnot(isTRUE(all.equal(f(rep(1, n)), ifelse(Y1 > dn, dn / (Y1 * (Y1 - dn)), 0), tolerance = 1e-12))); cat("OK 4\n")'

# 5. seed probe: the weighted-variance findings are gone; only the default RMST rule still differs from survfit
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2b.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2b.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 5a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2b.txt)" = "3" ] && echo "OK 5b" || echo "FAIL: controls"

# 6. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2b.txt 2>&1 && echo "OK 6a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2b.log | grep -Evq 'ERROR|WARNING' && echo "OK 6b" || echo "FAIL: R CMD check"
grep -q 'testthat (>= 3.1.7)' DESCRIPTION && echo "OK 6c" || echo "FAIL: testthat floor"

# 7. vignettes: only the two named articles may differ
git diff --name-only "$START_SHA" HEAD -- vignettes/ | grep -Ev '^vignettes/articles/(weightedkm_methodology|weightedrmst_methodology)\.Rmd$' && echo "FAIL: other vignette files changed" || echo "OK 7"

# 8. deliverables and record
for f in golden_s2.rds api_s2.rds probe_2b.txt 00check_2b.log km_plots_after2.pdf km_plot_compare_after2.csv; do
  [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2b.md ] && echo "OK 8 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD
```

```bash
# 9. after the final commit
[ -z "$(git status --porcelain)" ] && echo "OK 9" || echo "FAIL: uncommitted changes"
```

Notes on the block:
- Check 2c: adapt the copied `api.R` so that it exits 0 when the differences from `api_s2.rds` are exactly the two allowed ones, and prints the list either way.
- Check 5a: finding `I3` stays `REPRODUCED` because the probe uses the default rule, which is the trapezoid by Larry's decision. It is a documented default, not a defect.
- If an item is NOT DONE, the checks tied to it fail; report them as failed for that reason and name the item.
