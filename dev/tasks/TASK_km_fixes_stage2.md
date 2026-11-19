# TASK — weightedsurv 0.2.0: Stage-2 fixes from the KM audit, on a branch

**Stage 2 of the CRAN update.** Written in chat 2026-11-19 from `dev/audit/REPORT_km_audit.md` (HEAD `166919e`).
Executed by Claude Code in the local `weightedsurv` repo.

| | |
|---|---|
| Touches `R/` | **Yes.** Every change is listed by class in §3. |
| Branch | `km-fixes-0.2.0`, created from `master`. Every commit goes there. No merge. **Never push.** |
| Authorised — changes behaviour | **B1, B2, B3, B4, B5, B6, B7**, and **B10** if the kickoff says `B10: GO`. |
| Authorised — moves existing code | **V** |
| **Not authorised** | **B8** (RMST by rectangles), **B9** (`show_resamples` default), **M** (weighted KM variance). Do not implement them, even partially. Do not change the weighted variance formula or the trapezoid rule. |
| Compute | M1b: 6 cells × 2,000 replicates, about 20 s — only if the kickoff says `M1b: GO`. Test and regression suites: minutes. forestsearch reverse-dependency check: time-boxed at 30 min. |
| Installs | Into a temporary library only (the new build, and what the forestsearch check needs). Nothing into the user library. |
| Gates | Stop-on-failure, not stop-to-ask. See §2 for what stops the task and what only drops one package. |

---

## 1. Ground rules

1. **Read the current file before editing it.** The `file:line` references in the audit report are for `edd61a2`; confirm them.
2. **One commit per work package**, message prefixed with its ID (`B7: import survival::strata`). Tests for a package go in the same commit.
3. **The exported API is frozen.**
   - No export added or removed.
   - On every exported function: no argument removed, renamed, reordered or given a new default. An argument may be appended at the end, with a default that reproduces today's behaviour, only if a contract below cannot be met otherwise; list each one in the report.
   - No element removed or renamed in any returned list. Elements may be added.
   - New helpers are internal: not exported, `@noRd` or `@keywords internal`.
4. **Results on valid input do not move** unless a package below says they do. This is enforced by the golden snapshot (Step 2).
5. **No changes to `vignettes/`.**
6. `devtools::document()` is allowed. Afterwards the set of `export()` lines in `NAMESPACE` must be unchanged.
7. Stage-1 files in `dev/audit/` are a record. Do not rewrite them. New material goes in `dev/audit/stage2/`. Run the Stage-1 scripts with `--no-write`.
8. Error messages you add must not contain the strings `TRUE/FALSE`, `sig2w_multiplier` or `subscript out of bounds` (the seed probe keys on them).

---

## 2. Steps

### Step 0 — Intake and branch (Gate G0)

1. Note HEAD as `START_SHA`. Confirm the tree is clean and the branch is `master`.
2. `git switch -c km-fixes-0.2.0`
3. `mkdir -p dev/audit/stage2`; copy this file to `dev/tasks/TASK_km_fixes_stage2.md`; commit.

**G0 fails (stop) if** the tree is not clean, `dev/audit/REPORT_km_audit.md` is missing, or `pkgload::load_all(".")` errors.

### Step 1 — M1b (only if `M1b: GO`; touches `dev/audit/` only)

Skip this step if `dev/audit/REPORT_km_m1b.md` already exists.

In a copy of `dev/audit/km_m1_montecarlo.R` (save as `dev/audit/stage2/km_m1b_montecarlo.R`) add a fourth estimator, "weighted Greenwood":

`SE(t) = S(t) * sqrt( sum over s <= t of [ sum_i w_i^2 dN_i(s) ] / [ Y_w(s) * (Y_w(s) - dN_w(s)) ] )`, a term being 0 where `Y_w(s) - dN_w(s) <= 0`.

Same design and seeds as M1. Write `dev/audit/stage2/km_m1b_results.csv` and `dev/audit/REPORT_km_m1b.md`: the Stage-1 section-D table with the fourth estimator added, then three lines of reading. Findings only.

### Step 2 — Baseline capture, before any `R/` edit (Gate G1)

Write and run, against the unmodified code:

**a. Golden snapshot** — `dev/audit/stage2/golden.R make` → `golden_before.rds`.

Datasets (all seeded):
- `survival::veteran` with `treat = trt - 1`;
- `survival::gbsg` with `tte = rfstime/30.4375`, `event = status`, `treat = hormon`; and its subsets `er <= 0`, `er > 0`;
- simulated: continuous unweighted n = 200; ties (times rounded to 0.5) n = 200; weighted (U(0.3, 3), no ties) n = 300; three strata n = 300.

Calls (store every numeric, character and logical element of each result, and the element names):
- `df_counting()` — defaults; with `weight.name`; with `strata.name`; with `scheme = "MB"` and `mb_tstar` inside the follow-up range; with `qprob = 0.25`;
- `KM_diff()` — defaults; `draws.band = 200, show_resamples = FALSE`; the same with `modify_tau = TRUE`; weighted with `draws = 200`;
- `cox_rhogamma()` — FH(0,0) and FH(0,1), `draws = 0` and `draws = 200`;
- `wlr_dhat_estimates()` — `tzero` at the median observed time;
- `cumulative_rmst_bands(plot = FALSE)` on the `draws.band = 200` fit;
- `km_quantile()` on each arm of each dataset.

**b. API snapshot** — `dev/audit/stage2/api.R make` → `api_before.rds`: for every export, the argument names in order and each default, deparsed.

**c. Plot comparison, "before"** — Step P below, mode `before`.

`golden.R compare <whitelist>` and `api.R compare` must exist too (`compare final` uses the union of the whitelists of the packages that were done): they exit 0 when everything outside the whitelist is identical (numeric tolerance 1e-12 absolute; NA pattern, names, lengths and strings exact) and print each difference otherwise.

**G1 fails (stop) if** `git diff --quiet "$START_SHA" -- R/ NAMESPACE` is false at this point, or any snapshot call errors on these datasets.

### Step 3 — Work packages, in this order

`B7 → B6 → V → B4 → B2 → B1 → B3 → B5 → B10`

For each: implement to the contract in §3, add its tests, run `devtools::test()`, run `golden.R compare` with that package's whitelist and `api.R compare`, then commit.

**Package gate.** If a package cannot be made to meet its contract with the golden and API comparisons clean, do **not** commit it: `git restore` the tree to the last commit, record the package as NOT DONE with the reason, and go on to the next. Packages are independent, except that B10 needs B2.

### Step 4 — Package infrastructure (one commit)

- `tests/testthat/` (edition 3) holding the tests written in Step 3, plus `tests/testthat.R`.
- `DESCRIPTION`: `Version: 0.2.0`; add `testthat (>= 3.0.0)` to Suggests and `Config/testthat/edition: 3`; remove `pkgdown` from Suggests. No new Imports.
- `NEWS.md`: `# weightedsurv 0.2.0`, one bullet per package in user language, under "Bug fixes" and "Changes in behaviour". State what moves for users.
- Test suite runtime under 60 s. Tests need only `survival` and `testthat`. No test asserts a weighted-SE value (M is pending).

### Step 5 — Verification (Gate G2, stop on failure)

Run the block in §5. In addition:

1. **Pressure test.** Re-run `dev/audit/km_pressure_test.R` writing to `dev/audit/stage2/` (add an output-directory argument to a *copy* of the script if it has none). Tabulate before and after by family T1–T10. Every remaining failure must be attributed to a held item (M, B8, B9) or explained individually.
2. **`R CMD check --as-cran`** on the built tarball: 0 ERROR, 0 WARNING. Quote every NOTE verbatim; the expected ones are the two local-tooling notes seen in Stage 1 and any informational CRAN-incoming note. Save as `dev/audit/stage2/00check_after.log`.
3. **forestsearch.** In a temp library install the new build. Download the forestsearch source from CRAN and run `R CMD check --no-manual --ignore-vignettes` on it with that library first on the library path (time box 30 min). If it shows an ERROR or a test failure, repeat with CRAN weightedsurv 0.1.0 to see whether it is new. A **new** failure fails G2. If the check cannot be run, say exactly why; that does not fail G2.
4. **Plot comparison, "after"** — Step P, mode `after`.

### Step P — Kaplan–Meier plot comparison against `survfit()` on GBSG

Write `dev/audit/stage2/km_plot_compare.R <before|after>`. Run it in mode `before` during Step 2 and in mode `after` during Step 5. This step runs whatever the kickoff says about B10.

Data: `survival::gbsg`, `tte = rfstime/30.4375`, `event = status`, `treat = hormon`.

Subgroups (R expressions on the GBSG columns; estrogen first):

```
er <= 0 | er > 0 | er <= 5 | er > 5 | er <= 20 | er > 20 | er <= 50 | er > 50 | er > 100
er <= 0 & meno == 0 | er <= 0 & meno == 1 | er <= 0 & grade == 3 | er <= 0 & nodes > 3
er <= 0 & pgr <= 0 | er <= 0 & age <= 50 | er > 0 & pgr > 0 | er > 0 & grade == 3 | er > 0 & age <= 40
pgr <= 0 | pgr > 0 | pgr <= 20 | pgr > 20
meno == 0 | meno == 1 | grade == 1 | grade == 2 | grade == 3
nodes <= 3 | nodes > 3 | nodes > 9 | size <= 20 | size > 20 | size > 50
age <= 40 | age > 40 & age <= 60 | age > 60 | (all patients)
```

**KM-curve sheets** — `dev/audit/stage2/km_plots_<mode>.pdf`, one page per subgroup, four panels:
1. `plot_weighted_km(df_counting(<subset>, by.risk = 12), conf.int = TRUE, show.cox = FALSE, show.med = FALSE)`;
2. `plot(survfit(Surv(tte, event) ~ treat), conf.int = TRUE, mark.time = TRUE)`, same colours and y-range;
3. overlay without bands: `survfit` curves as thick light lines, the vertices the package draws as thin dark lines on top;
4. the same overlay zoomed on the left tail (x from 0 to 12 months).

Title each page with the subgroup, n and events per arm. If the package call fails, print its message on the page and carry on.

**Difference-curve sheets** — appended to the same PDF: `plotKM.band_subgroups(df, sg_labels = <four subgroups at a time>, draws.band = 0)` beside the same differences computed from `survfit`.

**Numbers** — `dev/audit/stage2/km_plot_compare_<mode>.csv`, one row per subgroup × arm:
`subgroup, arm, n, events, first_event_time, censored_before_first_event, first_drawn_x, first_drawn_y, band_upper_max, max_abs_vertex_diff_curve, max_abs_vertex_diff_band, max_abs_dhat_diff, status`

- "Drawn" means what the package's plotting code puts on the device. In mode `before` that is `at_points[idv]`, `surv[idv]`; in mode `after`, the output of the internal vertex helper from B10 (if B10 was not done, the same as `before`).
- `max_abs_vertex_diff_*` compare the drawn step function with `survfit`'s (`c(0, time)`, `c(1, surv)`, and `lower` / `upper` with `conf.type = "log"`) over `[0, last time]`, evaluated at every jump and at interval midpoints.
- `max_abs_dhat_diff` compares `Dhat_subgroups` at `xpoints` with `S1 − S0` from `survfit`.

### Step 6 — Report

`dev/audit/REPORT_km_fixes_stage2.md`, format in §4. Commit. Then run the reporting line from the kickoff.

---

## 3. Work packages — contracts

Class **B** = changes behaviour: output changes only where it is wrong or failing today. The golden whitelist names every field allowed to move.

### B7 — import `strata` (R6)
- Add `@importFrom survival strata` where `strata()` is used in a formula.
- Check: `"strata" %in% getNamespaceImports("weightedsurv")$survival`.
- Golden whitelist: none.

### B6 — side effects (P2, S7)
- After a call to any exported function, the caller's `.Random.seed` is exactly what it was before — including when it did not exist.
- Results for a given `seedstart` are unchanged: keep the internal `set.seed(seedstart)` and restore the caller's state on exit. Do not change how seeds are shared between arms or strata (S9 is recorded only).
- `get_validated_weights()`: no `cat()`; the information goes into the error message.
- Golden whitelist: none.

### V — housekeeping (S10). Class: moves existing code
- Rename the inner `check_km_curve()` defined inside `df_counting()` so it no longer shadows the export; leave the exported function alone.
- Remove the duplicate `ans$strata` assignment, the unused `dJ` and `n` in `calculate_risk_event_counts()`, and the `requireNamespace()` guards on packages that are in Imports.
- Golden whitelist: none. Output must be identical.

### B4 — input validation and degenerate data (R1, R2, R3a, R3b, R3c, R4a, R4b, R5, R7, S8, N1, N2, N3, N5, N6, N7, N9)

One internal validator shared by `KM_diff()`, `df_counting()` and `plotKM.band_subgroups()`. "Error" below means `stop()` with a message that names the column and the problem.

| Input | Required outcome |
|---|---|
| `NA` in time, event, treatment, weight or strata | Error. Rows are never dropped silently. |
| Non-finite or negative time | Error. Zero is valid. |
| Weight negative | Error. |
| Weight zero | Accepted. The row contributes nothing. KM is as today (already correct). `df_counting()`: the self-check compares on positive-weight rows and raises no false discrepancy; Cox and log-rank are fitted on positive-weight rows. Per-observation output keeps its length. |
| Event | Accepted codings: 0/1, logical, and 1/2 as `survival::Surv` reads it — in all three functions. Anything else: error. |
| Treatment | As today (0/1; factor, logical and character forms that work today keep working). |
| Missing column, including `strata.name` | Error naming the column. |
| One arm empty | Error naming the arm. In `plotKM.band_subgroups()` name the subgroup. |
| A stratum with one arm empty, in stratified `df_counting()` | No error. The stratum contributes nothing to the stratified statistics. |
| An arm with no events; all censored; n = 1 per arm | `KM_diff()`: result, as today. `df_counting()`: returns a result whenever both arms are non-empty; quantities that cannot be computed are `NA` with one warning that says why; the self-checks never fail on such data. A `survdiff` error is caught the same way. |
| `plotKM.band_subgroups()` when an arm has no events and `Maxtau` is `NULL` | Error saying `max_tau` is undefined and `Maxtau` can be supplied. |
| `cox_rhogamma()` when an arm has no events | Value unchanged; add a warning that the estimate is not finite / not identified. |
| `KM_diff(draws.band < 20, show_resamples = TRUE)` | Draws `min(20, draws.band)` resamples. No error. |

Golden whitelist: none (all golden inputs are valid).

### B2 — exact curve values and CI limits (K6, S5)
- `KM_estimates()` itself returns values in [0, 1]; when everything still at risk fails at a time (`dN` equals `ybar` up to relative tolerance 1e-10) survival is exactly 0 from there on, and its variance 0. Same for `survP`, `survG` and the stratum curves.
- Log-transform limits in `km_quantile()` / `kmq_calculations()`: where S = 0 the limits are `NA`, as `survfit` has them, and no "NaNs produced" warning is raised on valid data.
- `plot_km_confint_polygon()`: upper limit capped at 1; where S = 0 the band closes onto the curve; no `NaN` passed to `polygon()`.
- `check_km_curve()` (exported) is unchanged.
- Golden whitelist: none (differences are below 1e-12).

### B1 — Kaplan–Meier at arbitrary time points (K1, K1b, K7, S3, S6)
- `at_points` must be numeric, finite and non-decreasing; otherwise error. Duplicates stay allowed and are returned as given. **Do not change the default expression for `at_points`.**
- For any valid `at_points` — with or without the event times, starting anywhere — `KM_diff()` returns the right-continuous KM value at each point: `surv0`, `surv1` equal `summary(survfit, times = at_points, extend = TRUE)$surv` within 1e-10, and unweighted `sig2_surv*` equal Greenwood at those points.
- How: build, per arm, the grid `G = sort(unique(c(<that arm's observed times>, at_points)))`; run the existing calculations on `G`; return the rows at `match(at_points, G)`. The same for the resampled processes (`surv0_star`, `surv1_star`, `dhat_star`, `Zdhat_star`) and for the `draws > 0` variance. The random draws themselves must not change.
- This is **not** a change of estimator. On a grid that already holds every event time — the default, and everything `df_counting()` uses — output is identical. The weighted variance formula is untouched: it is evaluated on `G` and read off.
- Do it in an internal helper. `calculate_risk_event_counts()` and `KM_estimates()` keep their signatures; add to their documentation that they are Kaplan–Meier only when `at_points` holds every event time of the group.
- `modify_tau = TRUE` builds its grid as today. `plotKM.band_subgroups()` and `cumulative_rmst_bands()` inherit the fix through `KM_diff()`.
- Golden whitelist: every numeric output of the `KM_diff(modify_tau = TRUE)` calls, and `cumulative_rmst_bands()` output only if it was fed such a fit. New requirement on those: `surv0`, `surv1` equal `survfit` at the grid.

### B3 — quantiles (K2, K2b, N4, I4, I5)
- `qprob` keeps its meaning — the time at which S drops to `qprob` — and the documentation says so.
- Point estimate equals `quantile(survfit(...), probs = 1 - qprob)$quantile` in every Stage-1 T4 case, including exact plateaus, plateaus that run to the end of follow-up, and curves that never reach the level (`NA`). This applies to `km_quantile()` / `kmq_calculations()` themselves, weighted and unweighted.
- Unweighted confidence limits equal `survfit(conf.type = "log", conf.int = conf_level)`'s in every T4 case, including non-monotone limit curves (N4).
- The cross-check in `df_counting()` asks `survfit` for `probs = 1 - qprob` at `conf.int = conf_level`. A real discrepancy raises a warning whatever `verbose` is (an error under `stop.onerror`). After this package none is raised on valid data.
- Golden whitelist: `quantile_results`, `quantile_check` and `km_quantile()` output — and each value that moves must equal `survfit`'s.

### B5 — inference details (I1, I2a, I2b, S1, S2, S4)
- **I1.** `c_alpha_band` is the `1 - alpha` quantile of the resampled suprema. At `alpha = 0.05` nothing changes.
- **I2a.** In `wlr_dhat_estimates()` the pooled survival is taken at `tzero`: the value at the last grid point `<= tzero`. Remove the branch that can never be true.
- **I2b.** `tzero` before the first grid point: `dhat = 0`, variance 0. `tzero` at or after the last: values at the last grid point. Never a length-0 result.
- **S1.** `wt.rg.S(scheme = "MB")`: weights are `1 / max(S(t-), S(t*))` with `S(t*)` carried forward when `mb_tstar` is at or beyond the last time; always finite.
- **S2.** Fix `requires_cols`; no partial matching with `$` on the `df_counting` object; Schemper weights without the censoring curve give an error that says what is missing.
- **S4.** `cumulative_rmst_bands()` given a fit made with `draws.band = 0`: error saying a fit with `draws.band > 0` is needed. Leave the integration rule and the hard-coded levels as they are (B8 is held).
- Golden whitelist: in `wlr_dhat_estimates()` output, only quantities that depend on the pooled survival at `tzero` (`sig2_dhat`, its covariance with the weighted log-rank, and anything computed from them). Report the old/new ratio.

### B10 — KM curve rendering (new; only if `B10: GO`)

What chat saw on GBSG with 0.1.0: the curve is drawn from the group's first event time, not from time zero. For `er <= 0` the treated curve starts at (5.55, 0.958) and the control at (6.41, 0.982); nothing is drawn before that, the first drop is missing, and two censoring marks sit at S = 1 with no line under them. The confidence band is joined by straight lines between event times, rises above 1, and has a wall at the first event. The estimates are identical to `survfit`; only the drawing differs.

Contract, for `plot_km_curves_counting()` and therefore `KM_plot_2sample_weighted_counting()` and `plot_weighted_km()`:
- Each curve is drawn as the right-continuous KM step function from `time.zero` (where S = 1) to the group's last observed time. As a step function it equals `survfit`'s — `c(0, time)`, `c(1, surv)` — within 1e-10.
- The confidence band is a step band over the same range: both limits are 1 before the first event; limits as today (log transform), upper capped at 1, closing onto the curve where S = 0. Unweighted, its vertices equal `survfit`'s `lower` / `upper`.
- Censoring marks, axes, risk table, legends, median annotation, colours and every argument: unchanged. `df_counting()` output, including `idv0` / `idv1`, is unchanged.
- The plotting code takes its vertices from one internal helper, so tests can check geometry without a graphics device.
- `plotKM.band_subgroups()` is not changed by B10; Step P records how it compares.
- Golden whitelist: none.

---

## 4. Report format — `dev/audit/REPORT_km_fixes_stage2.md`

Tables first, then a short plain reading.

- **Top line.** Branch, `START_SHA`, final SHA. Which gates passed. Any package NOT DONE, with the reason.
- **A. Packages.** One table: `| ID | class | status (DONE / NOT DONE) | commit | files touched | golden fields that moved | arguments appended | tests added |`.
- **B. Finding register after the fixes.** One row per Stage-1 ID (probe, S-leads, N1–N9): `| ID | before | after | how verified |`.
- **C. Pressure test before / after**, by family T1–T10. Then each remaining failure with its attribution.
- **D. Plot comparison.** One table from `km_plot_compare_before.csv` and `_after.csv`: per subgroup × arm — n, events, first event time, marks stranded before it, first drawn vertex before → after, band maximum before → after, max vertex difference vs `survfit` before → after. Estrogen rows first. Then five lines of reading, and anything in the difference-curve sheets that departs from `survfit`.
- **E. CRAN check and forestsearch.** Status lines verbatim; test count and runtime.
- **F. Policy choices made inside B4** that a user could notice, one line each.
- **G. Gate G2 output.**

Findings only. No follow-up task proposals. If something blocked, say so in the top line.

---

## 5. Gate G2 — post-conditions

Run from the repo root with `START_SHA` set, after the last package commit. Every line must print `OK`. Checks 5 and 6 write `probe_after.txt` and `tests_after.txt`: paste the block's output into the report, commit the report and those two files, then run check 9.

```bash
S2=dev/audit/stage2

# 1. branch, nothing pushed, master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$START_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. held items untouched: trapezoid rule and show_resamples default
git diff "$START_SHA" HEAD -- R/rmst_calculations.R | grep -E '^[-+][^-+].*mid_dhat' && echo "FAIL: RMST integration edited" || echo "OK 2a"
Rscript -e 'pkgload::load_all(".", quiet=TRUE); stopifnot(isTRUE(formals(KM_diff)$show_resamples)); cat("OK 2b\n")'

# 3. vignettes untouched; exports unchanged
git diff --quiet "$START_SHA" HEAD -- vignettes/ && echo "OK 3a" || echo "FAIL: vignettes changed"
diff <(git show "$START_SHA":NAMESPACE | grep '^export' | sort) <(grep '^export' NAMESPACE | sort) > /dev/null && echo "OK 3b" || echo "FAIL: exports changed"

# 4. API and golden snapshots
Rscript $S2/api.R compare    && echo "OK 4a" || echo "FAIL: API"
Rscript $S2/golden.R compare final && echo "OK 4b" || echo "FAIL: golden"

# 5. seed probe: only the held findings still reproduce; controls pass
Rscript dev/audit/km_audit_probe.R --no-write > $S2/probe_after.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S2/probe_after.txt | sort | tr '\n' ' ')" = "I3 K3 K4 K5 " ] && echo "OK 5a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S2/probe_after.txt)" = "3" ] && echo "OK 5b" || echo "FAIL: controls"

# 6. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S2/tests_after.txt 2>&1 && echo "OK 6a" || echo "FAIL: tests"
grep -E '^Status: ' $S2/00check_after.log | grep -Evq 'ERROR|WARNING' && echo "OK 6b" || echo "FAIL: R CMD check"

# 7. deliverables
for f in golden_before.rds api_before.rds km_plots_before.pdf km_plots_after.pdf \
         km_plot_compare_before.csv km_plot_compare_after.csv probe_after.txt 00check_after.log; do
  [ -f "$S2/$f" ] || echo "FAIL: missing $S2/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2.md ] && [ -f NEWS.md ] && [ -d tests/testthat ] && echo "OK 7 (if no FAIL above)"

# 8. record
git log --oneline "$START_SHA"..HEAD
```

```bash
# 9. after the final commit
[ -z "$(git status --porcelain)" ] && echo "OK 9" || echo "FAIL: uncommitted changes"
```

If B10 was not done, say so in the top line; check 5a is unaffected. If a package is NOT DONE, the matching IDs will still reproduce in check 5a: report the check as failed for that reason and name the IDs.
