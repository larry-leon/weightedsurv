# TASK — weightedsurv 0.2.0: Stage 2c, robust variance and resampling for weighted Kaplan–Meier analyses

Written in chat 2026-11-19 from `dev/audit/REPORT_km_m2.md`.
Executed by Claude Code in the local `weightedsurv` repo, on branch `km-fixes-0.2.0`.

| | |
|---|---|
| Touches `R/` | **Yes.** Two items, both **change the method**, both approved by Larry after M2. |
| Branch | `km-fixes-0.2.0`. No merge. **Never push.** |
| Authorised | **R1** — analytic weighted KM variance becomes the robust (infinitesimal-jackknife) variance. **R2** — resampling of weighted KM curves uses per-subject influence contributions. |
| **Not authorised** | Any change to results when weights are absent or all equal. The log-rank and Cox code (`z.score`, `sig2_lr`, `cox_rhogamma()`, their resampling). The RMST integration rule. Any edit under `vignettes/`. The exported API. |
| Compute | **M2b**: the M2 run again, under 1 min. **M3**: band coverage, 4 cells × 2 code versions × 1,000 replicates, N = 800, 500 draws; expected 10–20 min on 12 cores. **M4**: null calibration of the weighted log-rank and Cox SE, 4 cells × 2,000 replicates; expected under 10 min. Each runs only if the kickoff says `GO`. Stop and record if a projection passes 60 min. |
| Gates | Stop-on-failure, not stop-to-ask. An item that cannot meet its contract is restored and recorded as NOT DONE; the task continues. |

---

## 1. Why

M2 showed that, when weights are related to the hazard (as propensity weights are):
- the weighted Greenwood SE adopted in Stage 2b is off in both directions (SE/SD 0.86 to 1.16 with ordinary propensity weights, 0.63 to 0.73 with extreme ones);
- the package's resampling SE has the same error, because it perturbs the events only;
- the robust variance stays calibrated with true weights and is conservative with estimated ones. An explicit formula reproduces `survfit(robust = TRUE)` to 1e-15 and equals Greenwood at unit weights.

Larry's decision: address this with the most robust and accurate method, for both the analytic SE and the resampling, without disturbing unweighted results.

---

## 2. The method

Notation for one group with times `U_i`, event indicators `D_i`, weights `w_i`; `S` is the weighted KM.

- `Y_w(s) = sum_i w_i 1(U_i >= s)`, `dN_w(s) = sum_i w_i D_i 1(U_i = s)`
- Per-subject influence contribution at time `t`:

  `psi_i(t) = D_i 1(U_i <= t) / (Y_w(U_i) - dN_w(U_i))  -  sum over event times s <= min(t, U_i) of dN_w(s) / ( Y_w(s) * (Y_w(s) - dN_w(s)) )`

  Any term with `Y_w - dN_w <= 1e-10 * Y_w` is 0.

**R1 — analytic variance:** `Var[S(t)] = S(t)^2 * sum_i w_i^2 psi_i(t)^2`

**R2 — one resampled curve,** for multipliers `G_i` (the same standard-normal draws the code generates today): `S*(t) - S(t) = - S(t) * sum_i G_i w_i psi_i(t)`

Its variance over draws is the R1 variance, so the analytic and resampled standard errors agree.

**Scope switch, used by both items:** the formulas above apply when the group's weights are **not all equal**. When weights are absent or all equal, the code path of today is used unchanged — Greenwood for the variance, the present events-only perturbation for the resampling. (At equal weights R1 equals Greenwood to rounding anyway.)

---

## 3. Ground rules

1. Read the current file on the branch before editing it.
2. One commit per item, message prefixed with its ID. Tests go in the same commit.
3. **Exported API frozen**: no export, argument, default or returned-element name changes.
4. **Nothing moves for unweighted or equal-weight input.** Enforced by the golden comparison against the branch as it stands at the start of this task.
5. No edits under `vignettes/`. Step V only inventories what would need revising.
6. Earlier files under `dev/audit/` are a record. New material goes in `dev/audit/stage2c/`. Copy a script before changing it.
7. Tests that encode the weighted Greenwood formula from Stage 2b are replaced by tests of R1; list each one. Do not weaken any other test.
8. Error and warning messages you add must not contain `TRUE/FALSE`, `sig2w_multiplier` or `subscript out of bounds`.

---

## 4. Steps

### Step 0 — Intake (Gate G0)
Confirm the branch is `km-fixes-0.2.0` and the tree is clean. Note HEAD as `START_SHA` and `git rev-parse master` as `MASTER_SHA`. `mkdir -p dev/audit/stage2c`; copy this file to `dev/tasks/TASK_km_fixes_stage2c.md`; commit.

**G0 fails (stop) if** the branch is wrong, the tree is not clean, or `pkgload::load_all(".")` errors.

### Step 1 — Baseline (Gate G1)
With copies of the Stage-2b scripts: `golden.R make` → `stage2c/golden_s2b.rds` (same datasets and calls as Stage 2b, warnings recorded) and `api.R make` → `stage2c/api_s2b.rds`.
Also create a git worktree of `START_SHA` in a temporary directory; M3 and Step V use it as the "before" code.

**G1 fails (stop) if** `git diff --quiet "$START_SHA" -- R/ NAMESPACE` is false at this point, or a snapshot call errors.

### Step 2 — Items, in order `R1 → R2`
Implement to the contracts in §5, add tests, run the tests and `golden.R compare`, commit.

### Step 3 — Documentation (one commit)
- Help pages of `calculate_risk_event_counts()`, `KM_diff()`, `df_counting()`, `resampling_survival()` and `cumulative_rmst_bands()`: replace the Stage-2b text on weighted Greenwood. State: the robust variance and its formula; that it equals `survfit(robust = TRUE)`; that resampled curves for weighted data are built from the same per-subject contributions; that results with no weights or equal weights are unchanged; that weights are treated as known, so the SE is conservative with estimated propensity weights; and that it is too small when weights are extreme (low effective sample size).
- `NEWS.md` under 0.2.0: rewrite the weighted-variance bullet for R1 and add one for R2. Say plainly that weighted standard errors, weighted bands and weighted RMST bands change, and unweighted results do not.

### Step 4 — Evaluations (each only if the kickoff says `GO`)
**M2b.** Run a copy of `dev/audit/m2/km_m2_montecarlo.R` that differs only in its output path (`stage2c/`). Required:
- `pkg` equals `ij`: maximum relative difference below 1e-8 over replicates in which `survfit` merged no times;
- on the first 500 replicates of each arm-cell, `mean(pkg_resampling SE) / mean(pkg SE)` is within 0.97–1.03 at each of the three times.

**M3 — band coverage.** Design in §6.

**M4 — weighted log-rank and Cox, measurement only.** Design in §7. No `R/` change follows from it in this task.

### Step V — Vignette impact inventory (no edits)
For `vignettes/weightedsurv_examples.Rmd`, every file in `vignettes/articles/`, and `README.md`:
1. **Calls.** List every chunk that passes weights to a package function: file, line, function, and what it shows (SE comparison, pointwise band, simultaneous band, RMST, log-rank, Cox).
2. **Numbers.** For each such chunk, compute the displayed quantities under the "before" worktree and under HEAD, reusing the vignette's own code (`knitr::purl`), and tabulate before → after.
3. **Figures.** Draw each affected figure under both versions side by side into `stage2c/vignette_weighted_figures.pdf`.
4. **Text.** List every sentence or formula that is no longer accurate, with file and line — including the two article passages edited in Stage 2b and any description of the resampling that does not hold for weighted data.
5. **Pending scope.** List separately the chunks that display a subject-weighted log-rank z-score or a `cox_rhogamma()` result with weights; these are not changed by this task.

### Step 5 — Verification (Gate G2) and report
Run the block in §9. Re-run the Stage-2b copy of the pressure test into `stage2c/` and tabulate Stage 2b → Stage 2c by family, explaining each change. `R CMD check --as-cran` on a tarball from a clean clone (`stage2c/00check_2c.log`). forestsearch check as before. Write `dev/audit/REPORT_km_fixes_stage2c.md` (§8), commit, then run the reporting line from the kickoff.

---

## 5. Contracts

### R1 — analytic weighted KM variance. Class: **changes the method**

Where: the analytic path of `calculate_risk_event_counts()` (`draws = 0`), for groups whose weights are not all equal.

- `KM_estimates()` is unchanged. It forms `S^2 * cumsum(sig2w_multiplier)`, so return the multiplier as the successive differences of `V(t_j) = sum_i w_i^2 psi_i(t_j)^2` along the grid. Individual differences may be negative; only the cumulative sum is used.
- As before, the result is valid when the grid holds every event time of the group (the Stage-2 internal helper guarantees this for `KM_diff()`).

| Property | Requirement |
|---|---|
| Reference | Through `KM_diff()` and `df_counting()`, the weighted SE equals `survfit(weights = w, robust = TRUE)`'s `std.err` within 1e-8 relative wherever that is finite and positive — with and without tied times. Compare on the time-fixed data (`aeqSurv`), or with `timefix = FALSE` when calling the helpers directly. |
| No weights / equal weights | Byte-identical to `START_SHA`. |
| Scale | `w`, `10·w`, `0.01·w`, `w / sum(w)` give the same SE within 1e-9 relative. |
| Estimates | Survival estimates unchanged everywhere. |
| Warnings | No golden call raises a warning it did not raise at baseline. |

### R2 — resampling of weighted KM curves. Class: **changes the method**

Where: every place a resampled KM curve or a resampling variance is produced for a group whose weights are not all equal —
- `resampling_survival()`, which feeds `surv0_star`, `surv1_star`, `dhat_star` and `Zdhat_star` in `KM_diff()`, and through them the simultaneous bands and `cumulative_rmst_bands()`;
- the `draws > 0` variance path of `calculate_risk_event_counts()`.

Requirements:
- Resampled deviation at grid point `t_j`, draw `b`: `- S(t_j) * sum_i G_ib w_i psi_i(t_j)`. Use the multiplier matrix exactly as it is generated today (same seeds, same dimensions, same order), so that only the way the draws enter changes.
- `draws > 0` variance path: the variance over draws of `sum_i G_ib w_i psi_i(t_j)` at each grid point, returned as successive differences so that `KM_estimates()` reproduces it. Do not sum per-interval variances: these increments are not independent.
- Everything downstream is untouched: the standardisation, the sup statistic, the quantile, `alpha`, the two draw sets and the integration rule in `cumulative_rmst_bands()`.
- Evaluate on the full event-time grid and read off at the requested points, as the Stage-2 helper already does.

| Property | Requirement |
|---|---|
| No weights / equal weights | Every output byte-identical to `START_SHA`, including all resampled curves. |
| Agreement with R1 | With 4,000 draws on a weighted dataset, the variance across draws of each row of `surv0_star` / `surv1_star` is within 0.85–1.15 of the analytic `sig2_surv0` / `sig2_surv1` at every grid point where the survival estimate is between 0.1 and 0.9. |
| `draws > 0` path | With `draws = 4000`, the SE is within 0.93–1.07 of the analytic SE over the same range. |
| Reproducibility | Same seed, same output, twice. The caller's RNG state is unchanged (the Stage-2 behaviour). |

### Golden whitelist for this task

A field may move only in a call that passes non-constant weights, and only if it is: a variance or computed from one (`sig2_*`, standard errors, pointwise limits, weighted quantile confidence limits); a resampled curve (`surv0_star`, `surv1_star`, `dhat_star`, `Zdhat_star`); a band critical value or band limit; or, in `cumulative_rmst_bands()` output, anything computed from the resampled curves.

Must not move anywhere: survival estimates, `dhat`, the RMST point estimate, counts, risk tables, quantile point estimates, **all log-rank and Cox output**, and every field of every unweighted or equal-weight call.

For each analytic `sig2_surv0` / `sig2_surv1` that moves, `golden.R` checks the new value against `survfit` robust.

---

## 6. M3 — simultaneous-band coverage (compute; only if `M3: GO`)

Two arms, N = 800, censoring `U(0, 40)`. `X ~ N(0, 1)`; `T ~ Exp(0.08 * exp(0.7 X - 0.4 A))`.

| Cell | Treatment assignment | Weights |
|---|---|---|
| U | `A ~ Bernoulli(0.5)` | none |
| D | `A ~ Bernoulli(expit(-0.5 + X))` | true propensity weights |
| E | as D | estimated by `glm(A ~ X)` in each replicate |
| F | `A ~ Bernoulli(expit(-0.5 + 2X))` | true propensity weights |

Each replicate (seed `200000 * k + r`), under each code version ("before" = the `START_SHA` worktree, "after" = HEAD):
- `fit <- KM_diff(df, ..., weight.name, draws.band = 500, show_resamples = FALSE)` on the default grid;
- `rm <- cumulative_rmst_bands(df, fit, ..., weight.name, draws_sb = 500, plot = FALSE)`.

Truth: `delta(t) = S_1(t) - S_0(t)` with `S_a(t) = E_X[exp(-0.08 exp(0.7X - 0.4a) t)]`, and `R(t)` its integral from 0, both from a fine grid (step 0.01) by numerical integration.

Record per cell × version:
- pointwise coverage of `lower` / `upper` for `delta` at t = 5, 10, 20;
- simultaneous coverage of `sb_lower` / `sb_upper` for `delta` — over the whole grid, and over grid points in [2, 25];
- for the cumulative RMST difference: pointwise coverage at t = 20, and simultaneous-band coverage over the whole grid and over [2, 25];
- mean band critical value; mean half-width of the simultaneous band at t = 10; mean analytic SE of `dhat` at t = 10 over the empirical SD of `dhat`.

Required: in cell U, "before" and "after" output is identical in every replicate. Everything else is evidence, not pass/fail. Monte Carlo error on coverage at 1,000 replicates is about ±0.007.

---

## 7. M4 — weighted log-rank and Cox under the null (measurement only; only if `M4: GO`)

Chat found on one dataset with propensity weights (0.1.0 code, unchanged since): `df_counting()$z.score` 4.48, against a robust `coxph` z of 2.66; the z-score becomes 14.2 when the weights are multiplied by 10 and `NaN` when they are normalised (`sig2_lr` negative); the `cox_rhogamma()` SE of the log hazard ratio is 0.072 / 0.023 / 3.54 for the same three scalings. This step measures how that plays out. It changes nothing.

Cells: U, D, E, F as in §6 but with **no treatment effect** (`T ~ Exp(0.08 * exp(0.7 X))`), 2,000 replicates each, seed `300000 * k + r`.

Per replicate: `dc <- df_counting(df, ..., weight.name)`, then `cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0)`; and the reference `coxph(Surv(time, status) ~ treat, weights = w, robust = TRUE)`. In cells D, E, F repeat the two package calls with the weights multiplied by 10 and with the weights divided by their sum.

Report per cell (and per weight scaling):
- mean and SD of `z.score`; rejection rate of `|z| > 1.96`; number of non-finite values;
- SD of the `cox_rhogamma()` log-hazard-ratio estimate against the mean of its reported SE;
- the same three quantities for `coxph`: robust z and robust SE; model-based z and SE.

---

## 8. Report format — `dev/audit/REPORT_km_fixes_stage2c.md`

Tables first, then a short plain reading.

- **Top line.** Branch, SHAs, gates, items DONE / NOT DONE, which evaluations ran and their wall-clock.
- **A. Items.** `| ID | status | commit | files touched | golden fields that moved | tests added or replaced |`
- **B. R1 and R2 in numbers.** One table, rows = golden datasets with weights: weighted SE at the 25th / 50th / 75th percentile of the grid, Stage 2b → Stage 2c, with the ratio to `survfit` robust; band critical value before → after; RMST pointwise limits and band limits at the last time, before → after.
- **C. M2b.** The M2 Table 1 again with the new `pkg` and `pkg_resampling` columns, beside `robust`.
- **D. M3.** One table, rows = cell × version, columns = the recorded coverages and widths. At most six lines of reading.
- **E. M4.** One table, rows = cell × weight scaling. At most six lines of reading.
- **F. Vignette impact inventory**, in the five parts of Step V.
- **G. Seed probe, pressure test, API, CRAN check, forestsearch.**
- **H. Gate G2 output.**

Findings only. No recommendation and no follow-up task proposals.

---

## 9. Gate G2 — post-conditions

Run from the repo root with `START_SHA` and `MASTER_SHA` set, after the last item commit. Every line must print `OK`. Checks 6 and 7 write two files: paste the output into the report, commit, then run check 10.

```bash
S=dev/audit/stage2c

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. API identical; vignettes untouched
Rscript $S/api.R compare && echo "OK 2a" || echo "FAIL: API changed"
git diff --quiet "$START_SHA" HEAD -- vignettes/ && echo "OK 2b" || echo "FAIL: vignettes changed"

# 3. golden: identical outside the Stage-2c whitelist; contract checks pass
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. R1: weighted SE equals survfit robust (ties included), is scale-free, and is Greenwood at equal weights
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(1); n <- 150; U <- ceiling(rexp(n, .1) * 2) / 2; D <- rbinom(n, 1, .7); w <- runif(n, .3, 3); at <- sort(unique(U))
  se <- function(w) { re <- calculate_risk_event_counts(U, D, w, at); sqrt(pmax(KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)$sig2_KM, 0)) }
  y  <- survival::Surv(U, D)
  rb <- summary(survival::survfit(y ~ 1, weights = w, robust = TRUE), times = at, extend = TRUE)$std.err
  gw <- summary(survival::survfit(y ~ 1), times = at, extend = TRUE)$std.err
  k  <- is.finite(rb) & rb > 0 & is.finite(gw) & gw > 0
  stopifnot(max(abs(se(w)[k] / rb[k] - 1)) < 1e-8, max(abs(se(w / sum(w))[k] / rb[k] - 1)) < 1e-8, max(abs(se(10 * w)[k] / rb[k] - 1)) < 1e-8,
            max(abs(se(rep(1, n))[k] / gw[k] - 1)) < 1e-8, max(abs(se(rep(2, n))[k] / gw[k] - 1)) < 1e-8); cat("OK 4\n")'

# 5. R2: variance of the resampled weighted curves agrees with the analytic variance
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(2); N <- 300; x <- rnorm(N); e <- plogis(-0.5 + x); a <- rbinom(N, 1, e)
  T <- rexp(N, 0.08 * exp(0.7 * x - 0.4 * a)); cc <- runif(N, 0, 40)
  d <- data.frame(time = pmin(T, cc), status = as.numeric(T <= cc), treat = a, w = ifelse(a == 1, 1 / e, 1 / (1 - e)))
  f <- KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = 4000, show_resamples = FALSE)
  for (g in 0:1) { S <- f[[paste0("surv", g)]]; v <- f[[paste0("sig2_surv", g)]]; vs <- apply(f[[paste0("surv", g, "_star")]], 1, var)
    k <- S > 0.1 & S < 0.9 & v > 0; r <- vs[k] / v[k]; stopifnot(sum(k) > 20, min(r) > 0.85, max(r) < 1.15) }
  cat("OK 5\n")'

# 6. seed probe: only the default RMST rule still differs from survfit
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2c.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2c.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 6a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2c.txt)" = "3" ] && echo "OK 6b" || echo "FAIL: controls"

# 7. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2c.txt 2>&1 && echo "OK 7a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2c.log | grep -Evq 'ERROR|WARNING' && echo "OK 7b" || echo "FAIL: R CMD check"

# 8. log-rank and Cox code untouched
git diff --quiet "$START_SHA" HEAD -- R/weightedCox.R && echo "OK 8" || echo "FAIL: weightedCox.R changed"

# 9. deliverables and record
for f in golden_s2b.rds api_s2b.rds probe_2c.txt 00check_2c.log; do [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2c.md ] && echo "OK 9 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD
```

```bash
# 10. after the final commit
[ -z "$(git status --porcelain)" ] && echo "OK 10" || echo "FAIL: uncommitted changes"
```

Notes:
- Check 2a: the API must be identical to `api_s2b.rds`; the script exits 0 only then.
- Check 8: if a help-page edit requires touching roxygen text in `R/weightedCox.R`, do not make that edit; record it instead.
- If an item is NOT DONE, the checks tied to it fail; report them as failed for that reason and name the item.
