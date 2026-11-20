# TASK — weightedsurv 0.2.0: Stage 2d, robust variance for the subject-weighted log-rank test and Cox estimate

Written in chat 2026-11-19 from `dev/audit/REPORT_km_oc.md` and `dev/audit/REPORT_km_fixes_stage2c.md`.
Executed by Claude Code in the local `weightedsurv` repo, on branch `km-fixes-0.2.0`.

| | |
|---|---|
| Touches `R/` | **Yes.** L1–L4 **change the method**; B12 **changes behaviour**. |
| Branch | `km-fixes-0.2.0`. No merge. **Never push.** |
| Authorised | **L1** log-rank with subject weights. **L2** `cox_rhogamma()` SE and interval. **L3** `cox_rhogamma()` resampling. **B12** `plotKM.band_subgroups()` default `draws = 0`. **L4** `wlr_dhat_estimates()` — only if the kickoff says `L4: GO`. |
| **Not authorised** | Any change to results when subject weights are absent or all equal. Any change to a point estimate (`lr`, `score`, `bhat`, survival, RMST). The Kaplan–Meier and RMST variance and resampling of Stage 2c. The default of `draws.band`. An option for estimated weights. Any edit under `vignettes/`. The exported API, apart from the one default in B12. |
| Compute | **M6**: 11 cells × 2,000 replicates, log-rank and Cox only; expected 3–8 min on 12 cores. Runs only if the kickoff says `GO`. Stop and record if the projection passes 30 min. |
| Gates | Stop-on-failure, not stop-to-ask. An item that cannot meet its contract is restored and recorded as NOT DONE; the task continues. |

---

## 1. Why

The operating-characteristics run showed, for the log-rank test and Cox estimate with subject weights:

- with propensity weights the released variance gives 0.73–0.75 coverage for the log hazard ratio and rejects a true null 27–29% of the time; the score-residual variance gives 0.94–0.95 and 6–7%;
- with weights unrelated to outcome, 0.895 against 0.950;
- the released z-score changes when the weights are rescaled, and is `NaN` when they are normalised;
- with no weights the two agree (0.955 and 0.954), so nothing is gained by changing unweighted results.

That evidence came from script code that did not handle tied times or strata. This task puts the method into the package, with ties and strata.

Separately, `plotKM.band_subgroups()` passes its default `draws = 20` to `KM_diff()`. With weights that makes the pointwise SE, and the standardisation of the band, a 20-draw variance: band 8% wider (critical value 3.35 against 3.07) with no coverage gain.

---

## 2. The method

One sample (or one stratum). Subject `i`: time `U_i`, event `D_i`, arm `Z_i` (0/1), subject weight `a_i`. `K(s)` is the time weight of the scheme at time `s`, as the package computes it now. Sums over `s` run over **distinct** times.

- `S0(s, b) = sum_j a_j 1(U_j >= s) exp(b Z_j)`
- `zbar(s, b) = sum_j a_j 1(U_j >= s) Z_j exp(b Z_j) / S0(s, b)`
- `dL(s, b) = sum_j a_j D_j 1(U_j = s) / S0(s, b)`

**Per-subject score contribution**

`psi_i(b) = a_i * [ D_i K(U_i) (Z_i − zbar(U_i, b))  −  exp(b Z_i) * sum over s <= U_i of K(s) (Z_i − zbar(s, b)) dL(s, b) ]`

**What follows from it**

| Quantity | Formula |
|---|---|
| Score (package sign) | `− sum_i psi_i(b)`. At `b = 0` this is `lr` / `score` / `u.zero`; it is already what the package computes. |
| Variance of the score | `sum_i psi_i(b)^2` |
| Log-rank z | `score(0) / sqrt(sum_i psi_i(0)^2)` |
| Information | unchanged: `i_bhat` |
| Variance of `bhat` | `sum_i psi_i(bhat)^2 / i_bhat^2` |
| Strata | compute `psi_i` within each stratum; the score and its variance are sums over strata |
| One resampled score | `− sum_i G_ib psi_i(b)`, with the multipliers the code generates today |

**Cumulative version** (`wlr_cumulative(return_cumulative = TRUE)`). The contribution up to time `t` is `psi_i` itself when `U_i <= t`, and `− a_i (Z_i cK(t) − cKz(t))` when `U_i > t`, where `cK(t) = sum_{s <= t} K(s) dL(s)` and `cKz(t) = sum_{s <= t} K(s) zbar(s) dL(s)` at `b = 0`. So

`V(t) = sum_{i: U_i <= t} psi_i^2 + (cK(t) − cKz(t))^2 * sum_{i: U_i > t, Z_i = 1} a_i^2 + cKz(t)^2 * sum_{i: U_i > t, Z_i = 0} a_i^2`

It is a step function that changes at observed times only, and its last value is `sum_i psi_i^2`.

**Tied times.** Events at a tied time are added up (`dL` above); every subject at that time gets the same `S0`, `zbar` and `K`. Where the package holds one time weight per observation, the weight of a tied time is the one at its first observation — the one the score already uses.

**Scope switch.** These formulas apply when the subject weights of the analysis sample are **not all equal** (`length(unique(w)) > 1`), in every stratum. When weights are absent or all equal, today's code path runs unchanged.

**Reference implementation.** Checked in chat on 0.1.0 objects, with and without tied times: it reproduces the package's `u.zero`, `lr` and `i_bhat` to rounding, and equals `coxph(weights = , robust = TRUE, ties = "breslow")` for the SE and the robust score test, stratified included. Put it in `tests/testthat/` as a helper and test against it. The package code may be organised differently.

```r
# time, delta, z (0/1), a = subject weights, K = time weight per observation (aligned with `time`), beta
ws_score_resid_ref <- function(time, delta, z, a, K = rep(1, length(time)), beta = 0) {
  o <- order(time); tt <- time[o]; D <- delta[o]; Z <- z[o]; w <- a[o]; Kk <- K[o]; n <- length(tt)
  first_idx <- which(!duplicated(tt)); grp <- cumsum(!duplicated(tt))        # distinct times and tie-group id
  eb <- exp(beta * Z)
  S0 <- rev(cumsum(rev(w * eb)))[first_idx]                                   # risk-set sums at distinct times
  zb <- rev(cumsum(rev(w * eb * Z)))[first_idx] / S0                          # weighted proportion treated in the risk set
  dNw <- as.vector(rowsum(w * D, grp))                                        # weighted events at each distinct time
  Kd <- Kk[first_idx]                                                         # time weight at each distinct time
  dL <- dNw / S0; cK <- cumsum(Kd * dL); cKz <- cumsum(Kd * zb * dL)
  psi <- numeric(n)
  psi[o] <- w * (D * Kd[grp] * (Z - zb[grp]) - eb * (Z * cK[grp] - cKz[grp]))
  list(psi = psi, U = sum(psi), info = sum(Kd * zb * (1 - zb) * dNw))        # package score = -U ; info = i_bhat
}
```

Guard: a term whose `S0` is 0 is 0.

---

## 3. Ground rules

1. Read the current file on the branch before editing it.
2. One commit per item, message prefixed with its ID. Tests go in the same commit.
3. **Exported API frozen**, apart from the default in B12: no export, argument, other default or returned-element name changes. If `df_counting()` has to hand subject-level data to `wlr_cumulative()`, do it without changing a signature.
4. **Nothing moves for unweighted or equal-weight input.** Enforced by the golden comparison against the branch as it stands at the start of this task.
5. **No point estimate moves anywhere.**
6. No edits under `vignettes/`. Step V only inventories.
7. Earlier files under `dev/audit/` are a record. New material goes in `dev/audit/stage2d/`. Copy a script before changing it.
8. A test that encodes the earlier weighted log-rank or Cox variance is replaced by a test of the new one; list each. Do not weaken any other test.
9. No new package dependency. Error and warning messages you add must not contain `TRUE/FALSE` or `subscript out of bounds`.

---

## 4. Steps

### Step 0 — Intake (Gate G0)
Confirm the branch is `km-fixes-0.2.0` and the tree is clean. Note HEAD as `START_SHA` and `git rev-parse master` as `MASTER_SHA`. `mkdir -p dev/audit/stage2d`; copy this file to `dev/tasks/TASK_km_fixes_stage2d.md`; commit.

**G0 fails (stop) if** the branch is wrong, the tree is not clean, `dev/audit/REPORT_km_oc.md` or `dev/audit/oc/oc_common.R` is missing, or `pkgload::load_all(".")` errors.

### Step 1 — Baseline (Gate G1)
Copy `stage2c/golden.R` and `stage2c/api.R` to `stage2d/`. Before making the baseline, **add to the golden script** any of these calls it does not already make — with weights `wv`, with no weights, and with a constant weight column — so that both what must move and what must not are recorded:

- `df_counting()` with `scheme = "fh"`, `(rho, gamma) = (0, 0)` and `(0, 1)`; and, on the dataset with strata, with `strata.name`;
- `wlr_cumulative(<df_counting result>, ..., return_cumulative = TRUE)`; and `wlr_cumulative()` on a bare data frame of counts (no subject-level columns);
- `wlr_dhat_estimates(<df_counting result>, tzero = <median of the grid>)`;
- `cox_rhogamma()` with `draws = 0` for FH(0,0), FH(0,1) and one `"MB"` and one `"schemper"` call; and with `draws = 200` for FH(0,0);
- `plotKM.band_subgroups(..., draws.band = 200)` with `draws` not supplied, drawing to a null device.

Then `golden.R make` → `stage2d/golden_s2c.rds` (warnings recorded) and `api.R make` → `stage2d/api_s2c.rds`. Create a git worktree of `START_SHA` in a temporary directory for Step V.

**G1 fails (stop) if** `git diff --quiet "$START_SHA" -- R/ NAMESPACE` is false at this point, or a snapshot call errors.

### Step 2 — Items, in order `L1 → L2 → L3 → L4 (if GO) → B12`
Implement to the contracts in §5, add tests, run the tests and `golden.R compare`, commit.

### Step 3 — Documentation (one commit)
- Help pages of `df_counting()`, `wlr_cumulative()`, `wlr_dhat_estimates()` (if L4), `cox_rhogamma()`, `cox_score_rhogamma()`, `cox_rhogamma_resample()` and `plotKM.band_subgroups()`. State: with subject weights that are not all equal, the variance of the weighted score is the sum of squared per-subject score contributions (give the formula); that for FH(0,0) it is the robust variance of `coxph(robust = TRUE)`; that results no longer depend on how the weights are scaled; that with no weights or equal weights results are as before; that weights are treated as known, so intervals are conservative with estimated propensity weights; and that they are too short when weights are extreme.
- `NEWS.md` under 0.2.0: one bullet each for L1–L3 (and L4), saying plainly that z-scores, p-values, standard errors and intervals for subject-weighted analyses change and unweighted ones do not; one bullet for B12.

### Step 4 — M6 (only if the kickoff says `M6: GO`)
Design in §6.

### Step V — Vignette impact inventory (no edits)
For `vignettes/weightedsurv_examples.Rmd`, every file in `vignettes/articles/`, and `README.md`, as in Stage 2c (knit each file with weighted calls under the `START_SHA` worktree and under HEAD):
1. **Calls.** Every chunk that shows a subject-weighted log-rank z or p-value, a `cox_rhogamma()` result with weights, or a weighted `plotKM.band_subgroups()` figure that does not pass `draws`.
2. **Numbers.** Displayed quantities before → after.
3. **Figures.** Each figure that differs, side by side, in `stage2d/vignette_figures_2d.pdf`.
4. **Text.** Every sentence or formula no longer accurate, with file and line.
5. **Weights.** For each weighted example: the code that creates the weight column — the model, its formula, and the form of the weight.

### Step 5 — Verification (Gate G2) and report
Run the block in §8. Re-run a copy of the Stage-2c pressure test into `stage2d/` and tabulate Stage 2c → Stage 2d by family, explaining each change. `R CMD check --as-cran` on a tarball from a clean clone (`stage2d/00check_2d.log`), run as in Stage 2c (`RSTUDIO_PANDOC` set, `DISPLAY` unset). forestsearch check as before. Write `dev/audit/REPORT_km_fixes_stage2d.md` (§7), commit, then run the reporting line from the kickoff.

---

## 5. Contracts

### L1 — log-rank with subject weights. Class: **changes the method**

Where: `wlr_cumulative()`, and through it `df_counting()` (unstratified and per stratum).

- Applies when subject-level data are present and the subject weights are not all equal. `df_counting()` always has them. A direct call to `wlr_cumulative()` has them when its first argument is a `df_counting()` result (elements `time`, `delta`, `z`, `w_hat`). On a bare table of counts the present formula runs unchanged.
- `score` is unchanged. `sig2.score` becomes `sum_i psi_i(0)^2`, and `z.score` follows. With `return_cumulative = TRUE`, `sig2.score` is `V(t)` on the grid.
- `df_counting()` fields that follow: `z.score`, `sig2_lr`, `get_score`, `zlogrank_text`; stratified: `sig2_lr_stratified`, `z.score_stratified`.

| Property | Requirement |
|---|---|
| Reference, FH(0,0) | `z.score^2` equals `coxph(Surv(time, status) ~ treat, weights = w, robust = TRUE, ties = "breslow")$rscore` within 1e-8 relative, with and without tied times. Stratified: `z.score_stratified^2` equals the same from `~ treat + strata(s)`. |
| Other schemes | FH(0,1), MB and schemper: `sig2.score` equals `sum(ws_score_resid_ref(...)$psi^2)` with the scheme's time weights within 1e-8, with and without ties. |
| Cumulative | Equals the direct subject-by-time computation within 1e-10 relative; its last value equals the non-cumulative `sig2.score`. |
| No weights / equal weights | Byte-identical to `START_SHA`. |
| Bare counts | `wlr_cumulative()` on a data frame without subject-level columns: byte-identical to `START_SHA`. |
| Scale | `w`, `10·w`, `w / sum(w)` give the same `z.score` within 1e-9 relative, all finite. |
| Estimates | `lr`, `score`, `lr_stratified` unchanged everywhere. |

### L2 — `cox_rhogamma()` SE and interval. Class: **changes the method**

Where: wherever `cox_rhogamma()` obtains the variance of the score at 0 and at `bhat` (today `score_calculation()` via `cox_score_rhogamma(score_only = FALSE)`), when `w_hat` is not all equal. The root search is unchanged.

- Fields that follow: `fit$sig2_score`, `fit$z.score`, `fit$sig_bhat_asy`, `hr_ci_asy`, `zlogrank_text`, `cox_text_asy`, and what `cox_score_rhogamma(score_only = FALSE)` returns as `sig2_score` and `sig2_beta_asy`.
- Unchanged: `bhat`, `u.beta`, `u.zero`, `i_bhat`, `wt_rg`, `K_wt_rg`.

| Property | Requirement |
|---|---|
| Reference, FH(0,0) | `sig_bhat_asy` equals `sqrt(coxph(..., weights = w, robust = TRUE, ties = "breslow")$var)` within 1e-6 relative; `fit$z.score^2` equals `$rscore` within 1e-8 — with and without tied times. |
| Other schemes | `sig_bhat_asy` equals `sqrt(sum(psi^2)) / info` from `ws_score_resid_ref()` at `bhat` within 1e-8, FH(0,1), MB, schemper, with and without ties. |
| Consistency | `fit$z.score` equals `df_counting()$z.score` for the same scheme within 1e-8, as it does today. |
| No weights / equal weights; scale; estimates | As L1. |

### L3 — `cox_rhogamma()` resampling. Class: **changes the method**

Where: `cox_rhogamma_resample()`, when `w_hat` is not all equal.

- `score_star[b] = − sum_i G_ib psi_i(bhat)` and `score_star_null[b] = − sum_i G_ib psi_i(0)`, where row `i` of `G` is that subject's row in the multiplier matrices generated today (`G0.draws` for controls, `G1.draws` for treated; same seed, dimensions and order).
- Everything computed from them is untouched: `bhat_center_star`, `Z_bstar`, `sig_bhat_star`, `bhat_debiased`, the `wald` fields, `hr_ci_star`, `sig2_score_star`, `z.score_debiased`, `cox_text_star`.

| Property | Requirement |
|---|---|
| Agreement with L2 | With `draws = 4000` on a weighted dataset, `sig_bhat_star / sig_bhat_asy` is within 0.95–1.05 and `var(score_star_null) / sig2_score` within 0.90–1.10. |
| No weights / equal weights | Every output byte-identical to `START_SHA`. |
| Reproducibility | Same seed, same output, twice; the caller's RNG state is treated as it is today. |

### L4 — `wlr_dhat_estimates()`. Class: **changes the method**. Only if `L4: GO`

Where: `wlr_dhat_estimates()` on a `df_counting()` result whose subject weights are not all equal. It reports the log-rank statistic, the survival difference at `tzero`, and their covariance; all three variances have to rest on the same contributions.

- `score` unchanged; `sig2.score` and `z.score` as L1.
- With `d_i` the contribution of subject `i` to the survival difference at the grid point the function reads — `d_i = − S_1 a_i psi^KM_i` in the treated arm, `+ S_0 a_i psi^KM_i` in the control arm, `psi^KM` being the Stage 2c Kaplan–Meier quantity (`ws_km_psi_parts()`):
  `sig2_dhat = sum_i d_i^2`, `cov_wlr_dhat = sum_i (− psi_i(0)) d_i`, `cor_wlr_dhat = cov_wlr_dhat / sqrt(sig2.score * sig2_dhat)`.

| Property | Requirement |
|---|---|
| Identity | `sig2_dhat` equals `KM_diff(draws = 0)$sig2_dhat` at the same grid point within 1e-8 relative. |
| Range | `abs(cor_wlr_dhat) <= 1`. |
| No weights / equal weights | Byte-identical to `START_SHA`. |
| Evidence | M6 compares the estimated correlation with the empirical one (chat check, 2,000 replicates: 0.79 against 0.80 for FH(0,0), 0.51 against 0.52 for FH(0,1)). |

If `L4: HOLD`, leave the function untouched and state in the report that with subject weights it still returns the earlier formulas.

### B12 — `plotKM.band_subgroups()` default `draws`. Class: **changes behaviour**

- The default of `draws` becomes `0` (it is 20). Nothing else in the function changes. Update its help text.
- Effect: with non-constant weights and `draws` not supplied, the pointwise SE and the band's standardisation use the analytic variance. Calls without weights are numerically unaffected.

| Property | Requirement |
|---|---|
| API | The API comparison shows exactly one difference: this default. |
| Unweighted | Golden output of unweighted and equal-weight calls byte-identical. |

Record without changing: the current default of `draws.band` in that function, and whether `forestsearch` passes `draws` or `draws.band` in its calls to it.

### Golden whitelist for this task

A field may move only in a call that passes non-constant subject weights, and only if it is:
- **L1:** `sig2.score`, `z.score`, `sig2_lr`, `sig2_lr_stratified`, `z.score_stratified`, `zlogrank_text`, and the same inside `get_score`;
- **L2:** `sig2_score`, `sig2_beta_asy`, `sig_bhat_asy`, `z.score`, `hr_ci_asy` (`sig_bhat`, `lower`, `upper`), `zlogrank_text`, `cox_text_asy`, `wald`;
- **L3:** `fit_resamples`, `sig_bhat_star`, `bhat_debiased`, `wald_debiased1`, `wald_debiased2`, `hr_ci_star`, `sig2_score_star`, `z.score_debiased`, `cox_text_star`;
- **L4:** `sig2.score`, `z.score`, `sig2_dhat`, `cov_wlr_dhat`, `cor_wlr_dhat` of `wlr_dhat_estimates()`;
- **B12:** in a weighted `plotKM.band_subgroups()` call that does not supply `draws`: variances, pointwise limits, the band critical value and band limits.

Must not move anywhere: `lr`, `score`, `lr_stratified`, `bhat`, `u.beta`, `u.zero`, `i_bhat`, time weights, every Kaplan–Meier and RMST field outside B12, counts, risk tables, and every field of every unweighted or equal-weight call.

`golden.R` also checks each moved FH(0,0) `z.score` and `sig_bhat_asy` against `coxph` as in L1 and L2.

---

## 6. M6 — Monte Carlo with the package code (compute; only if `M6: GO`)

The evaluation's `proposed` columns for the log-rank test and Cox estimate came from script code. This run repeats that part with the package, and adds tied times and strata.

Use the data-generating model, truth and seeds of `dev/audit/oc/oc_common.R` (copy what is needed into `stage2d/`; do not edit files under `oc/`). N = 400, 2,000 replicates per cell, HEAD only.

| k | Cell | Data | Evaluated |
|---|---|---|---|
| 1, 2, 5, 8, 10, 14, 16 | U-400, R-400, D-400, E-400, F-400, D-400-null, F-400-null | exactly the evaluation's cells (seed `1e6 * k + r`) | log HR, log-rank |
| 12 | Dt-400 | the evaluation's cell 12 (tied times) | log HR, log-rank — not run before |
| 21 | Dt-400-null | `gen(400, "null", "conf")`, true weights, `ceiling(time)`; seed `1e6 * 21 + r` | log HR, log-rank |
| 22 | D-400-null-strata | the data of cell 14, with `s = ifelse(X > 0, "hi", "lo")` as `strata.name` | stratified log-rank only |
| 23 | Dt-400-null-strata | the data of cell 21, with the same strata | stratified log-rank only |

Per replicate, for FH(0,0) and FH(0,1): `df_counting()` and `cox_rhogamma(draws = 0)` as in the evaluation (Appendix A5 of its report). First 500 replicates, FH(0,0): `cox_rhogamma(draws = 200)` after `set.seed(3e7 + 1e6 * k + r)`. Cells 5 and 14, if L4: `wlr_dhat_estimates(dc, tzero = 18)`.

Report per cell and scheme: bias (Monte Carlo SE), SE/SD, coverage and mean length for the log HR (asymptotic interval; resampling interval on 500); rejection rate of the two-sided 5% test in null cells; beside them the evaluation's `proposed` values for the same cell from `oc/oc_summary.csv`. For L4: the mean estimated correlation against the empirical correlation of `score` and `dhat`, and the mean of each variance against the empirical variance.

**Required:**
- in cells 2, 5, 8, 10, 14, 16, in every replicate with no tied time after the time fix, the package's `sig_bhat_asy` and `abs(z.score)` equal the evaluation script's `cox_robust()` values (`se`, `abs(z)` at 0) within 1e-8 relative, both schemes;
- in cells 12 and 21, in every replicate, the package's FH(0,0) `sig_bhat_asy` equals `coxph(robust = TRUE, ties = "breslow")` within 1e-6 relative.

Everything else is evidence. Monte Carlo error: coverage ±0.005 at 2,000 and ±0.010 at 500.

---

## 7. Report format — `dev/audit/REPORT_km_fixes_stage2d.md`

Tables first, then a short plain reading.

- **Top line.** Branch, SHAs, gates, items DONE / NOT DONE, whether M6 ran and its wall-clock.
- **A. Items.** `| ID | status | commit | files touched | golden fields that moved | tests added or replaced |`
- **B. In numbers.** One table, rows = golden datasets with weights: log-rank z (FH(0,0) and FH(0,1)) and the `cox_rhogamma()` SE, before → after, with the `coxph` robust value beside the FH(0,0) ones; the stratified z on the strata dataset; the resampling SE before → after.
- **C. M6.** One table, rows = cell × scheme, columns as in §6. At most eight lines of reading.
- **D. Vignette impact inventory**, in the five parts of Step V.
- **E. Recorded, not changed.** The default of `draws.band` and how `forestsearch` calls `plotKM.band_subgroups()`. Every output that, for weighted data, is still computed without the weights or compared with an unweighted reference (for example `df_counting()$logrank_results` and `check_results()`).
- **F. Seed probe, pressure test, API, CRAN check, forestsearch.**
- **G. Gate G2 output.**

Findings only. No recommendation and no follow-up task proposals.

---

## 8. Gate G2 — post-conditions

Run from the repo root with `START_SHA` and `MASTER_SHA` set, after the last item commit. Every line must print `OK`. Checks 6 and 7 write two files: paste the output into the report, commit, then run check 9.

```bash
S=dev/audit/stage2d

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. API: the B12 default and nothing else; vignettes untouched
Rscript $S/api.R compare && echo "OK 2a" || echo "FAIL: API"
git diff --quiet "$START_SHA" HEAD -- vignettes/ && echo "OK 2b" || echo "FAIL: vignettes changed"

# 3. golden: identical outside the Stage-2d whitelist; coxph checks pass
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. L1, L2: against coxph (ties and strata included), scale-free, unchanged at equal weights
Rscript -e 'pkgload::load_all(".", quiet = TRUE); suppressMessages(library(survival))
  set.seed(4); n <- 300; x <- rnorm(n); e <- plogis(-0.5 + x); a <- rbinom(n, 1, e)
  T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60)
  mk <- function(ties) data.frame(time = if (ties) ceiling(pmin(T, cc)) else pmin(T, cc), status = as.numeric(T <= cc), treat = a,
                                  w = ifelse(a == 1, 1 / e, 1 / (1 - e)), s = ifelse(x > 0, "hi", "lo"), one = 1)
  q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
  for (ties in c(FALSE, TRUE)) { d <- mk(ties); d$w10 <- 10 * d$w; d$wn <- d$w / sum(d$w)
    dc  <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
    cx  <- q(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0))
    dcs <- q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s"))
    cp  <- coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
    cps <- coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, robust = TRUE, ties = "breslow")
    z10 <- q(df_counting(d, "time", "status", "treat", weight.name = "w10"))$z.score
    zn  <- q(df_counting(d, "time", "status", "treat", weight.name = "wn"))$z.score
    u   <- q(df_counting(d, "time", "status", "treat")); c1 <- q(df_counting(d, "time", "status", "treat", weight.name = "one"))
    stopifnot(abs(dc$z.score^2 / cp$rscore - 1) < 1e-8, abs(cx$fit$z.score^2 / cp$rscore - 1) < 1e-8,
              abs(cx$fit$sig_bhat_asy / sqrt(cp$var) - 1) < 1e-6, abs(dcs$z.score_stratified^2 / cps$rscore - 1) < 1e-8,
              is.finite(zn), abs(z10 / dc$z.score - 1) < 1e-9, abs(zn / dc$z.score - 1) < 1e-9,
              identical(u$z.score, c1$z.score), identical(u$sig2_lr, c1$sig2_lr)) }
  cat("OK 4\n")'

# 5. L3: the resampling SE agrees with the analytic SE
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(5); n <- 300; x <- rnorm(n); e <- plogis(-0.5 + x); a <- rbinom(n, 1, e)
  T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60)
  d <- data.frame(time = pmin(T, cc), status = as.numeric(T <= cc), treat = a, w = ifelse(a == 1, 1 / e, 1 / (1 - e)))
  q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
  dc <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
  cx <- q(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 4000))
  r <- cx$fit$sig_bhat_star / cx$fit$sig_bhat_asy; stopifnot(is.finite(r), r > 0.95, r < 1.05); cat("OK 5\n")'

# 6. seed probe: statuses as after Stage 2c
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2d.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2d.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 6a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2d.txt)" = "3" ] && echo "OK 6b" || echo "FAIL: controls"

# 7. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2d.txt 2>&1 && echo "OK 7a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2d.log | grep -Evq 'ERROR|WARNING' && echo "OK 7b" || echo "FAIL: R CMD check"

# 8. deliverables and record
for f in golden_s2c.rds api_s2c.rds probe_2d.txt 00check_2d.log; do [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2d.md ] && echo "OK 8 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD
```

```bash
# 9. after the final commit
[ -z "$(git status --porcelain)" ] && echo "OK 9" || echo "FAIL: uncommitted changes"
```

Notes:
- Check 2a: `stage2d/api.R compare` exits 0 only when the sole difference from `api_s2c.rds` is the default of `draws` in `plotKM.band_subgroups()` (0 instead of 20).
- Check 3: `golden.R compare final` uses the union of the whitelists of the items that are DONE.
- Check 4 on the branch as it stands today fails at its first comparison (z 2.75 against a robust 1.74 on that dataset; `NaN` with normalised weights). The reference implementation passes every line of it, with and without tied times.
