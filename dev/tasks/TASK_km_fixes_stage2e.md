# TASK — weightedsurv 0.2.0: Stage 2e, an option for estimated weights, and four riders

Written in chat 2026-11-19 from `dev/audit/REPORT_km_fixes_stage2d.md`, `dev/audit/REPORT_km_oc.md`, the approved interface proposal (choices A1, B1, C as written) and the branch source at `8f820bb`.
Executed by Claude Code in the local `weightedsurv` repo, on branch `km-fixes-0.2.0`.

| | |
|---|---|
| Touches `R/` | **Yes.** **W1** adds an export. **W2–W4** **change the method**, only when the new option is asked for. **Rd, Rb, Rc** **change behaviour**, each in a stated class of calls. **Ra** adds an output. |
| Branch | `km-fixes-0.2.0`. No merge. **Never push.** |
| Authorised | Riders **Rd** (weights test on the positive weights), **Rb** (two `df_counting()` outputs for weighted data), **Rc** (`draws.band` default), **Ra** (effective sample size). Option **W1** `ps_weights()`, **W2** Kaplan–Meier, bands, RMST, **W3** log-rank, **W4** `cox_rhogamma()`. Each group only if the kickoff says `GO`. |
| **Not authorised** | Any change to a point estimate. Any change to a result of a call that does not ask for `weight_variance = "estimated"`, apart from the three rider changes as whitelisted in §5. Any edit under `vignettes/` or to `README.md`. A warning or message triggered by a low effective sample size. Truncating or rescaling weights. Propensity models other than a logistic `glm`. Subgroup analyses that reuse a full-sample model under `"estimated"`. Changing how a constant weight other than 1 is read. A new package dependency. |
| Compute | **M7**: 15 cells with the package at HEAD — 2,000 replicates in 10 cells, 500 in 5 — and a bootstrap benchmark in 2 cells (500 replicates × 200 resamples). Expected 35–70 min on 12 cores. Runs only if the kickoff says `M7: GO`. Pilot first; do not start if the projection passes 120 min. |
| Gates | Stop-on-failure, not stop-to-ask. A rider that cannot meet its contract is restored and recorded as NOT DONE; the task continues. W1–W4 stand or fall together (§4, Step 2). |

---

## 1. Why

**The option.** Propensity weights are estimated, and the package treats them as known numbers. In the operating-characteristics run (cells with a fitted propensity model, N = 200–800; nominal 0.95):

| Coverage | Weights treated as known | Estimation carried through (script code) |
|---|---|---|
| Survival difference, month 18 | 0.970–0.975 | 0.938–0.948 |
| Simultaneous band, months 3–42 | 0.959–0.986 | 0.929–0.957 |
| RMST difference at 36 | 0.975–0.982 | 0.940–0.950 |
| Log hazard ratio | 0.969–0.981 | 0.939–0.945 |

Intervals were 11–18% shorter and matched a bootstrap that refits the model (0.934 against 0.936). That came from script code, for one weight form, without tied times or strata. This task puts it into the package behind an explicit option whose default is today's behaviour.

**What is not measured yet.** Every weighted example in the vignettes has a treated-arm effective sample size near 0.3 of the arm (Stage 2d report, D5). The evaluation has no cell with weights both estimated and that concentrated. A chat pilot in a setting of that shape (2,000 replicates, script code): a true null is rejected 8.8% of the time with weights treated as known and 10.2% with estimation carried through, at a nominal 5%. M7 measures this with the package, before any vignette uses the option.

**The riders** (Stage 2d report, section E and the zero-weight finding).

- **Rd.** A 0/1 weight column counts as "weights not all equal": z = −0.7569, against −0.7554 for the same rows analysed without weights.
- **Rb.** For weighted data `df_counting()` returns two numbers that do not belong to the weighted analysis: `cox_results$score` (model-based, the weights read as counts: 12.184 against the robust 5.849 on `sim_weighted`) and `logrank_results` (an unweighted `survdiff`). `check_results()` therefore warns on every weighted analysis.
- **Rc.** `plotKM.band_subgroups()` has `draws.band = 20`: a band critical value from 20 draws.
- **Ra.** Nothing reports how concentrated the weights are — the diagnostic that separated the cells where intervals held from the cells where they did not.

---

## 2. The method (W1–W4)

`n` subjects, the rows analysed. Subject `i`: time `U_i`, event `D_i`, arm `A_i` (0/1), covariate row `x_i` of the propensity model (`p` columns), fitted probability `e_i`, weight `w_i`. In the score formulas of Stage 2d the same weight is written `a_i`.

**Propensity model.** A logistic `glm` for `A` fitted on exactly the rows analysed.

- Influence of the coefficients, `n × p`: `IF_i = { sum_j x_j x_j' e_j (1 − e_j) }^{-1} x_i (A_i − e_i)`.
- Weights, `type = "ipw"`: `w_i = 1 / e_i` (treated), `1 / (1 − e_i)` (control).
- Derivative of each weight with respect to the coefficients, `n × p`: `dw_i = −((1 − e_i) / e_i) x_i` (treated), `(e_i / (1 − e_i)) x_i` (control).
- `type = "stabilized"`: with `p1 = mean(A)`, weights and `dw` are those of `"ipw"` multiplied by `p1` (treated) and `1 − p1` (control); `IF` is the same. The arm proportion is treated as a fixed number.

**Kaplan–Meier of a group `g`** (an arm; an arm within a stratum; all subjects). `psi_i(t)` is the Stage 2c influence value of a member of the group (`ws_km_psi_parts()`), computed from the group's own weighted counts.

- `Dg(t) = sum_{k in g} psi_k(t) dw_k` (a `p`-vector for each `t`)
- `Phi_i(t) = 1(i in g) w_i psi_i(t) + IF_i' Dg(t)`, for **every** subject `i = 1..n`
- `var S_g(t) = S_g(t)^2 * sum_{i=1}^n Phi_i(t)^2`

**Difference of the two arms.** `var(S_1(t) − S_0(t)) = sum_i ( S_1(t) Phi^1_i(t) − S_0(t) Phi^0_i(t) )^2`. It is **not** the sum of the two arm variances: the model's coefficients are shared.

**Resampled curves** (`draws.band > 0`). `G` is the `n × B` matrix holding, for each subject, that subject's row of the multipliers generated today (`G0.draws` for controls, `G1.draws` for treated; same seed, dimensions and order). Then

`surv_a_star  <-  surv_a_star − (S_a * Dg_a) %*% (t(IF) %*% G)`

for each arm `a`, where `surv_a_star` is today's matrix, `S_a * Dg_a` is the `m × p` matrix with row `j` equal to `S_a(t_j) Dg_a(t_j)'`, and `t(IF) %*% G` is `p × B`. Everything computed from the resampled curves (difference curves, standardised curves, the band, the RMST variance and band) then follows unchanged.

**Score of the log-rank test and Cox estimate.** `psi_i(b)` is the Stage 2d contribution (`ws_score_psi()`), computed within each stratum.

- `B(b) = sum_{j=1}^n (psi_j(b) / a_j) dw_j` — all subjects, all strata
- `psi^ew_i(b) = psi_i(b) + IF_i' B(b)`
- variance of the score: `sum_{i=1}^n psi^ew_i(b)^2`; log-rank z, variance of `bhat` (`/ i_bhat^2`) and one resampled score (`− sum_i G_ib psi^ew_i(b)`) as in Stage 2d with `psi^ew` in place of `psi`.
- With strata the variance is this single sum over all subjects. It is **not** a sum of per-stratum variances.

**Cumulative log-rank** (`wlr_cumulative(return_cumulative = TRUE)`). With `c_i(t)` the Stage 2d contribution of subject `i` up to `t` (`psi_i` when `U_i <= t`, `− a_i (Z_i cK(t) − cKz(t))` otherwise):

`V^ew(t) = sum_i ( c_i(t) + IF_i' B(t) )^2`,  `B(t) = sum_j (c_j(t) / a_j) dw_j`

`= V(t) + 2 B(t)' C(t) + B(t)' { sum_i IF_i IF_i' } B(t)`, with `C(t) = sum_i c_i(t) IF_i`. `B(t)` and `C(t)` are cumulative sums over the grid:
`B(t) = sum_{U_j <= t} (psi_j / a_j) dw_j − cK(t) sum_{U_j > t, Z_j = 1} dw_j + cKz(t) sum_{U_j > t} dw_j`,
`C(t) = sum_{U_i <= t} psi_i IF_i − cK(t) sum_{U_i > t, Z_i = 1} a_i IF_i + cKz(t) sum_{U_i > t} a_i IF_i`.
Its last value is the non-cumulative variance; its value at `t` is the non-cumulative variance of the data censored at `t`.

**Covariance of the log-rank statistic with the survival difference** (`wlr_dhat_estimates()`). With `d_i` the Stage 2d contribution to the difference at the grid point read for `tzero`: `d^ew_i = d_i + IF_i' ( S_0 Dg_0 − S_1 Dg_1 )` at that grid point; `sig2_dhat = sum_i (d^ew_i)^2`, `cov_wlr_dhat = sum_i (− psi^ew_i(0)) d^ew_i`.

**Scope and limits of this version** (to be stated in the help pages).

- A logistic `glm` for a binary treatment, fitted on exactly the rows analysed. Any other weight model: `"known"`.
- Time weights of the scheme (FH and the others) depend on the weighted pooled curve; that dependence on the estimated weights is ignored.
- With `"stabilized"` the arm proportion is treated as fixed. The Kaplan–Meier quantities do not depend on that (they are unchanged when an arm's weights are multiplied by a constant).
- With `"estimated"` the influence-based formulas always run; the test of whether the weights are all equal is not applied.

**Reference implementation.** Checked in chat against the branch, as extracted from this document: the Kaplan–Meier part equals the package's known-weight variance when the second term is dropped (3e-15) and, for arms, a stratum group and all subjects, the helper-based form of the implementation note below; the score part equals the same adjustment applied to `coxph()`'s own residuals, tied times (Breslow) and strata included (8 digits); the cumulative form equals the direct sum (1e-15). Put it in `tests/testthat/` as a helper and test against it. The package code may be organised differently.

```r
# Propensity model: weights, influence of the coefficients (n x p), derivative of each weight (n x p)
ws_ps_ref <- function(formula, data, type = c("ipw", "stabilized")) {
  type <- match.arg(type)
  g  <- stats::glm(formula, family = stats::binomial(), data = data)
  xm <- stats::model.matrix(g); e <- as.numeric(stats::fitted(g)); A <- as.numeric(g$y)
  IF <- (xm * (A - e)) %*% solve(crossprod(xm * sqrt(e * (1 - e))))        # coefficient estimate minus its limit ~ colSums(IF)
  w  <- ifelse(A == 1, 1 / e, 1 / (1 - e))
  dw <- ifelse(A == 1, -(1 - e) / e, e / (1 - e)) * xm
  if (type == "stabilized") { p1 <- mean(A); sc <- ifelse(A == 1, p1, 1 - p1); w <- sc * w; dw <- sc * dw }
  list(weights = w, IF = IF, dw = dw, treat = A)
}

# Kaplan-Meier of the group `rows` (row numbers), on a sorted grid that holds the group's event times.
# time, delta, w: all n subjects; IF, dw: n x p in the same row order. Returns S, Phi (n x m) and the variance.
ws_km_ew_ref <- function(time, delta, w, IF, dw, rows, grid) {
  U <- time[rows]; D <- delta[rows]; W <- w[rows]
  ybar <- colSums(outer(U, grid, ">=") * W)
  dN   <- diff(c(0, colSums(outer(U[D == 1], grid, "<=") * W[D == 1])))
  ok   <- ybar > 0 & (ybar - dN) > 1e-10 * ybar
  psi  <- t(apply(sweep(outer(U, grid, "==") * D, 2, ifelse(ok, 1 / (ybar - dN), 0), "*") -
                  sweep(outer(U, grid, ">="), 2, ifelse(ok, dN / (ybar * (ybar - dN)), 0), "*"), 1, cumsum))   # n_g x m
  S    <- cumprod(1 - ifelse(ybar > 0, dN / ybar, 0))
  Dg   <- crossprod(psi, dw[rows, , drop = FALSE])                          # m x p
  Phi  <- IF %*% t(Dg); Phi[rows, ] <- Phi[rows, ] + W * psi                # n x m
  list(S = S, Dg = Dg, Phi = Phi, sig2 = S^2 * colSums(Phi^2))
}
# Difference of the arms, k1 and k0 from ws_km_ew_ref(), n = nrow(IF):
ws_km_ew_diff_ref <- function(k1, k0) { n <- nrow(k1$Phi); colSums((k1$Phi * rep(k1$S, each = n) - k0$Phi * rep(k0$S, each = n))^2) }

# Score contributions with estimated weights: psi from ws_score_resid_ref() (within stratum), a = weights, all n subjects
ws_score_ew_ref <- function(psi, a, IF, dw) as.numeric(psi + IF %*% colSums((psi / a) * dw))
```

**The same through `survival`** (used as the independent check, and for the `survival`-based outputs of `df_counting()`): for a fitted `coxph(..., weights = w, robust = TRUE)` and a vector `r` of its weighted residuals, `adj(r) = r + IF %*% colSums((r / w) * dw)`.

- Log-rank chi-square at 0: `r = residuals(fit0, type = "score", weighted = TRUE)` from a fit with `init = 0, iter.max = 0`; chi-square `= sum(r)^2 / sum(adj(r)^2)`. With `ties = "breslow"` it equals the package's `z.score^2` for FH(0,0), strata included.
- Standard error of the log hazard ratio: `r = residuals(fit, type = "dfbeta", weighted = TRUE)`; `se = sqrt(sum(adj(r)^2))`.

**Implementation note** (not a contract). The existing helpers give everything from `m × p` matrices, with no `n × m` matrix: for a group with `parts <- ws_km_psi_parts(U, D, grid, ybar, nbar)` and weights `W`, `Dg = ws_km_ij_draws(parts, W, dw_g / W)`, `A = ws_km_ij_draws(parts, W, IF_g)`, `V = ws_km_ij_var(parts, W)`; then `sum_i Phi_i^2 = V + 2 rowSums(A * Dg) + rowSums((Dg %*% crossprod(IF)) * Dg)`, and for two arms `sum_i Phi^1_i Phi^0_i = rowSums(A1 * Dg0) + rowSums(A0 * Dg1) + rowSums((Dg1 %*% crossprod(IF)) * Dg0)`. Checked in chat against the reference (3e-15). A stand-in built this way reproduced the evaluation's script values in 12 replicates of cells 7, 8 and 15 (largest relative difference 4e-13), bands and RMST included.

---

## 3. Ground rules

1. Read the current file on the branch before editing it.
2. One commit per item, message prefixed with its ID. Tests go in the same commit.
3. **Exported API**: only the differences listed in §5 ("API").
4. **Default behaviour is frozen.** A call that does not pass `weight_variance = "estimated"` returns what it returns at `START_SHA`, apart from the rider whitelists of §5 and the new elements listed there. Enforced by the golden comparison.
5. **No point estimate moves anywhere**, under either option.
6. No edits under `vignettes/` or to `README.md`. Step V only inventories.
7. Earlier files under `dev/audit/` are a record. New material goes in `dev/audit/stage2e/`. Copy a script before changing it.
8. A test that encodes behaviour a rider changes is replaced by a test of the new behaviour; list each. Do not weaken any other test.
9. No new package dependency. Error and warning messages you add must not contain `TRUE/FALSE` or `subscript out of bounds`.
10. New arguments go at the **end** of each argument list, in the order `weight_variance`, `weight_model`.

---

## 4. Steps

### Step 0 — Intake (Gate G0)
Confirm the branch is `km-fixes-0.2.0` and the tree is clean. Note HEAD as `START_SHA` and `git rev-parse master` as `MASTER_SHA`. `mkdir -p dev/audit/stage2e`; copy this file to `dev/tasks/TASK_km_fixes_stage2e.md`; commit.

**G0 fails (stop) if** the branch is wrong; the tree is not clean; `git diff --quiet 8f820bb HEAD -- R/ NAMESPACE` is false (this document was written against that source); `dev/audit/REPORT_km_fixes_stage2d.md`, `dev/audit/oc/oc_common.R`, `dev/audit/oc/oc_sim.R`, `dev/audit/oc/oc_boot.R` or `dev/audit/stage2d/golden.R` is missing; or `pkgload::load_all(".")` errors.

### Step 1 — Baseline (Gate G1)
Copy `stage2d/golden.R` and `stage2d/api.R` to `stage2e/`. Before making the baseline, **add to the golden script** these calls, so that what must move and what must not are both recorded:

- on one weighted dataset, with a 0/1 weight column and with a column of unequal positive weights and about 30% zeros: `df_counting()` FH(0,0) (and with `strata.name` on a dataset with strata), `KM_diff(draws.band = 200)`, `cox_rhogamma()` with `draws = 0` and `draws = 200`;
- `check_results(<df_counting result>, verbose = FALSE)` for a weighted and an unweighted result;
- `plotKM.band_subgroups()` with neither `draws` nor `draws.band` supplied, weighted and unweighted, to a null device.

Then `golden.R make` → `stage2e/golden_s2d.rds` (warnings recorded) and `api.R make` → `stage2e/api_s2d.rds`. Create a git worktree of `START_SHA` in a temporary directory for Step V.

**G1 fails (stop) if** `git diff --quiet "$START_SHA" -- R/ NAMESPACE` is false at this point, or a snapshot call errors.

### Step 2 — Items, in order `Rd → Rb → Rc → Ra → W1 → W2 → W3 → W4`
Implement to the contracts in §5, add tests, run the tests and `golden.R compare`, commit.

**W1–W4 stand or fall together.** If W2, W3 or W4 cannot meet its contract, `git revert` the W commits already made (no history rewriting), record W1–W4 as NOT DONE, skip M7 and Gate checks 6–9, and continue with the riders' documentation, Step V parts 1–3 and the report.

### Step 3 — Documentation (one commit)
- A help page for `ps_weights()` with an example that runs in under 2 seconds (for instance on `survival::gbsg`, `hormon` as treatment). Add it to the reference index in `_pkgdown.yml`.
- One shared section, "Estimated weights", inherited by the help pages of `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()`, `plotKM.band_subgroups()`, `wlr_cumulative()`, `wlr_dhat_estimates()` and `cox_rhogamma()`. It states: what `"known"` and `"estimated"` mean; the formulas of §2 in words and symbols; that estimates never change, only standard errors, intervals, bands and tests; the scope and limits of §2; that with `"estimated"` the variance of a difference is not the sum of the arm variances; that `draws > 0` is not available with it; that subgroup analyses reusing weights from a full-sample model must use `"known"`; and the reference — León, L. F., Cai, T. and Wei, L. J. (2009), Robust inferences for covariate effects on survival time with censored linear regression models, *Statistics in Biosciences* 1, 50–64, doi:10.1007/s12561-009-9002-5 — whose device of carrying the estimation of a nuisance quantity into the inference by perturbing its estimate is used here in linearised form.
- `df_counting()`: the new elements (`ess`, `weight_variance`, `ew_parts`, `cox_results$se`, `logrank_results$method`); what `ess` is and how to read it — the ratio `ess / n` near 1 means the weights are even, a low ratio means a few subjects carry most of the weight, and standard errors, intervals and tests can then be too optimistic under either option; state no threshold. What `logrank_results` and `cox_results$score` are with and without weights (Rb).
- `plotKM.band_subgroups()`: the default of `draws.band` (Rc). Wherever a help page says that zero-weight rows contribute nothing, say that such rows no longer affect which variance formula runs (Rd).
- Replace the sentence "with estimated propensity-score weights the intervals are conservative" wherever it stands by a pointer to the new option.
- `NEWS.md` under 0.2.0: one bullet for the option (new, default unchanged), one for each rider, each saying plainly which calls change and which do not.

### Step 4 — M7 (only if the kickoff says `M7: GO` and W1–W4 are DONE)
Design in §6.

### Step V — Vignette impact inventory (no edits)
For `vignettes/weightedsurv_examples.Rmd`, every file in `vignettes/articles/`, and `README.md`:

1. **Riders.** Every chunk that (a) calls `plotKM.band_subgroups()` without `draws.band`; (b) shows `cox_results$score`, `logrank_results` or `check_results()` output for weighted data; (c) analyses a weight column that contains a zero. For each: the displayed quantities before → after (knit the file under the `START_SHA` worktree and under HEAD, as in Stage 2d), and each figure that differs, side by side.
2. **Where the option could apply.** For each weighted example: the code that makes the weight column; whether `ps_weights()` with the same formula and `type` reproduces that column (largest absolute difference); and, for each analysis call that uses the column, whether it analyses exactly the model's rows.
3. **Effective sample size** by arm for each weighted example (`df_counting()$ess`).
4. **Known against estimated.** For each weighted example where part 2 holds, one table: log-rank z and p-value (`df_counting()`), hazard ratio and interval (`cox_rhogamma()`, FH(0,0)), standard error of the survival difference at the median follow-up time, band critical value, and RMST difference with its interval at the last time — under `"known"` and under `"estimated"`. And the weighted `plotKM.band_subgroups()` and `cumulative_rmst_bands()` figures under the two options, side by side, in `stage2e/vignette_known_vs_estimated.pdf`. Skip part 4 if W1–W4 are NOT DONE.
5. **Text.** Every sentence or formula no longer accurate after this task, with file and line.

### Step 5 — Verification (Gate G2) and report
Run the block in §8. Re-run a copy of the pressure test into `stage2e/` and tabulate Stage 2d → Stage 2e by family, explaining each change. `R CMD check --as-cran` on a tarball from a clean clone (`stage2e/00check_2e.log`), run as in Stage 2d. forestsearch check as before. Write `dev/audit/REPORT_km_fixes_stage2e.md` (§7), commit, then run the reporting line from the kickoff.

---

## 5. Contracts

### Rd — the weights test looks at the positive weights. Class: **changes behaviour** (calls whose weight column contains a zero)

Where: the seven places that test `length(unique(<weights>)) > 1` — `df_counting()` (`lr_unequal`), `ws_lr_subject_data()`, `calculate_risk_event_counts()` (twice), `resampling_survival()`, `cox_score_rhogamma()`, `cox_rhogamma_resample()`. Use one internal helper: the weights are "not all equal" when the **positive** weights are not all equal.

| Property | Requirement |
|---|---|
| Zero-weight rows are inert | With a 0/1 weight column, and with zeros among unequal positive weights, results equal those of the positive-weight rows analysed with the same weights (no weights for 0/1): Kaplan–Meier variances at common time points, `lr`, `sig2_lr`, `z.score`, the stratified versions, `KM_diff()$sig2_dhat`, and `cox_rhogamma()` `bhat`, `sig_bhat_asy`, `z.score` — within 1e-10 relative, with and without tied times. (Chat check on a patched copy: exactly equal.) |
| Resampled quantities | Not required to equal the subset's draw by draw — the multiplier matrices keep a row per subject. Record one example of each (band critical value, `sig_bhat_star`). |
| No zero weight | Byte-identical to `START_SHA`. |

### Rb — `cox_results$score` and `logrank_results` with subject weights. Class: **changes behaviour** (calls whose positive weights are not all equal)

Where: `df_counting()`. Both outputs come from the `survival` package and stay that way.

- `logrank_results`: the robust score test at 0 of the weighted Cox model — `coxph()` on the positive-weight rows with the weights, the same formula as `cox_fit` (strata included), `init = 0, iter.max = 0, robust = TRUE, ties = "breslow"`; `chisq` is its `rscore`, `pval = 1 − pchisq(chisq, 1)`, `logrank_text` as today from that p-value. A failure is handled as the `survdiff` failure is today (`NA`, one reason).
- New element `logrank_results$method`: `"survdiff"` or `"coxph_robust_score"`.
- `cox_results$score`: `summary(cox_fit)$robscore[1]`.
- New element `cox_results$se`: the standard error of the log hazard ratio behind `hr_ci` and `pval` (`sqrt(cox_fit$var)`; `NA` when `cox_results` is the `NA` placeholder).
- With subject weights `survdiff()` has no counterpart, so `logrank_results` is then the `rho = 0` test whatever the scheme. Document it.
- No weights, or all positive weights equal: unchanged (`survdiff`, `sctest`). The comparisons at the end of `df_counting()` and the code of `check_results()` are unchanged.

| Property | Requirement |
|---|---|
| Reference | `logrank_results$chisq` equals `coxph(..., weights = w, robust = TRUE, ties = "breslow")$rscore` within 1e-8 relative, stratified too, with and without tied times; and equals `z.score^2` (`z.score_stratified^2`) for FH(0,0) within 1e-8. `cox_results$score` equals `summary(coxph(..., weights = w, robust = TRUE))$robscore[1]`. |
| `check_results()` | On a weighted FH(0,0) result its three numbers agree and it does not warn. Report its output before → after on `sim_weighted`. |
| No weights / equal weights | `cox_results$score`, `logrank_results$chisq`, `$pval`, `$logrank_text` byte-identical to `START_SHA`. |

### Rc — `plotKM.band_subgroups()` default `draws.band`. Class: **changes behaviour** (calls that do not supply `draws.band`)

- The default of `draws.band` becomes `1000` (it is 20). Nothing else in the function changes. Update its help text.
- List the examples and tests that rely on the default, with their run time before → after. If an example then takes more than 5 seconds, pass `draws.band` explicitly in that example and say so.

| Property | Requirement |
|---|---|
| API | This default is one of the listed API differences. |
| Explicit `draws.band` | Golden output of every call that supplies `draws.band` is byte-identical. |

### Ra — effective sample size. Class: **adds an output**

Where: `df_counting()`. New element `ess`, a data frame with one row per arm in the order of `quantile_results` (`arms[1]` = treated, then `arms[2]`): `group`, `n` (rows of the arm with positive weight), `ess` = `(sum w)^2 / sum(w^2)` over those rows, `ratio` = `ess / n`. Without weights `ess = n`. Present in every result. No warning or message depends on it.

| Property | Requirement |
|---|---|
| Values | Equal the direct computation within 1e-10; `ess = n` without weights and with constant weights. |
| Nothing else | No existing element changes. |

### W1 — `ps_weights()`. Class: **adds an export**

`ps_weights(formula = NULL, data = NULL, type = c("ipw", "stabilized"), model = NULL)`; new file `R/ps_weights.R`.

- Either `formula` and `data` — fitted by `stats::glm(formula, family = stats::binomial(), data = data)` with default control, so that weights computed by hand from the same call agree — or `model`, an existing fit. Not both.
- Errors, each naming the problem: not a `glm`; family not binomial or link not logit; prior weights other than 1 or an offset; not converged; a fitted probability within 1e-10 of 0 or 1; rows dropped for missing values; a response that does not take both values 0 and 1.
- Returns an object of class `"ps_weights"`: `weights` (length `n`), `type`, `treat` (0/1, length `n`), `ps` (fitted probabilities), `IF` (`n × p`), `dw` (`n × p`), `model`, `call`. `p` counts the estimated coefficients (aliased columns dropped).
- `print.ps_weights()`: type, `n`, number of coefficients, and by arm the size, the range of the weights and the effective sample size. Returns its argument invisibly.

| Property | Requirement |
|---|---|
| Reference | `weights`, `IF`, `dw` equal `ws_ps_ref()` within 1e-10, both types. The `model` route equals the `formula` route. |
| Derivatives | `dw` equals the central-difference derivative of the weights in the coefficients within 1e-5 relative. `IF` rows equal the derivative of the fitted coefficients in the case weight of that subject (refit with weight `1 + 1e-4` on one subject, fits iterated to convergence) within 1e-3 relative, for five subjects. |
| Stabilised | `weights` and `dw` are the `"ipw"` ones times the arm proportion; `IF` identical. |

### W2 — Kaplan–Meier, bands, RMST. Class: **changes the method** (only with `weight_variance = "estimated"`)

New arguments `weight_variance = c("known", "estimated")`, `weight_model = NULL` on `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()`, `plotKM.band_subgroups()`.

**Validation** (one internal helper), with `"estimated"`: `weight.name` is given; `weight_model` inherits from `"ps_weights"`; `nrow(df)` equals the length of its `weights`; the treatment column equals its `treat`; the weight column equals its `weights` (`all.equal`, tolerance 1e-8); `draws` is 0. Each failure is an error naming the argument and the problem. With `"known"`, `weight_model` is ignored. The rows of `IF` and `dw` follow every reordering of the data (each function sorts by time after the time fix).

**`KM_diff()`.** Fields that follow §2: `sig2_surv0`, `sig2_surv1`, `sig2_dhat`, `lower`, `upper`, `surv0_star`, `surv1_star`, `dhat_star`, `Zdhat_star`, `c_alpha_band`, `sb_lower`, `sb_upper`. New element `weight_variance`, in every result.

**`df_counting()`.** Every Kaplan–Meier variance it returns: `sig2_surv0`, `sig2_surv1`, `sig2_survP`, and `sig2_surv0_stratum`, `sig2_surv1_stratum`, `sig2_survP_stratum` (group = the arm, or all subjects, within the stratum). `quantile_results` limits follow. With `check.seKM = TRUE` the comparison plot is still drawn; the discrepancy warning for the standard error is not raised (the two differ by design). New elements: `weight_variance` in every result; with `"estimated"`, `ew_parts = list(IF, dw)` with rows in the order of `time`. Confirm by reading that `plot_weighted_km()` takes its standard errors from the result; report.

**`cumulative_rmst_bands()`.** Passes both arguments to its own `KM_diff()` call. Error if `fit` was made with the other option (a fit without the element counts as `"known"`). New element `weight_variance`.

**`plotKM.band_subgroups()`.** The overall `KM_diff()` call gets both arguments. The subgroup curves are point estimates and are computed as today.

| Property | Requirement |
|---|---|
| Reference | `sig2_surv0`, `sig2_surv1`, `sig2_dhat` of `KM_diff()` equal `ws_km_ew_ref()` / `ws_km_ew_diff_ref()` within 1e-8 relative, with and without tied times; the same for `df_counting()` at common time points, for `sig2_survP`, and for each stratum group. |
| Resampling | With `draws.band = 4000`, the variance over draws of `dhat_star` over `sig2_dhat`, median over grid points with positive variance, is within 0.95–1.05. `dhat_star` equals `− t(S_1 Phi^1 − S_0 Phi^0) %*% G` within 1e-8 relative. |
| Weight form | `"ipw"` and `"stabilized"` objects from the same model give the same `sig2_surv0`, `sig2_surv1`, `sig2_dhat`, `dhat_star`, `c_alpha_band` within 1e-8. |
| Reduces to known | A `ps_weights` object with its `IF` set to zero gives the `"known"` results within 1e-10. |
| Estimates | `surv0`, `surv1`, `dhat`, `rmst_time` and every count identical to the `"known"` call. |
| Default | Every call without `weight_variance = "estimated"`: as ground rule 4. |
| Misuse | Each validation failure is an error; so is a `fit` of the other option in `cumulative_rmst_bands()`. |

### W3 — log-rank. Class: **changes the method** (only with `"estimated"`)

Where: `df_counting()` (unstratified and stratified); `wlr_cumulative()` and `wlr_dhat_estimates()` when their first argument is a `df_counting()` result made with `"estimated"` (they read `weight_variance` and `ew_parts` from it; their signatures do not change; on a bare table of counts nothing changes).

- `df_counting()` fields that follow §2: `sig2_lr`, `z.score`, `get_score`, `zlogrank_text`; stratified `sig2_lr_stratified`, `z.score_stratified`.
- The `survival`-based outputs carry the same adjustment (§2, "The same through `survival`"): `logrank_results$chisq`, `$pval`, `$logrank_text` from the Breslow fit at 0; `cox_results$score` from the fit at 0 with `cox_fit`'s ties; `cox_results$se`, `hr_ci`, `pval`, `cox_text` from `cox_fit`'s `dfbeta` residuals (`hr_ci = exp(coef ± qnorm(0.975) se)`, same shape and names as today; `pval` two-sided normal). `hr` and `cox_fit` are unchanged.
- `wlr_cumulative()`: `sig2.score` and `z.score`, cumulative or not. `wlr_dhat_estimates()`: `sig2.score`, `z.score`, `sig2_dhat`, `cov_wlr_dhat`, `cor_wlr_dhat`.

| Property | Requirement |
|---|---|
| Reference, FH(0,0) | `z.score^2` and `logrank_results$chisq` equal the chi-square from `coxph()`'s score residuals (§2) within 1e-8 relative, with and without tied times; stratified likewise with `strata(s)`. |
| Other schemes | FH(0,1), MB, schemper: `sig2_lr` equals `sum(ws_score_ew_ref(psi, ...)^2)` with `psi` from `ws_score_resid_ref()` and the scheme's time weights, within 1e-8. |
| `cox_results` | `cox_results$se` equals `sqrt(sum(adj(dfbeta)^2))` from `coxph(..., robust = TRUE)` within 1e-8; `hr_ci` and `pval` follow from it. |
| Cumulative | Equals the direct subject-by-time sum within 1e-10 relative; its last value equals `sig2_lr`; its value at a time `t` equals `sig2_lr` of the data censored at `t` (same weights and model) within 1e-8. |
| `wlr_dhat_estimates()` | `sig2_dhat` equals `KM_diff(weight_variance = "estimated")$sig2_dhat` at the same grid point within 1e-8; `abs(cor_wlr_dhat) <= 1`. |
| Estimates | `lr`, `score`, `lr_stratified`, `hr` identical to the `"known"` call. |
| Reduces to known | `IF` set to zero: the `"known"` results within 1e-10. |

### W4 — `cox_rhogamma()`. Class: **changes the method** (only with `"estimated"`)

Where: `cox_rhogamma()` when `dfcount` was made with `"estimated"`, and `cox_rhogamma_resample()` through the `fit` it is handed. No signature changes. A direct call to `cox_score_rhogamma()` is unchanged.

- Fields that follow §2: `fit$sig2_score`, `fit$z.score`, `z.score`, `fit$sig_bhat_asy`, `hr_ci_asy`, `zlogrank_text`, `cox_text_asy`; with `draws > 0`, `score_star`, `score_star_null` and everything computed from them. Unchanged: `bhat`, `u.beta`, `u.zero`, `i_bhat`, time weights.
- New element `weight_variance` in every result.

| Property | Requirement |
|---|---|
| Reference, FH(0,0) | `fit$sig_bhat_asy` equals the `dfbeta` route of §2 with `ties = "breslow"` within 1e-6 relative; `fit$z.score^2` equals the score-residual route within 1e-8 — with and without tied times. |
| Other schemes | `fit$sig_bhat_asy` equals `sqrt(sum(ws_score_ew_ref(psi(bhat), ...)^2)) / info` within 1e-8, FH(0,1), MB, schemper. |
| Consistency | `fit$z.score` equals `df_counting()$z.score` for the same scheme within 1e-8. |
| Resampling | With `draws = 4000`, `sig_bhat_star / sig_bhat_asy` within 0.95–1.05. |
| Estimates; reduces to known | As W3. |

### API — the only differences allowed

1. One export added: `ps_weights` (and the S3 method `print.ps_weights`).
2. `weight_variance = c("known", "estimated")`, `weight_model = NULL` appended, in that order, to `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()`, `plotKM.band_subgroups()`.
3. `plotKM.band_subgroups()`: default of `draws.band`, 20 → 1000.

### Golden whitelist for this task

A field of a call in the baseline may move only as follows.

- **Rd** — in a call whose weight column contains a zero: variances, standard errors, z-scores, p-value texts, confidence limits, bands and resampled quantities.
- **Rb** — in a call whose positive weights are not all equal: `cox_results$score`, `logrank_results$chisq`, `$pval`, `$logrank_text`, and the output of `check_results()`.
- **Rc** — in a `plotKM.band_subgroups()` call that does not supply `draws.band`: the resampled matrices, the band critical value and band limits.
- **Ra, W1–W4** — nothing.

New elements allowed, and no others: `df_counting()`: `ess`, `weight_variance`, `cox_results$se`, `logrank_results$method`; `KM_diff()` (and so `plotKM.band_subgroups()`), `cumulative_rmst_bands()`, `cox_rhogamma()`: `weight_variance`.

Must not move anywhere: every point estimate, count, risk table and time weight, and every field of a call outside the three classes above.

---

## 6. M7 — Monte Carlo with the package code (compute; only if `M7: GO`)

Package at HEAD. Use the data-generating model, truth, time grid, exact calls and summary definitions of `dev/audit/oc/` (copy `oc_common.R`, `oc_sim.R`, `oc_summarise.R`, `oc_boot.R` into `stage2e/` and adapt the copies; do not edit files under `oc/`). Per-replicate results go to `stage2e/raw/`, ignored by git as `oc/raw/` is.

One generator is added:

```r
gen_v <- function(N, sc) {            # shaped like the vignette examples: about 14% treated, treated-arm ESS about a third of the arm
  X <- rnorm(N); H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))
  e <- plogis(-2.3 + 1.3 * X); A <- rbinom(N, 1, e)
  T <- ifelse(A == 1, LamInv(H, 1, sc), LamInv(H, 0, sc)); C <- runif(N, 0, 72)
  data.frame(X = X, e = e, treat = A, time = pmin(T, C), status = as.numeric(T <= C))
}
```

The propensity model is always `treat ~ X`, fitted on the replicate's data with `ps_weights()`. Seed of replicate `r` of cell `k`: `set.seed(1e6 * k + r)` (cells 31–33 take the data of the evaluation's cells 8 and 15 with their seeds).

| k | Cell | Data | Weights | Replicates | Evaluated |
|---|---|---|---|---|---|
| 7, 8, 9, 15, 19 | E-200, E-400, E-800, E-400-null, E-400-delay | the evaluation's cells | `"ipw"` (from the cell's own `glm`, `ps_weights(model = )`) | 1–500 | reproduction only (required, below) |
| 31 | Es-400 | data of cell 8 | `"stabilized"` | 2,000 | all |
| 32 | Es-400-null | data of cell 15 | `"stabilized"` | 2,000 | all |
| 33 | Es-400-null-strata | data of cell 15, `s = ifelse(X > 0, "hi", "lo")` as `strata.name` | `"stabilized"` | 2,000 | stratified log-rank only |
| 34 | Et-400-null | `gen(400, "null", "conf")`, `time <- ceiling(time)` | `"ipw"` | 2,000 | all except RMST |
| 35 | Fe-400 | `gen(400, "ph", "strong")` | `"ipw"` | 2,000 | all |
| 36 | Fe-800 | `gen(800, "ph", "strong")` | `"ipw"` | 2,000 | all |
| 37 | Fe-400-null | `gen(400, "null", "strong")` | `"ipw"` | 2,000 | all |
| 38 | V-1450-null | `gen_v(1450, "null")` | `"stabilized"` | 2,000 | all, FH(0,0) only |
| 39 | V-1450-ph | `gen_v(1450, "ph")` | `"stabilized"` | 2,000 | all, FH(0,0) only |
| 40 | Sub-800-null | `gen(800, "null", "conf")`; the analysis is of the subgroup `X > 0` | see below | 2,000 | all except bands and RMST |

**Per replicate**, for each of `"known"` and `"estimated"` (same data, same weight column): the evaluation's `KM_diff()` call (`draws.band = 500`, `seedstart = 8316951 + r`), `cumulative_rmst_bands()` (`draws_sb = 500`), and for FH(0,0) and FH(0,1) `df_counting()` then `cox_rhogamma(draws = 0)`. Record the effective-sample-size ratio by arm from `df_counting()$ess`. Cells 31 and 32 also: `wlr_dhat_estimates(dc, tzero = 18)` for FH(0,0) under `"estimated"`.

**Cell 40** has three analyses of the subgroup rows: (A) weights from the model fitted on all 800 rows, `"known"` — the route §2 prescribes for subgroups; (B) the model refitted on the subgroup rows, `"estimated"`; (C) the refitted weights, `"known"`. Truth: zero difference, log hazard ratio 0.

**Bootstrap benchmark** (cells 38 and 39, replicates 1–500, 200 resamples of subjects): as `oc_boot.R`, the propensity model refitted and the stabilised weights rebuilt in every resample; survival difference at 18, RMST difference to 36, log hazard ratio; normal interval from the bootstrap SD, and percentile interval.

**Pilot.** Before the full run, time 10 replicates of cells 31, 36, 38 and 40 on one core and print the projected wall-clock for the cores to be used. If it passes 120 min, do not start; record and go to Step V.

**Required** (write one line per requirement, `OK` or `FAIL` with the count of replicates failing and the largest difference, to `stage2e/m7_required.txt`):

1. Cells 7, 8, 9, 15, 19, replicates 1–500: with `"estimated"`, the package's `sig2_surv0`, `sig2_surv1`, `sig2_dhat` at the evaluation's time points, the band critical value and limits, and the RMST variance, critical value and band limits equal the script's (`ew_km_all()` in `oc_sim.R`, computed as that script computes them) within 1e-8 relative; `fit$sig_bhat_asy` and `abs(z.score)` for FH(0,0) and FH(0,1) equal the script's `se_ew` and `abs(z_ew)` within 1e-8, in every replicate where the script's alignment confirmation passes.
2. Cells 31 and 32, every replicate: `sig2_dhat`, the band critical value and the RMST variance under `"stabilized"` equal those under `"ipw"` from the same model (same seeds), within 1e-8, under both options.
3. Every cell, every replicate: `dhat`, `rmst_time`, `bhat` and `lr` are identical under `"known"` and `"estimated"`.
4. Cell 34, every replicate: with `"estimated"`, FH(0,0) `z.score^2` equals the `coxph()` score-residual route (Breslow) within 1e-8 and `fit$sig_bhat_asy` the `dfbeta` route (Breslow) within 1e-6.

For orientation only — a chat pilot with script code and other seeds, 2,000 replicates, known → estimated: ordinary propensity weights under the null (as cell 32), rejection 3.6% → 6.4%, stratified 5.3% → 5.9%; under proportional hazards (as cell 31), log hazard ratio coverage 0.976 → 0.943; the vignette-shaped null (as cell 38), rejection 8.8% → 10.2%, log hazard ratio coverage 0.919 → 0.898, survival difference at 36 months 0.880 → 0.874.

Everything else is evidence. Report, per cell and option, as the evaluation defined them: coverage and mean length for the survival difference at month 18, the band on [3, 42], RMST difference at 36 and the log hazard ratio (FH(0,0), FH(0,1)); SE/SD for the log hazard ratio; rejection rate of the two-sided 5% log-rank test in null cells; mean effective-sample-size ratio by arm; failures (including `ps_weights()` errors). For cells 31 and 32, the mean estimated correlation of the log-rank score with the difference at 18 against the empirical one. The full time profile (11 months) goes to `stage2e/m7_summary.csv`. Monte Carlo error: coverage ±0.005 at 2,000 and ±0.010 at 500.

---

## 7. Report format — `dev/audit/REPORT_km_fixes_stage2e.md`

Tables first, then a short plain reading.

- **Top line.** Branch, SHAs, gates, items DONE / NOT DONE, whether M7 ran, its pilot projection and wall-clock.
- **A. Items.** `| ID | status | commit | files touched | golden fields that moved | tests added or replaced |`
- **B. In numbers.** (1) Riders, on the golden datasets: the zero-weight cases before → after beside the positive-weight subset; `cox_results$score`, `logrank_results$chisq` and `check_results()` before → after on the weighted datasets; `ess` by arm. (2) The option, on one propensity-weighted dataset with tied times and one without: `"known"` → `"estimated"` for the log-rank z, the `cox_rhogamma()` SE, the SE of the survival difference at the median time, the band critical value, with the `coxph()`-route value beside each FH(0,0) number.
- **C. M7.** The `m7_required.txt` lines. Then **one** table: rows = cell × option (and the two bootstrap rows in cells 38, 39; analyses A, B, C in cell 40); columns = ESS ratio control / treated, coverage of the difference at 18, of the band, of RMST at 36, of the log hazard ratio FH(0,0) and FH(0,1), SE/SD FH(0,0), rejection rate FH(0,0) and FH(0,1), and mean length under `"estimated"` over mean length under `"known"` for the difference at 18, RMST at 36 and the log hazard ratio. `—` where not evaluated. At most ten lines of reading.
- **D. Vignette impact inventory**, in the five parts of Step V.
- **E. Recorded, not changed.** With all weights equal to one constant other than 1 the equal-weight formulas are unchanged: show `z.score` and the `cox_rhogamma()` SE on one dataset with no weights, with weights all 1 and with weights all 2. Anything else noticed and left alone.
- **F. Seed probe, pressure test, API, CRAN check, forestsearch.**
- **G. Gate G2 output.**

Findings only. No recommendation and no follow-up task proposals.

---

## 8. Gate G2 — post-conditions

Run from the repo root with `START_SHA` and `MASTER_SHA` set, after the last item commit. Every line must print `OK`. A check for an item recorded as NOT DONE is not run, and the report says so: check 4 is Rd, check 5 is Rb and Ra, checks 6–9 are W1–W4; if M7 did not run, check 9 is not run. Checks 10 and 11 write two files: paste the output into the report, commit, then run check 13.

```bash
S=dev/audit/stage2e

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. API: the listed differences and nothing else; vignettes and README untouched
Rscript $S/api.R compare && echo "OK 2a" || echo "FAIL: API"
git diff --quiet "$START_SHA" HEAD -- vignettes/ README.md && echo "OK 2b" || echo "FAIL: vignettes or README changed"

# 3. golden: identical outside the Stage-2e whitelist; new elements only those listed
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. Rd: zero-weight rows are inert
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
set.seed(3); n <- 300; x <- rnorm(n); a <- rbinom(n, 1, plogis(-0.3 + x))
T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60); keep <- runif(n) < 0.7
for (ties in c(FALSE, TRUE)) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, cc)) else pmin(T, cc), status = as.numeric(T <= cc), treat = a, s = ifelse(x > 0, "hi", "lo"))
  d$w01 <- as.numeric(keep); d$wu <- ifelse(keep, ifelse(a == 1, 1 / plogis(-0.3 + x), 1 / (1 - plogis(-0.3 + x))), 0)
  sub <- d[keep, ]
  for (wn in c("w01", "wu")) {
    f <- q(df_counting(d, "time", "status", "treat", weight.name = wn, strata.name = "s"))
    p <- q(df_counting(sub, "time", "status", "treat", weight.name = if (wn == "w01") NULL else wn, strata.name = "s"))
    cf <- q(cox_rhogamma(f, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5)))
    cp <- q(cox_rhogamma(p, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5)))
    kf <- q(KM_diff(d, "time", "status", "treat", weight.name = wn, at_points = p$at_points))
    kp <- q(KM_diff(sub, "time", "status", "treat", weight.name = if (wn == "w01") NULL else wn, at_points = p$at_points))
    j <- match(p$at_points, f$at_points)
    stopifnot(!anyNA(j), max(abs(f$sig2_surv1[j] - p$sig2_surv1)) < 1e-12, max(abs(f$sig2_surv0[j] - p$sig2_surv0)) < 1e-12,
              abs(f$z.score - p$z.score) < 1e-10, abs(f$z.score_stratified - p$z.score_stratified) < 1e-10,
              abs(cf$fit$sig_bhat_asy / cp$fit$sig_bhat_asy - 1) < 1e-10, abs(cf$z.score - cp$z.score) < 1e-10,
              max(abs(kf$sig2_dhat - kp$sig2_dhat)) < 1e-12)
  }
}
cat("OK 4\n")'

# 5. Rb, Ra: the survival-based outputs with weights; effective sample size
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
set.seed(4); n <- 300; x <- rnorm(n); e <- plogis(-0.5 + x); a <- rbinom(n, 1, e)
T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60)
ess <- function(w) sum(w)^2 / sum(w^2)
for (ties in c(FALSE, TRUE)) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, cc)) else pmin(T, cc), status = as.numeric(T <= cc), treat = a,
                  w = ifelse(a == 1, 1 / e, 1 / (1 - e)), s = ifelse(x > 0, "hi", "lo"), two = 2)
  dc  <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
  dcs <- q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s"))
  u   <- q(df_counting(d, "time", "status", "treat")); c2 <- q(df_counting(d, "time", "status", "treat", weight.name = "two"))
  cpb <- coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
  cps <- coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, robust = TRUE, ties = "breslow")
  cpe <- summary(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE))
  stopifnot(abs(dc$logrank_results$chisq / cpb$rscore - 1) < 1e-8, abs(dc$logrank_results$chisq / dc$z.score^2 - 1) < 1e-8,
            abs(dcs$logrank_results$chisq / cps$rscore - 1) < 1e-8, abs(dcs$logrank_results$chisq / dcs$z.score_stratified^2 - 1) < 1e-8,
            abs(dc$cox_results$score / cpe$robscore[1] - 1) < 1e-8,
            abs(u$logrank_results$chisq / survdiff(Surv(time, status) ~ treat, data = d)$chisq - 1) < 1e-10,
            identical(u$logrank_results$chisq, c2$logrank_results$chisq),
            abs(dc$ess$ess[1] / ess(d$w[d$treat == 1]) - 1) < 1e-10, abs(dc$ess$ess[2] / ess(d$w[d$treat == 0]) - 1) < 1e-10,
            all(abs(u$ess$ess - u$ess$n) < 1e-9), all(abs(c2$ess$ess - c2$ess$n) < 1e-9), all(abs(dc$ess$ratio - dc$ess$ess / dc$ess$n) < 1e-12))
}
cat("OK 5\n")'

# 6. W2: Kaplan-Meier with estimated weights against the n x m reference; weight form; reduction to known; misuse; resampling
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
set.seed(11); n <- 400; x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); A <- rbinom(n, 1, plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
T <- rexp(n, 0.04 * exp(0.5 * x1 - 0.4 * A)); C <- runif(n, 0, 72)
# n_g x m matrix of Kaplan-Meier influence values of one group on a sorted grid holding its event times
kmpsi <- function(U, D, W, g) {
  ybar <- colSums(outer(U, g, ">=") * W); dN <- diff(c(0, colSums(outer(U[D == 1], g, "<=") * W[D == 1])))
  ok <- ybar > 0 & (ybar - dN) > 1e-10 * ybar
  t(apply(sweep(outer(U, g, "==") * D, 2, ifelse(ok, 1 / (ybar - dN), 0), "*") - sweep(outer(U, g, ">="), 2, ifelse(ok, dN / (ybar * (ybar - dN)), 0), "*"), 1, cumsum))
}
for (ties in c(FALSE, TRUE)) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, C)) else pmin(T, C), status = as.numeric(T <= C), treat = A, x1 = x1, x2 = x2)
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = "ipw"); pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized")
  d$w <- pw$weights; d$wst <- pst$weights; g <- sort(unique(d$time))
  fk <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 300))
  fe <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 300, weight_variance = "estimated", weight_model = pw))
  fs <- q(KM_diff(d, "time", "status", "treat", weight.name = "wst", at_points = g, draws.band = 300, weight_variance = "estimated", weight_model = pst))
  pz <- pw; pz$IF <- pw$IF * 0
  f0 <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 300, weight_variance = "estimated", weight_model = pz))
  de <- q(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw)); jd <- match(g, de$at_points)
  # reference: Phi_a[i, j] = 1(i in arm a) w_i psi_i(g_j) + IF_i (sum_{k in a} psi_k(g_j) dw_k)
  Phi <- list()
  for (a in c(0, 1)) { ia <- which(d$treat == a); P <- kmpsi(d$time[ia], d$status[ia], d$w[ia], g)
    Fa <- pw$IF %*% t(crossprod(P, pw$dw[ia, , drop = FALSE])); Fa[ia, ] <- Fa[ia, ] + d$w[ia] * P; Phi[[a + 1]] <- Fa }
  v0 <- fe$surv0^2 * colSums(Phi[[1]]^2); v1 <- fe$surv1^2 * colSums(Phi[[2]]^2)
  vd <- colSums((Phi[[2]] * rep(fe$surv1, each = n) - Phi[[1]] * rep(fe$surv0, each = n))^2)
  stopifnot(identical(fe$dhat, fk$dhat), identical(fe$surv0, fk$surv0), identical(fe$surv1, fk$surv1),
            rel(fe$sig2_surv0, v0) < 1e-8, rel(fe$sig2_surv1, v1) < 1e-8, rel(fe$sig2_dhat, vd) < 1e-8,
            rel(de$sig2_surv0[jd], v0) < 1e-8, rel(de$sig2_surv1[jd], v1) < 1e-8,
            rel(fs$sig2_dhat, fe$sig2_dhat) < 1e-8, rel(fs$dhat_star, fe$dhat_star) < 1e-8, abs(fs$c_alpha_band / fe$c_alpha_band - 1) < 1e-8,
            rel(f0$sig2_dhat, fk$sig2_dhat) < 1e-10, rel(f0$dhat_star, fk$dhat_star) < 1e-10,
            identical(fe$weight_variance, "estimated"), identical(fk$weight_variance, "known"), identical(de$weight_variance, "estimated"),
            err(KM_diff(d, "time", "status", "treat", weight_variance = "estimated", weight_model = pw)),
            err(KM_diff(d, "time", "status", "treat", weight.name = "wst", weight_variance = "estimated", weight_model = pw)),
            err(KM_diff(d[-1, ], "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw)),
            err(KM_diff(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw$model)),
            err(KM_diff(d, "time", "status", "treat", weight.name = "w", draws = 50, weight_variance = "estimated", weight_model = pw)),
            err(df_counting(d, "time", "status", "treat", weight.name = "wst", weight_variance = "estimated", weight_model = pw)))
  if (!ties) {
    f4 <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 4000, weight_variance = "estimated", weight_model = pw))
    pos <- f4$sig2_dhat > 0; r <- median(apply(f4$dhat_star[pos, ], 1, var) / f4$sig2_dhat[pos])
    rm <- q(cumulative_rmst_bands(d, f4, "time", "status", "treat", weight.name = "w", draws_sb = 300, plot = FALSE, weight_variance = "estimated", weight_model = pw))
    stopifnot(r > 0.95, r < 1.05, identical(rm$weight_variance, "estimated"),
              err(cumulative_rmst_bands(d, fk, "time", "status", "treat", weight.name = "w", draws_sb = 300, plot = FALSE, weight_variance = "estimated", weight_model = pw)))
    cat(sprintf("resampled / analytic variance of the difference (median over the grid): %.3f\n", r))
  }
}
cat("OK 6\n")'

# 7. W3, W4: log-rank and Cox against the coxph residual routes (ties, strata); cumulative; resampling
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
set.seed(11); n <- 400; x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); A <- rbinom(n, 1, plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
T <- rexp(n, 0.04 * exp(0.5 * x1 - 0.4 * A)); C <- runif(n, 0, 72)
P00 <- list(rho = 0, gamma = 0)
for (ties in c(FALSE, TRUE)) for (type in c("ipw", "stabilized")) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, C)) else pmin(T, C), status = as.numeric(T <= C), treat = A, x1 = x1, x2 = x2, s = ifelse(x2 == 1, "hi", "lo"))
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = type); d$w <- pw$weights
  adj <- function(r) { r <- as.numeric(r); as.numeric(r + pw$IF %*% colSums((r / d$w) * pw$dw)) }
  sres <- function(f) { r <- as.numeric(residuals(f, type = "score", weighted = TRUE)); sum(r)^2 / sum(adj(r)^2) }
  chi  <- sres(coxph(Surv(time, status) ~ treat, data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
  chis <- sres(coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
  se_b <- sqrt(sum(adj(residuals(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow"), type = "dfbeta", weighted = TRUE))^2))
  se_e <- sqrt(sum(adj(residuals(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE), type = "dfbeta", weighted = TRUE))^2))
  dk  <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
  de  <- q(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw))
  des <- q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s", weight_variance = "estimated", weight_model = pw))
  ck  <- q(cox_rhogamma(dk, scheme = "fh", scheme_params = P00)); cx <- q(cox_rhogamma(de, scheme = "fh", scheme_params = P00))
  one <- q(wlr_cumulative(de, scheme = "fh", scheme_params = P00)); cum <- q(wlr_cumulative(de, scheme = "fh", scheme_params = P00, return_cumulative = TRUE))
  t0 <- median(d$time); j0 <- max(which(de$at_points <= t0)); dcen <- d; dcen$status <- d$status * (d$time <= t0); dcen$time <- pmin(d$time, t0)
  dn <- q(df_counting(dcen, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw))
  stopifnot(identical(de$lr, dk$lr), identical(cx$fit$bhat, ck$fit$bhat), identical(de$cox_results$hr, dk$cox_results$hr),
            abs(de$z.score^2 / chi - 1) < 1e-8, abs(des$z.score_stratified^2 / chis - 1) < 1e-8,
            abs(de$logrank_results$chisq / chi - 1) < 1e-8, abs(des$logrank_results$chisq / chis - 1) < 1e-8,
            abs(cx$fit$z.score^2 / chi - 1) < 1e-8, abs(cx$fit$sig_bhat_asy / se_b - 1) < 1e-6, abs(de$cox_results$se / se_e - 1) < 1e-8,
            abs(one$sig2.score / de$sig2_lr - 1) < 1e-10, abs(cum$sig2.score[length(cum$sig2.score)] / de$sig2_lr - 1) < 1e-10,
            abs(cum$sig2.score[j0] / dn$sig2_lr - 1) < 1e-8,
            identical(de$weight_variance, "estimated"), identical(cx$weight_variance, "estimated"), identical(ck$weight_variance, "known"))
  if (!ties && type == "ipw") {
    cs <- q(cox_rhogamma(de, scheme = "fh", scheme_params = P00, draws = 4000)); r <- cs$fit$sig_bhat_star / cs$fit$sig_bhat_asy
    stopifnot(is.finite(r), r > 0.95, r < 1.05, cx$fit$sig_bhat_asy < ck$fit$sig_bhat_asy)
    cat(sprintf("Cox: resampling SE / analytic SE %.3f; analytic SE estimated %.4f, known %.4f\n", r, cx$fit$sig_bhat_asy, ck$fit$sig_bhat_asy))
  }
}
cat("OK 7\n")'

# 8. W1: ps_weights() against numerical derivatives
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
set.seed(11); n <- 400; x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); A <- rbinom(n, 1, plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
d <- data.frame(treat = A, x1 = x1, x2 = x2)
g <- glm(treat ~ x1 + x2, family = binomial(), data = d); e <- as.numeric(fitted(g)); p1 <- mean(A); xm <- model.matrix(g)
pw <- ps_weights(treat ~ x1 + x2, data = d, type = "ipw"); pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized"); pm <- ps_weights(model = g, type = "ipw")
# numerical derivatives, from fits iterated to convergence
ctl <- glm.control(epsilon = 1e-13, maxit = 100); b <- coef(glm(treat ~ x1 + x2, family = binomial(), data = d, control = ctl))
wof <- function(bb) { pp <- as.numeric(plogis(xm %*% bb)); ifelse(A == 1, 1 / pp, 1 / (1 - pp)) }
num <- sapply(seq_along(b), function(j) { h <- 1e-6; bp <- b; bp[j] <- b[j] + h; bm <- b; bm[j] <- b[j] - h; (wof(bp) - wof(bm)) / (2 * h) })
ii <- c(1, 7, 50, 200, 399)
ifn <- t(sapply(ii, function(i) { wt <- rep(1, n); wt[i] <- 1 + 1e-4
  (coef(glm(treat ~ x1 + x2, family = quasibinomial(), data = d, weights = wt, control = ctl)) - b) / 1e-4 }))
stopifnot(inherits(pw, "ps_weights"), all(c("weights", "type", "treat", "IF", "dw", "model") %in% names(pw)),
          max(abs(pw$weights - ifelse(A == 1, 1 / e, 1 / (1 - e)))) < 1e-10, max(abs(pst$weights - ifelse(A == 1, p1 / e, (1 - p1) / (1 - e)))) < 1e-10,
          max(abs(pm$weights - pw$weights)) < 1e-10, max(abs(pm$IF - pw$IF)) < 1e-10, max(abs(pm$dw - pw$dw)) < 1e-10,
          all(dim(pw$IF) == c(n, 3)), all(dim(pw$dw) == c(n, 3)), all(pw$treat == A),
          max(abs(pw$dw - num)) / max(abs(num)) < 1e-5, max(abs(pw$IF[ii, ] - ifn)) / max(abs(ifn)) < 1e-3,
          max(abs(pst$dw - ifelse(A == 1, p1, 1 - p1) * pw$dw)) < 1e-12, max(abs(pst$IF - pw$IF)) < 1e-12,
          err(ps_weights(model = glm(treat ~ x1, family = binomial(link = "probit"), data = d))), err(ps_weights(model = lm(treat ~ x1, data = d))),
          err(ps_weights(treat ~ x1, data = transform(d, x1 = ifelse(seq_len(n) == 3, NA, x1)))))
invisible(utils::capture.output(print(pw)))
cat("OK 8\n")'

# 9. M7 requirements
[ -f $S/m7_required.txt ] && [ "$(grep -c '^OK' $S/m7_required.txt)" = "4" ] && ! grep -q FAIL $S/m7_required.txt && echo "OK 9" || echo "FAIL: M7 requirements"

# 10. seed probe: statuses as after Stage 2c
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2e.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2e.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 10a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2e.txt)" = "3" ] && echo "OK 10b" || echo "FAIL: controls"

# 11. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2e.txt 2>&1 && echo "OK 11a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2e.log | grep -Evq 'ERROR|WARNING' && echo "OK 11b" || echo "FAIL: R CMD check"

# 12. deliverables and record
for f in golden_s2d.rds api_s2d.rds probe_2e.txt 00check_2e.log; do [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2e.md ] && echo "OK 12 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD
```

```bash
# 13. after the final commit
[ -z "$(git status --porcelain)" ] && echo "OK 13" || echo "FAIL: uncommitted changes"
```

Notes:
- Check 2a: `stage2e/api.R compare` exits 0 only when the differences from `api_s2d.rds` are exactly those of §5 "API" (those of the items that are DONE).
- Check 3: `golden.R compare final` uses the union of the whitelists of the items that are DONE.
- Checks 4–8 were run in chat, as extracted from this document, against a stand-in (the branch with the seven Rd tests changed, plus wrappers applying §2 with the existing helpers); every line passed. On the branch as it stands, check 4 fails at the log-rank z (0/1 weights against the subset: a difference of 0.018), check 5 at its first comparison (`logrank_results$chisq` is the unweighted `survdiff` value), and checks 6–8 at the first call (`ps_weights()` does not exist).
- Check 6 prints the ratio of the resampled to the analytic variance of the difference; with the multipliers mapped as §2 says it is 0.993. Check 7 prints the resampling SE over the analytic SE of the log hazard ratio: 1.002.
