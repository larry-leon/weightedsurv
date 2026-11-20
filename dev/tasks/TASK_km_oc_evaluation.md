# TASK — weightedsurv 0.2.0: operating characteristics of weighted inference (bias, coverage, interval length)

Written in chat 2026-11-19 (revision 2: replaces the file of the same name written earlier that day, which was not run).
Executed by Claude Code in the local `weightedsurv` repo, on branch `km-fixes-0.2.0`.

| | |
|---|---|
| Touches `R/` | **No.** Evaluation only. |
| Writes | Only under `dev/tasks/` and `dev/audit/` (new folder `dev/audit/oc/`). |
| Depends on | Stage 2c (R1 and R2 committed; `dev/audit/REPORT_km_fixes_stage2c.md` in the tree). |
| Branch | `km-fixes-0.2.0`. No merge. **Never push.** |
| Compute | **Main run**: 19 cells × 2,000 replicates × 2 code versions, 500 multiplier draws; expected 20–45 min on 12 cores (scaled from the M3 and M4 timings). **Part E**: bootstrap benchmark, 3 cells × 500 replicates × 200 resamples; expected 5–10 min. Each runs only if the kickoff says `GO`. Stop and record if the projection for the main run passes 120 min. |
| Gates | Stop-on-failure, not stop-to-ask. A part that cannot meet its gate is recorded as NOT DONE and the task continues, where §11 says so. |

---

## 1. Why

One study with known truth, measuring bias, coverage and interval length across the range of event times, for:

- the released package (0.1.0) against the branch after Stage 2c (R1 analytic variance, R2 resampling);
- the matching variance for the subject-weighted log-rank test and Cox estimate, which is **not** in the package and is not authorised — coded in the evaluation script;
- a version that carries the **estimation of the weights** into the variance and the resampling (the perturbed-nuisance device of León, Cai & Wei 2009, eqs. 3.8–3.13) — also script-side only;
- the path `plotKM.band_subgroups()` takes by default with weights (`draws = 20`, a 20-draw pointwise SE), which the Stage 2c vignette inventory showed is what the vignette band figures use.

Settings where nothing should change (no weights; weights unrelated to outcome) and settings where the proposal is expected to fall short (extreme weights; small samples; the far tail) are included on purpose. The report states findings and recommends nothing.

Not repeated: identities already established (R1 equals `survfit(robust = TRUE)`; unweighted output unchanged by Stage 2c) and the Stage 2b-versus-2c comparisons of M2b and M3.

---

## 2. What is compared

| Name | What it is | Where it runs |
|---|---|---|
| `cran` | weightedsurv 0.1.0 as released: a git worktree of commit `166919e` (`MASTER_SHA` of Stage 2c; `R/` byte-identical to CRAN 0.1.0), loaded with `pkgload::load_all()` | "cran" process |
| `proposed` | Kaplan–Meier quantities: the branch at HEAD, `draws = 0`. Log-rank and Cox: the score-residual variance of §6, coded in the script. Weights treated as known | "proposed" process |
| `proposed_rect` | RMST only: HEAD with `mean_rule = "rectangle"` (`proposed` uses the default, trapezoid) | "proposed" process |
| `proposed_ew` | cells with estimated weights only (7, 8, 9, 15, 19): `proposed` with the weight model carried through, §7; script-side | "proposed" process |
| `proposed_d20` | cells 5 and 8 only: HEAD with `draws = 20` — the default path of `plotKM.band_subgroups()` | "proposed" process |
| `cran_mean1` | Part C only: `cran` with weights rescaled to mean 1 (`w * n / sum(w)`); FH(0,0) | "cran" process |
| `cran_star` | Part C only: the package's resampling SE, `cox_rhogamma(draws = 200)`; first 500 replicates, FH(0,0) | "cran" process |
| `boot`, `boot_pct` | Part E only: nonparametric bootstrap, normal and percentile intervals | bootstrap script |

One R process per code version. Both generate identical data from the seeds in §3; a checksum per replicate proves it.

Calls use the same arguments under both versions, except arguments 0.1.0 does not have. Read each version's signatures and return values from its source before writing the calls, and list in the report the exact calls and the returned fields used.

---

## 3. Data-generating model

Time is in months. Marginal survival in each arm is exact by construction, so every truth is closed-form.

```r
# cumulative hazard of the potential outcome T^a, and its inverse; sc = "ph", "null" or "delay"
Lam    <- function(t, a, sc) if (a == 0 || sc == "null") 0.04 * t else
                             if (sc == "ph") 0.04 * exp(-0.4) * t else 0.04 * pmin(t, 6) + 0.02 * pmax(t - 6, 0)
LamInv <- function(h, a, sc) if (a == 0 || sc == "null") h / 0.04 else
                             if (sc == "ph") h / (0.04 * exp(-0.4)) else ifelse(h <= 0.24, h / 0.04, 6 + (h - 0.24) / 0.02)

gen <- function(N, sc, assign) {                       # assign = "rand", "conf" or "strong"
  X <- rnorm(N)
  H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))           # standard exponential, larger X -> earlier event
  e <- switch(assign, rand = rep(0.5, N), conf = plogis(-0.5 + X), strong = plogis(-0.5 + 2 * X))
  A <- rbinom(N, 1, e)
  T <- ifelse(A == 1, LamInv(H, 1, sc), LamInv(H, 0, sc))
  C <- runif(N, 0, 72)
  data.frame(X = X, e = e, treat = A, time = pmin(T, C), status = as.numeric(T <= C))
}
```

- `X` is prognostic and, under "conf" and "strong", drives treatment. Marginally `P(T^a > t) = exp(-Lam(t, a, sc))`; control median 17.3 months.
- Scenarios: **ph** — hazard ratio `exp(-0.4)` throughout (true log HR −0.4); **null** — no effect; **delay** — no effect to month 6, hazard ratio 0.5 after.
- About 37% censoring. Under "conf" about 10% of each arm is still at risk at month 42 and 6% at month 48.

### Time points

`tt = c(1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48)`. Every pointwise quantity is evaluated at all eleven. Month 48 is the far tail by design. Headline times for the summary tables: 6, 18, 36. Band range for restricted coverage: `[3, 42]`.

A replicate is **evaluable at t** if each arm has at least one subject with observed time `>= t`. Metrics at t use evaluable replicates only; the number excluded is reported.

### Weights

| Code | Assignment | Weights `w` |
|---|---|---|
| U | rand | none (`weight.name = NULL`) |
| R | rand | `runif(N, 0.5, 2)`, drawn after `gen()` — unrelated to outcome |
| S | rand | `exp(0.5 * X - 0.125)` — standardises to a population with `X ~ N(0.5, 1)` |
| D | conf | true propensity weights: `ifelse(treat == 1, 1 / e, 1 / (1 - e))` |
| E | conf | the same with `e` replaced by fitted values of `g <- glm(treat ~ X, binomial)`, refitted in each replicate |
| F | strong | true propensity weights |
| Dt | conf | as D, with observed times rounded up to whole months (`ceiling(time)`) — tied times |

### Cells

| k | Cell | Scenario | Weights | N | Not evaluated |
|---|---|---|---|---|---|
| 1 | U-400 | ph | U | 400 | |
| 2 | R-400 | ph | R | 400 | |
| 3 | S-400 | ph | S | 400 | log HR bias and coverage |
| 4, 5, 6 | D-200, D-400, D-800 | ph | D | 200, 400, 800 | |
| 7, 8, 9 | E-200, E-400, E-800 | ph | E | 200, 400, 800 | |
| 10, 11 | F-400, F-800 | ph | F | 400, 800 | |
| 12 | Dt-400 | ph | Dt | 400 | RMST; all of Part C |
| 13 | U-400-null | null | U | 400 | |
| 14 | D-400-null | null | D | 400 | |
| 15 | E-400-null | null | E | 400 | |
| 16 | F-400-null | null | F | 400 | |
| 17 | U-400-delay | delay | U | 400 | log HR bias and coverage |
| 18 | D-400-delay | delay | D | 400 | log HR bias and coverage |
| 19 | E-400-delay | delay | E | 400 | log HR bias and coverage |

Where log HR bias and coverage are not evaluated (no single true hazard ratio), SE/SD and interval length are still recorded.

### Seeds

Replicate `r` of cell `k`: `set.seed(1e6 * k + r)` immediately before `gen()`; data and weights are complete before any package call. Multiplier draws: `seedstart = 8316951 + r` in `KM_diff()`. Before `cox_rhogamma(draws = 200)`: `set.seed(3e7 + 1e6 * k + r)`. Part E: `set.seed(6e7 + 1e6 * k + r)`. (0.1.0 resets the global RNG inside some calls; this ordering makes that harmless.)

Record per replicate the checksum `sprintf("%.10f", sum(time) + sum(status) + sum(w))` (`w = 1` in U cells).

### Truth

- `S_a(t) = exp(-Lam(t, a, sc))`. Cell 3 only: `S_a(t) = pnorm(qnorm(exp(-Lam(t, a, "ph"))) - 0.3)`.
- `delta(t) = S_1(t) - S_0(t)`; `R(t) = ` integral of `delta` over `[0, t]` by `stats::integrate` (`rel.tol = 1e-10`).
- Log hazard ratio: −0.4 in ph cells except cell 3; 0 in null cells.

---

## 4. Part A — pointwise (all 19 cells)

Per replicate and version:

`fit <- KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = <"w" or NULL>, alpha = 0.05, draws = 0, draws.band = 500, seedstart = 8316951 + r, show_resamples = FALSE)` on the default grid.

At every `t` in `tt`, reading each curve at the last grid point `<= t`:

| Quantity | Estimate | 95% interval |
|---|---|---|
| `S_0(t)`, `S_1(t)` | `surv0`, `surv1` | estimate `± qnorm(0.975) * sqrt(sig2_surv0 / sig2_surv1)` |
| `delta(t)` | `dhat` | the function's `lower`, `upper` |

`proposed_d20` (cells 5 and 8, HEAD): the same call with `draws = 20`; `delta(t)` intervals at every `t`, and the band of §5.

## 5. Part B — simultaneous bands and RMST

**Band for `delta`** (all 19 cells): `sb_lower`, `sb_upper` from the same `fit`. As in M3, a grid point counts only where the band has positive width (both limits finite, upper above lower). Coverage = the true `delta(t_j)` lies inside the band at every such grid point (i) in `[3, 42]`, and (ii) anywhere. Width = the mean of `sb_upper - sb_lower` over such points in `[3, 42]`. Also record `c_alpha_band`.

**RMST** (not cell 12): `cumulative_rmst_bands(df, fit, tte.name, event.name, treat.name, weight.name, draws_sb = 500, plot = FALSE)`; under HEAD a second call with `mean_rule = "rectangle"`.
- Pointwise at every `t` in `tt`, at the last grid point `t_j <= t`: estimate and the function's pointwise interval; truth `R(t_j)` at that same grid point.
- Simultaneous band: coverage of `R(t_j)` over grid points with positive width in `[3, 42]`, and anywhere; mean width over `[3, 42]`.

A replicate in which a call errors is caught, counted as a failure for that version, and counts as not covering. Non-finite limits are counted separately and count as not covering. Lengths and widths average finite values only.

## 6. Part C — log-rank and Cox (not cell 12)

Schemes FH(0,0) and FH(0,1). Per replicate:

`dc <- df_counting(df, tte.name, event.name, treat.name, weight.name, scheme = "fh", scheme_params = list(rho = 0, gamma = g))`, then `cx <- cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = g), draws = 0)`.

| Estimator | Log HR and interval | Log-rank z |
|---|---|---|
| `cran` | `cx$fit$bhat`; interval from `cx$hr_ci_asy` on the log scale; SE `cx$fit$sig_bhat_asy` | `dc$z.score` |
| `cran_mean1` | the same with rescaled weights; FH(0,0) only | the same |
| `cran_star` | `cox_rhogamma(draws = 200)`: SE `sig_bhat_star`, interval `hr_ci_star`; FH(0,0), first 500 replicates | — |
| `proposed` | `cx$fit$bhat` from the HEAD process; SE and z from `cox_robust()` below; interval `bhat ± qnorm(0.975) * se` | `cox_robust(cx$fit, 0)$z` |

The HEAD process also records the package's own Cox output; the report gives its maximum difference from `cran` (expected 0: that code is unchanged).

### The proposed variance

For subject `i` with time `U_i`, event `D_i`, arm `Z_i`, subject weight `a_i`, and time weight `K(s)`:

`psi_i(b) = a_i * [ D_i K(U_i) (Z_i − zbar(U_i))  −  exp(b Z_i) * sum over s <= U_i of K(s) (Z_i − zbar(s)) dL(s) ]`

with `S0(s) = sum_j a_j 1(U_j >= s) exp(b Z_j)`, `zbar(s) = sum_j a_j 1(U_j >= s) Z_j exp(b Z_j) / S0(s)`, `dL(s) = sum_j a_j D_j 1(U_j = s) / S0(s)`.

Score `U(b) = sum_i psi_i(b)`; information `I(b) = sum_s K(s) zbar(s) (1 − zbar(s)) sum_j a_j D_j 1(U_j = s)`; `Var(bhat) = sum_i psi_i(bhat)^2 / I(bhat)^2`; log-rank `z = U(0) / sqrt(sum_i psi_i(0)^2)`.

Reference implementation, checked in chat against 0.1.0 objects (root check 1e-10; information equal to the package's `i_bhat`; SE equal to `coxph(robust = TRUE)` to 1e-14 at FH(0,0)). Use it as given. It assumes no tied times, which holds wherever Part C runs.

```r
cox_robust <- function(fit, beta = fit$bhat) {         # fit = cox_rhogamma(...)$fit ; psi is returned in the order of fit$time
  o <- order(fit$time); D <- fit$delta[o]; Z <- fit$z[o]; a <- fit$w_hat[o]; K <- fit$wt_rg[o]
  eb <- exp(beta * Z); S0 <- rev(cumsum(rev(a * eb))); S1 <- rev(cumsum(rev(a * eb * Z))); zb <- S1 / S0
  dL <- a * D / S0; cK <- cumsum(K * dL); cKz <- cumsum(K * zb * dL)
  psi <- numeric(length(o)); psi[o] <- a * (D * K * (Z - zb) - eb * (Z * cK - cKz))
  U <- sum(a * D * K * (Z - zb)); info <- sum(a * D * K * zb * (1 - zb))
  list(U = U, info = info, se = sqrt(sum(psi^2)) / info, z = U / sqrt(sum(psi^2)), psi = psi)
}
```

The sign of `z` is opposite to the package's; only `|z|` is used.

---

## 7. Part D — estimated weights carried through (`proposed_ew`; cells 7, 8, 9, 15, 19; only if `Part D: GO`)

The weights in these cells come from the fitted model `g`. `proposed` treats them as known. `proposed_ew` adds, to each subject's contribution, the effect that subject has on the weights through the fitted coefficients. Checked in chat on this model (1,000 replicates, N = 400): SE/SD of the survival difference at months 3 to 42 is 1.06–1.16 with weights treated as known and 0.97–1.03 carried through; for the log HR 1.17 and 0.98.

```r
ew_parts <- function(g, treat) {                       # g = glm(treat ~ ..., binomial) fitted on the analysis rows
  xm <- model.matrix(g); e <- fitted(g)
  list(IF = (xm * (treat - e)) %*% solve(crossprod(xm * sqrt(e * (1 - e)))),   # n x p, influence of the coefficients
       dw = ifelse(treat == 1, -(1 - e) / e, e / (1 - e)) * xm)                # n x p, derivative of the weights
}
```

### Kaplan–Meier quantities

Let `P_a` be the `n_a × m` matrix of `psi_i(t_j)` for the subjects of arm `a` (rows in data order, `idx_a` their row numbers) over the grid of `fit` — the Stage 2c influence quantity. Take it from the package's internal helper (`ws_km_psi_parts()` and companions in `R/internal_helpers.R`) or from the `ij` code of `dev/audit/m2/km_m2_montecarlo.R`.

```r
ew_km <- function(P0, P1, S0, S1, idx0, idx1, w, IF, dw) {     # S0, S1 = surv0, surv1 on the grid
  n <- nrow(IF)
  Dg0 <- crossprod(P0, dw[idx0, , drop = FALSE]); Dg1 <- crossprod(P1, dw[idx1, , drop = FALSE])   # m x p
  Phi0 <- IF %*% t(Dg0); Phi0[idx0, ] <- Phi0[idx0, ] + w[idx0] * P0                                 # n x m
  Phi1 <- IF %*% t(Dg1); Phi1[idx1, ] <- Phi1[idx1, ] + w[idx1] * P1
  list(Dg0 = Dg0, Dg1 = Dg1, sig2_surv0 = S0^2 * colSums(Phi0^2), sig2_surv1 = S1^2 * colSums(Phi1^2),
       sig2_dhat = colSums((Phi1 * rep(S1, each = n) - Phi0 * rep(S0, each = n))^2))
}
```

- **Pointwise** (`S_0`, `S_1`, `delta` at every `t` in `tt`): the `proposed` estimates with these variances; intervals `± qnorm(0.975) * SE`.
- **Resampled curves.** With `G` the `n × B` matrix of the multipliers `KM_diff()` used, rows in data order:
  `surv_a_star_ew = surv_a_star − (S_a * Dg_a) %*% (t(IF) %*% G)`, and `dhat_star_ew = surv1_star_ew − surv0_star_ew`.
- **Band for `delta`.** Apply to `dhat_star_ew` and `sqrt(sig2_dhat)` from `ew_km()` exactly the steps `KM_diff()` applies to `dhat_star` and its SE (standardisation, time range, sup, quantile, limits); copy those lines from the HEAD source.
- **RMST, pointwise and band** (`mean_rule = "trapezoid"`). Apply to the adjusted resampled curves exactly the steps `cumulative_rmst_bands()` applies to the unadjusted ones, for both of its draw sets.

### Log-rank and Cox

With `ord <- order(df$time)` (confirm `cx$fit$time` equals `df$time[ord]` and `cx$fit$z` equals `df$treat[ord]`), `cr <- cox_robust(cx$fit)`:

`psi_ew = cr$psi + IF[ord, ] %*% colSums((cr$psi / cx$fit$w_hat) * dw[ord, ])`; `se_ew = sqrt(sum(psi_ew^2)) / cr$info`; the log-rank z likewise from `cox_robust(cx$fit, 0)`.

For FH(0,1) the dependence of the time weights on the estimated weights is ignored; say so in the report.

---

## 8. Part E — bootstrap benchmark (only if `Part E: GO`)

Cells 5, 8 and 10; the first 500 replicates; 200 resamples of subjects with replacement. In cell 8 the propensity model is refitted in every resample; in cells 5 and 10 the true weights travel with the subjects.

Computed on the full data and on each resample, with `survival` only:
- `delta(18)`: difference of weighted `survfit` curves;
- `R(36)`: difference of the weighted restricted means to 36 (`rmean = 36`);
- log HR: `coef(coxph(Surv(time, status) ~ treat, weights = w))`.

`boot`: full-data estimate `± qnorm(0.975) *` SD over resamples. `boot_pct`: the 2.5% and 97.5% quantiles. Compared on the same 500 replicates with `proposed` (`proposed_rect` for `R(36)`) and, in cell 8, `proposed_ew`.

## 9. Diagnostics (HEAD process, every replicate)

Per arm: the effective-sample-size ratio `(sum w)^2 / (n_a * sum w^2)`, and the largest single share of the influence sum at month 18, `max_i (w_i psi_i(18))^2 / sum_i (w_i psi_i(18))^2`.

---

## 10. Metrics

Per cell × quantity × time × estimator: replicates, evaluable replicates, failures, non-finite intervals; truth; **bias** (mean estimate − truth) with its Monte Carlo SE; empirical SD; mean SE; **SE/SD**; **coverage**; **mean interval length**. For log-rank z in null cells: rejection rate of `|z| > 1.96`.

Monte Carlo error (one SD): coverage ±0.005 at 2,000 replicates and ±0.010 at 500; SE/SD ±0.016 at 2,000. State this under the tables.

---

## 11. Steps

### Step 0 — Intake (Gate G0)
1. Branch `km-fixes-0.2.0`, clean tree. Note HEAD as `START_SHA`, `git rev-parse master` as `MASTER_SHA`.
2. `mkdir -p dev/audit/oc`; copy this file to `dev/tasks/TASK_km_oc_evaluation.md`; commit.
3. `git worktree add <temp dir outside the repo>/oc_cran 166919e`.

**G0 fails (stop) if** the branch is wrong; the tree is not clean; `dev/audit/oc/` or `dev/tasks/TASK_km_oc_evaluation.md` existed before step 2 (an earlier version of this task was started); `dev/audit/REPORT_km_fixes_stage2c.md` is not in the tree; commit `166919e` is missing or its `DESCRIPTION` version is not `0.1.0`; `pkgload::load_all()` errors for either version; or either identity check fails on one cell-5 dataset under HEAD:
- the weighted SE of `surv1` at month 18 equals `survfit(weights = w, robust = TRUE)` within 1e-6 relative;
- with `draws.band = 2000`, the variance across draws of `surv1_star` is within 0.85–1.15 of `sig2_surv1` wherever `surv1` is between 0.1 and 0.9.

### Step 1 — Scripts and sanity run (Gate G1)
Write in `dev/audit/oc/`:
- `oc_sim.R` — arguments `--version=cran|proposed --pkg=<path> --cells= --reps= --check`. `parallel::mclapply` allowed. Per-replicate results go to `dev/audit/oc/raw/`, listed in `dev/audit/oc/.gitignore`.
- `oc_boot.R` — Part E.
- `oc_summarise.R` — builds `oc_summary.csv` (long: one row per cell × quantity × time × estimator, with the §10 columns; `time` is `NA` for bands and log HR) and `oc_profiles.pdf`.

**G1 fails (stop) if** any of these holds:
1. **Truth.** On one sample of 10^6 per scenario with both potential outcomes: `|mean(T^a > t) − S_a(t)| > 0.002` for any arm and any `t` in `tt`; the cell-3 weighted mean differs from its formula by more than 0.003; or, under "conf" with true weights, the weighted mean of `1(T^a > t)` among subjects with `treat == a` differs from `S_a(t)` by more than 0.004.
2. **`cox_robust()`.** On one cell-5 dataset: `|U| / sqrt(sum(psi^2)) > 1e-6` at `bhat` for either scheme, or the FH(0,0) SE differs from `sqrt(coxph(Surv(time, status) ~ treat, weights = w, robust = TRUE, ties = "breslow")$var)` by more than 1e-6 relative.
3. **Harness**, on 200 replicates of cells 1 and 5 under both versions: checksums differ between versions; in cell 1 any `cran` or `proposed` call errors, `proposed` has `|bias|` of `delta(18)` above 0.01, or `proposed` pointwise coverage of `delta(18)` is outside 0.91–0.99.
4. The projected wall-clock of the main run exceeds 120 min.

**Part D gates** (one cell-8 dataset unless stated). Failing 5a makes all of Part D NOT DONE. Failing 5b, 5c or 5d makes only the band or RMST pieces that depend on it NOT DONE; pointwise and log HR `proposed_ew` are still delivered. Record which, and continue.
- 5a. `S_a(t_j)^2 * sum_i w_i^2 psi_i(t_j)^2` from `P_a` equals `sig2_surv0` / `sig2_surv1` of `KM_diff(draws = 0)` within 1e-8 relative wherever positive. And on 200 replicates of cell 8: the ratio of the mean `proposed_ew` SE to the mean `proposed` SE of `delta(18)` is within 0.75–0.98, and the `proposed_ew` SE/SD of `delta(18)` is within 0.80–1.15.
- 5b. The script's `− S_a(t_j) * sum_i G_ib w_i psi_i(t_j)` equals the package's `surv0_star` / `surv1_star` within 1e-8 (the script has the package's multipliers).
- 5c. With `IF` set to zero, the script's band reproduces `c_alpha_band`, `sb_lower`, `sb_upper` within 1e-8.
- 5d. With `IF` set to zero, the script's RMST pointwise limits and band limits reproduce `cumulative_rmst_bands()` within 1e-8.

### Step 2 — Main run (only if `Main run: GO`)
All 19 cells, 2,000 replicates, both versions. Then Part E if `GO`. Save the console output as `oc_stdout.txt`; record wall-clock per version.

### Step 3 — Summary, figures, report
Run `oc_summarise.R`. `git worktree remove` the 0.1.0 worktree. Write `dev/audit/REPORT_km_oc.md` (§12). Commit. Run the §13 block. Then the reporting line from the kickoff.

---

## 12. Report format — `dev/audit/REPORT_km_oc.md`

Tables first. `—` marks "not evaluated". Cells in the order of §3.

- **Top line.** Branch, SHAs, gates, parts DONE / NOT DONE, replicates run, wall-clock, failures by version.
- **Table H — headline coverage (nominal 0.95).** Rows = the 19 cells. Columns, each for `cran` | `proposed` | `proposed_ew`: `delta(18)` · band for `delta` on [3, 42] · `R(36)` · log HR FH(0,0).
- **Profile tables.** Rows = cell × estimator (`cran`, `proposed`, and where they exist `proposed_ew`, `proposed_d20`, `proposed_rect`); columns = the eleven time points, then the band on [3, 42] and the band anywhere.
  - **P1** coverage of `delta(t)`; **P2** mean interval length of `delta(t)` (band: mean width); **P3** bias of the `delta(t)` estimate, one row per cell, with the truth and the largest difference between the two versions' estimates.
  - **P4** coverage of `R(t)`; **P5** mean interval length of `R(t)`; **P6** bias of `R(t)` for `cran`, `proposed`, `proposed_rect`.
- **Table C — log HR.** Per cell: bias (Monte Carlo SE); then SE/SD, coverage and mean length for FH(0,0) [`cran`, `cran_mean1`, `cran_star`, `proposed`, `proposed_ew`] and FH(0,1) [`cran`, `proposed`, `proposed_ew`].
- **Table N — null cells (13–16).** Rejection rate of the two-sided 5% log-rank test: FH(0,0) [`cran`, `cran_mean1`, `proposed`, `proposed_ew`]; FH(0,1) [`cran`, `proposed`, `proposed_ew`].
- **Table B — bootstrap benchmark.** Cells 5, 8, 10 × `delta(18)`, `R(36)`, log HR: coverage and mean length for `proposed`, `proposed_ew` (cell 8), `boot`, `boot_pct`, same 500 replicates.
- **Table D — diagnostics and checks.** Per cell: effective-sample-size ratio and largest influence share, by arm; evaluable replicates at months 42 and 48; non-finite counts by version; mean band critical value by version; the G1 check values; package Cox output, HEAD against `cran`.
- **Appendix.** Coverage, length and SE/SD of `S_0(t)` and `S_1(t)` in the profile layout; SE/SD of `delta(t)`; the exact calls and returned fields used.
- **Figures — `oc_profiles.pdf`,** base graphics.
  - Page 1: coverage of `delta(t)` against t, one panel per cell, one line per estimator, a line at 0.95 and dashed lines at 0.95 ± 2 Monte Carlo SD.
  - Page 2: the same for `R(t)`.
  - Page 3: mean interval length of `delta(t)`, each estimator relative to `proposed`, against t.
  - Page 4: coverage by cell for the band for `delta`, the RMST band and log HR FH(0,0) — open circles `cran`, filled `proposed`, squares `proposed_ew`, triangles `boot`.
- **Reading.** At most twelve lines of plain language: where the two versions agree; where they separate, by how much, in which direction, and at which times; what `proposed` reaches in the F cells, at N = 200 and at months 42 and 48; what carrying the weight estimation through changes, in coverage and in length; what the 20-draw path costs; what the rectangle rule changes for `R(t)`; how these compare with the bootstrap.

Findings only. No recommendation and no follow-up task proposals.

---

## 13. Gate G2 — post-conditions

Run from the repo root with `START_SHA` and `MASTER_SHA` set, after the report is committed. Every line must print `OK`.

```bash
O=dev/audit/oc
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"
git diff --name-only "$START_SHA" HEAD | grep -Ev '^dev/(tasks|audit)/' && echo "FAIL: paths outside dev/tasks, dev/audit" || echo "OK 2"
for f in oc_sim.R oc_boot.R oc_summarise.R oc_summary.csv oc_profiles.pdf oc_stdout.txt .gitignore; do
  [ -f "$O/$f" ] || echo "FAIL: missing $O/$f"; done
[ -f dev/audit/REPORT_km_oc.md ] && [ -f dev/tasks/TASK_km_oc_evaluation.md ] && echo "OK 3 (if no FAIL above)" || echo "FAIL: report or task doc missing"
Rscript -e 'd <- read.csv("dev/audit/oc/oc_summary.csv")
  m <- d[d$quantity == "delta" & !is.na(d$time) & d$estimator %in% c("cran", "proposed"), ]
  stopifnot(length(unique(d$cell_id)) == 19, nrow(m) == 19 * 11 * 2, all(m$n_rep == 2000),
            all(is.finite(m$coverage[m$estimator == "proposed" & m$time <= 42])),
            all(is.finite(m$bias[m$estimator == "proposed" & m$time <= 42])))
  cat("OK 4\n")'
Rscript $O/oc_sim.R --version=proposed --pkg=. --check > /dev/null 2>&1 && echo "OK 5" || echo "FAIL: not reproducible"
git worktree list | grep -q oc_cran && echo "FAIL: 0.1.0 worktree not removed" || echo "OK 6"
[ -z "$(git status --porcelain)" ] && echo "OK 7" || echo "FAIL: uncommitted changes"
```

Notes:
- Check 4 assumes the summary names the columns `cell_id`, `quantity`, `time`, `estimator`, `n_rep`, `coverage`, `bias`, with `quantity == "delta"` for the survival difference.
- Check 5: `--check` re-runs replicates 1–10 of cells 1 and 5 and compares with the stored per-replicate results within 1e-10; it exits 0 only then.
- If Part D or Part E is `HOLD` or NOT DONE, its columns, tables and symbols are omitted and the report says so; no check depends on them.
