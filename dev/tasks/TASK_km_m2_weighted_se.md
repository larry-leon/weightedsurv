# TASK — weightedsurv: M2, weighted-SE evaluation when weights are related to the hazard

Written in chat 2026-11-19. Executed by Claude Code in the local `weightedsurv` repo, on branch `km-fixes-0.2.0`.

| | |
|---|---|
| Touches `R/` | **No.** Evaluation only. |
| Writes | Only under `dev/tasks/` and `dev/audit/` (new folder `dev/audit/m2/`). |
| Branch | `km-fixes-0.2.0`, so that "package" means the code Larry approved in Stage 2b. No merge. **Never push.** |
| Compute | **M2: 7 cells (10 arm-cells) × 2,000 replicates**; n = 400 for one-group cells, N = 800 for two-arm cells. Expected 10–15 min. Run only if the kickoff says `M2: GO`. Stop and record if the projection passes 30 min. |
| Gates | Stop-on-failure, not stop-to-ask. |

---

## 1. Why

Stage 2b made the package's analytic weighted Kaplan–Meier SE the weighted Greenwood formula. M1, M1b and M1c showed it calibrated — but in every one of those cells the weights were drawn independently of the outcome.

A large-sample calculation and a short simulation in chat indicate that when the weight is related to the hazard, weighted Greenwood is not consistent for the true variance: too large when high weights go with high hazard, too small when they go with low hazard. Propensity weights are of that kind. A robust (infinitesimal-jackknife) variance does not need the assumption.

This task measures that, so Larry can decide what the analytic weighted SE should be. It decides nothing itself.

---

## 2. Design

Censoring `C ~ U(0, 40)`, independent, in every cell. Evaluation times `t = 5, 10, 20`. Replicate `r` of cell number `k` uses `set.seed(100000 * k + r)`, so results do not depend on run order or parallelism.

### Cells

| k | Cell | Data | Weights |
|---|---|---|---|
| 1 | **A** — weight unrelated to hazard | one group, n = 400; `X ~ Bernoulli(0.5)`; `T ~ Exp(0.10)` | `w = 1 + 2X` |
| 2 | **B** — high weight, high hazard | as A, but hazard 0.05 if X = 0, 0.15 if X = 1 | `w = 1 + 2X` |
| 3 | **C** — high weight, low hazard | as A, but hazard 0.15 if X = 0, 0.05 if X = 1 | `w = 1 + 2X` |
| 4 | **Bt** — B with tied times | cell B, observed times rounded up to the next 0.5 | `w = 1 + 2X` |
| 5 | **D** — true propensity weights | two arms, N = 800; `X ~ N(0, 1)`; `A ~ Bernoulli(expit(-0.5 + X))`; `T ~ Exp(0.08 * exp(0.7 X - 0.4 A))` | treated `1 / e(X)`, control `1 / (1 - e(X))`, with the true `e(X)` |
| 6 | **E** — estimated propensity weights | as D | the same, with `e(X)` from `glm(A ~ X, binomial)` fitted in each replicate |
| 7 | **F** — extreme weights | as D, but `A ~ Bernoulli(expit(-0.5 + 2 X))` | true weights |

Cells D, E and F are evaluated separately in each arm, giving 10 arm-cells: A, B, C, Bt, D0, D1, E0, E1, F0, F1.

### Truth

The target in each arm-cell is the limit of the weighted KM, `E[w 1(T > t)] / E[w]`.
- A, B, C, Bt: closed form from the two-point mixture. (For Bt the three times are multiples of 0.5, so rounding does not change the target.)
- D, E, F: the counterfactual marginal survival `S_a(t) = E_X[ exp(-0.08 * exp(0.7 X - 0.4 a) * t) ]`, by numerical integration over `X ~ N(0, 1)`. E has the same target as D.

### Point estimate

The weighted KM of the group, from the package: `calculate_risk_event_counts(U, D, w, at)` then `KM_estimates()`, with `at` the sorted unique observed times of the group, read at the last grid point `<= t`. Confirm once per cell that it equals `survfit(weights = w)` within 1e-10.

### SE estimators

Notation for one group: `Y_w(s) = sum_i w_i 1(U_i >= s)`, `dN_w(s) = sum_i w_i D_i 1(U_i = s)`, `S` the weighted KM.

| Name | Definition |
|---|---|
| `pkg` | The package's analytic SE: `sqrt(sig2_KM)` from the two calls above with `draws = 0`. |
| `pkg_resampling` | The same with `draws = 200`, `seedstart = 816951 + r`. **First 500 replicates of each arm-cell only.** |
| `robust` | `summary(survfit(Surv(U, D) ~ 1, weights = w, robust = TRUE), times = t, extend = TRUE)$std.err` |
| `ij` | Explicit infinitesimal-jackknife formula, coded in the script: `psi_i(t) = D_i 1(U_i <= t) / (Y_w(U_i) - dN_w(U_i))  -  sum over s <= min(t, U_i) of dN_w(s) / (Y_w(s) * (Y_w(s) - dN_w(s)))`, then `SE(t) = S(t) * sqrt( sum_i w_i^2 psi_i(t)^2 )`. Terms with `Y_w - dN_w <= 1e-10 * Y_w` are 0. |
| `xie_liu` | `SE(t) = S(t) * sqrt( sum over s <= t of (1 - q(s)) / (M(s) * q(s)) )`, with `q(s) = 1 - dN_w(s) / Y_w(s)` and `M(s) = Y_w(s)^2 / sum_i w_i^2 1(U_i >= s)`. Terms with `q(s) <= 1e-10` are 0. |

Two cross-checks on the hand-coded formulas, on one dataset from cell D before the main run:
- `ij` against `robust`: report the maximum relative difference at the three times. Chat found 1e-15 with and without ties; a difference above 1e-6 means the formula was mis-coded — fix the script, not the definition.
- `xie_liu` against `adjustedCurves::adjustedsurv(method = "iptw_km", conf_int = TRUE)` on the same two-arm dataset and weights, if `adjustedCurves` is installed: report the maximum relative difference in `se`. If they differ by more than 1e-6, keep both columns in the output, named `xie_liu` and `adjustedCurves`, and say so in the report. If the package is not installed, say so and continue.

### Recorded per arm-cell × estimator × time

Mean SE; empirical SD of the point estimate over all 2,000 replicates; **SE/SD**; **coverage** of `S ± 1.96 * SE` for the truth; number of non-finite SEs. Per arm-cell × time also: bias of the point estimate, and the mean over replicates of the weights' effective-sample-size ratio `(sum w)^2 / (n * sum w^2)`.

Monte Carlo error at 2,000 replicates is about ±0.016 on SE/SD and ±0.005 on coverage (one SD). State this under the tables.

---

## 3. Steps

### Step 0 — Intake (Gate G0)
1. Confirm the branch is `km-fixes-0.2.0` and the tree is clean. Note HEAD as `START_SHA`.
2. `mkdir -p dev/audit/m2`; copy this file to `dev/tasks/TASK_km_m2_weighted_se.md`; commit.

**G0 fails (stop) if** the branch is wrong, the tree is not clean, or `pkgload::load_all(".")` errors.

### Step 1 — Script and sanity run (Gate G1)
Write `dev/audit/m2/km_m2_montecarlo.R`. It loads the package with `pkgload::load_all(".")`. `parallel::mclapply` is allowed. Give it a `--check` mode that re-runs the first 20 replicates of cell B serially and compares them with the stored results.

Run the two formula cross-checks, then 200 replicates of cell A.

**G1 fails (stop) if** any of these holds on the 200 replicates of cell A: absolute bias of the point estimate above 0.01 at any time; the point estimate differs from `survfit` by more than 1e-10; any of `pkg`, `robust`, `ij`, `xie_liu` has SE/SD outside 0.80–1.20; or the projected wall-clock for the full run exceeds 30 min.

### Step 2 — Full run
All 10 arm-cells. Write `dev/audit/m2/km_m2_results.csv` (one row per arm-cell × estimator × time: at least 150 rows) and `km_m2_stdout.txt`. Record the wall-clock.

### Step 3 — Report
`dev/audit/REPORT_km_m2.md`, format in §4. Commit. Run the §5 block. Then run the reporting line from the kickoff.

---

## 4. Report format — `dev/audit/REPORT_km_m2.md`

- **Top line.** Branch, SHA, gates, wall-clock, replicates actually run.
- **Table 1 — t = 10.** Rows = the 10 arm-cells; columns = `pkg`, `pkg_resampling`, `robust`, `ij`, `xie_liu`. Each entry: SE/SD (coverage).
- **Table 2 — t = 5 and t = 20.** Same layout, SE/SD only, the two times side by side.
- **Table 3 — checks.** Per arm-cell: truth and bias of the point estimate at the three times; effective-sample-size ratio; `ij` versus `robust` maximum relative difference; the `adjustedCurves` cross-check result.
- **Reading.** At most eight lines of plain language: where the estimators agree, where they separate, by how much and in which direction, and whether estimating the weights (E against D) or extreme weights (F against D) changes the picture.

Findings only. No recommendation and no follow-up task proposals: the choice of estimator is Larry's.

---

## 5. Gate G2 — post-conditions

Run from the repo root with `START_SHA` set, after the report is committed. Every line must print `OK`.

```bash
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1" || echo "FAIL: wrong branch"
git diff --name-only "$START_SHA" HEAD | grep -Ev '^dev/(tasks|audit)/' && echo "FAIL: paths outside dev/tasks, dev/audit" || echo "OK 2"
for f in audit/m2/km_m2_montecarlo.R audit/m2/km_m2_results.csv audit/m2/km_m2_stdout.txt audit/REPORT_km_m2.md tasks/TASK_km_m2_weighted_se.md; do
  [ -f "dev/$f" ] || echo "FAIL: missing dev/$f"; done; echo "OK 3 (if no FAIL above)"
Rscript -e 'd <- read.csv("dev/audit/m2/km_m2_results.csv"); stopifnot(nrow(d) >= 150); cat("OK 4\n")'
Rscript dev/audit/m2/km_m2_montecarlo.R --check > /dev/null 2>&1 && echo "OK 5" || echo "FAIL: not reproducible"
[ -z "$(git status --porcelain)" ] && echo "OK 6" || echo "FAIL: uncommitted changes"
```
