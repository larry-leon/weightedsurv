# weightedsurv 0.2.0: estimated-weights statistical appendix on branch `km-fixes-0.2.0`

Task: `dev/tasks/TASK_estimated_weights_appendix.md`. Start: `85c0abc` (contains `a15ff33`). Nothing pushed, nothing merged.

Deliverable: `vignettes/articles/estimated_weights_methodology.Rmd`, "Estimated propensity-score weights: statistical details and properties" (9 numbered sections, 34 displayed equations, 8 evaluated chunks). Rendered HTML: `dev/audit/estimated_weights_appendix/estimated_weights_methodology.html`.

Commits after `85c0abc`:

| commit | content |
|---|---|
| `c49d238` | task document copied to `dev/tasks/` |
| `69339cd` | the new article; three entries added to `vignettes/weightedsurv_references.bib` |
| `17c3469`, `34eb402`, `f91e3c5`, `f3a52ce` | one pointer sentence in `weightedsurv_examples.Rmd`, `weightedkm_methodology.Rmd`, `weightedcox_methodology.Rmd`, `weightedrmst_methodology.Rmd` |
| `6da5721` | `_pkgdown.yml`: one line in the "Methodology" articles index |
| `0c66cbe` | `dev/scratch/render_vignettes.R`: entry for the new article |
| `e95728a` | the article: one sentence added to section 1 (derivation steps labelled as background) |
| final commit | this report, gate scripts and outputs, rendered HTML |

Files outside the task's list of sources that were changed: `vignettes/weightedsurv_references.bib` (under `vignettes/`).

## 1. Gate table

Outputs: `dev/audit/estimated_weights_appendix/gates_output.txt`, `g9_numerals.txt`, `00check_0c66cbe.log`, `00check_a15ff33_same_session.log`, `testthat_0c66cbe.Rout`. Scripts: `dev/scratch/render_vignettes.R`, `gates_appendix.R`, `g9_appendix_numerals.py`. The package used for rendering is the build of `ac76f82` installed in the previous task; `R/` has not changed since (`git diff --name-only ac76f82..HEAD -- R/ man/ tests/ DESCRIPTION NAMESPACE` is empty).

| # | Result | Evidence |
|---|---|---|
| G1 | PASS | `git diff --name-only a15ff33..HEAD` lists files under `vignettes/`, `dev/` and `_pkgdown.yml` only. `git diff a15ff33..HEAD -- _pkgdown.yml`: one added line, `      - articles/estimated_weights_methodology`. |
| G2 | PASS | `rmarkdown::render()` with `options(warn = 2)`: `OK` for `estimated_weights_methodology.Rmd` (1 s), `weightedsurv_examples.Rmd` (8 s), `weightedkm_methodology.Rmd`, `weightedcox_methodology.Rmd` (2 s), `weightedrmst_methodology.Rmd`. Every citation key of the new article resolves (8 entries in its reference list). |
| G3 | PASS | The chunks of the article run in 1.0 s (limit 60). Its 115 code lines contain no `replicate(`, `boot`, `for (b in`, `mclapply`, `future` and no `glm(`; one call to `ps_weights()`. A `grep` of the whole file for `boot` matches 5 prose lines (the word "bootstrap" in the text of section 7 and the summary, and the "boot" rows of the copied table), none in code. |
| G4 | PASS | Largest absolute difference between a by-hand value and the package value: 3.55e-15 (limit 1e-8). Table in section 3. Nothing was omitted under the exported-only rule. |
| G5 | PASS | `grep -l "estimated_weights_methodology"` over `vignettes/*.Rmd` and `vignettes/articles/*.Rmd` returns the four files and no other (`grep -rl` over `vignettes/` gives the same four). Each of the four diffs is 1 line changed, with one added sentence and nothing removed (`git diff --word-diff`: one added run each). |
| G6 | PASS | `pkgdown::check_pkgdown()`: "No problems found." (pkgdown 2.2.1). The site was not built. |
| G7 | PASS | `R CMD check --as-cran --timings` on the tarball built from a clean `git clone` at `0c66cbe` (`RSTUDIO_PANDOC` set to the quarto tools directory, `DISPLAY` unset): `Status: 2 NOTEs`, 0 errors, 0 warnings. NOTEs: top-level files; HTML version of the manual. The identical procedure on `a15ff33` in the same session: `Status: 2 NOTEs`, the same two. Tests in both: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`. The commits after `0c66cbe` change `vignettes/articles/estimated_weights_methodology.Rmd` and files under `dev/` only; both paths are in `.Rbuildignore`, so the tarball is the one checked. |
| G8 | PASS | Section 2: 34 displayed equations, 34 rows. 20 rows point to code. 14 rows are derivation steps or background, each labelled in the article (section 1 states that the lemmas, the expansions and the projection argument are background; individual formulas carry "*background*" where they stand alone). No row is "not found". |
| G9 | PASS | `g9_numerals.txt`: 235 numerals in the prose of the new article outside the copied table (78 in plain text, 5 in code spans, 152 in math), 0 without a class; 1 numeral on the changed line of `weightedsurv_examples.Rmd` (`FH(0,1)`, text of `a15ff33`). Classes in plain text: section, lemma, stage and table references; list-item numbers; scheme names (`FH(0, 0)`); the values 0 and 1 in conditions; the design of the copied table (nominal 0.95, two-sided 5%, month 18, 2,000 and 500 replicates); the YAML header. In math: indices, exponents, the tolerance $10^{-10}$, the constants 1.96 and 0.95 of `cumulative_rmst_bands()`, and the Monte Carlo errors 0.005 and 0.010 quoted from the source of the table. The copied evaluation table has 20 lines with 75 numerals. Every result of section 8 is inline R or chunk output. |

## 2. Formula–code map

E-numbers follow the order of the displayed equations in the article (line of the `.Rmd` in brackets). "Statement" = a step of a derivation, labelled in the article; the quantities in it are computed where the other rows say.

| # | Formula (section) | File, function, object |
|---|---|---|
| E1 [45] | weights, `"ipw"` and `"stabilized"`; $\hat\pi_a$ (2) | `R/ps_weights.R`, `ps_weights()`: `w`, `p1 <- mean(A)`, `sc`, `w <- sc * w`; returned as `$weights` |
| E2 [52] | $\bar Y_a^w$, $d\bar N_a^w$, $\hat S^a$ (2) | `R/df_counting_helpers.R`, `calculate_risk_event_counts()`: `ybar`, `nbar` (cumulative; $d\bar N$ is `diff(c(0, nbar))`); `R/km_wlr_calculations_helpers.R`, `KM_estimates()`: `S_KM`; per arm through `km_arm_on_grid()` (`R/internal_helpers.R`) in `KM_diff()` |
| E3 [63] | $e_i(\gamma)$, $\hat e_i$ (3) | `ps_weights()`: `stats::glm(formula, family = stats::binomial(), data = data)`, `e <- fitted(model)`; returned as `$ps` |
| E4 [67] | score $s_i$, information $\mathcal I_n$ (3) | `ps_weights()`, inside the line for `IF`: `xm * (A - e)` and `crossprod(xm * sqrt(e * (1 - e)))` |
| E5 [73] | expansion of $\hat\gamma$ (statement); $IF_i$ (3) | `ps_weights()`: `IF <- (xm * (A - e)) %*% solve(crossprod(xm * sqrt(e * (1 - e))))`; returned as `$IF` |
| E6 [79] | $dw_i$ (3) | `ps_weights()`: `dw <- ifelse(A == 1, -(1 - e) / e, e / (1 - e)) * xm`, `dw <- sc * dw` for stabilized; returned as `$dw` |
| E7 [91] | $\psi_i(t)$ (4.1) | `R/internal_helpers.R`, `ws_km_psi_parts()`: `a`, `H`, `exit` ($\psi_i(g_j) = D_i 1(k_i \le j) a[k_i] - H[\min(j, k_i)]$); tolerance in `ok` |
| E8 [95] | $\psi_i(t)$ as a sum of residual increments (4.1) | the same quantity as E7 (identity proved in the article); not coded in this form |
| E9 [101] | Lemma 1 | statement. Checked numerically in `gates_appendix.R`: central differences of `KM_diff()$surv1`, `$surv0` in one subject's weight, 10 subjects, largest relative difference 4.4e-07 |
| E10 [107] | $\sum_{i \in a} w_i \psi_i(t) = 0$ (4.1) | statement. Shown in section 8 of the article (7.9e-17, 1.3e-16) |
| E11 [113] | representation under known weights; $\psi_i^0$ (4.1) | background (labelled); cited |
| E12 [119] | $\partial \hat S^a / \partial\gamma = -\hat S^a D_a$; $D_a(t)$ (4.1) | `R/estimated_weights.R`, `ws_km_ew_group()`: `Dg <- ws_km_ij_draws(parts, rep(1, length(U)), ew$dw[rows, ])`. The derivative identity itself is a statement |
| E13 [123] | expansion (statement); $\Phi_i^a(t)$ (4.1) | `ws_km_ew_group()`: not formed subject by subject; its sum of squares is `sumsq`, its resampled sums are in `KM_diff()` (`arm_star`). Subject by subject in `wlr_dhat_estimates()` (`d`) |
| E14 [129] | expansion (statement); $\xi_i(t)$ (4.2) | `R/km_wlr_calculations_helpers.R`, `wlr_dhat_estimates()`: `d` after `d <- d + as.numeric(ew$IF %*% dDg)` equals $-\xi_i$ at the grid point of `tzero` |
| E15 [133] | cross term; $E_a(t)$, $M$ (4.2) | `KM_diff()`: `cross <- rowSums(e1$A * Dg0) + rowSums(e0$A * Dg1) + rowSums((Dg1 %*% ew$CIF) * Dg0)`; $E_a$ is `A` of `ws_km_ew_group()`; $M$ is `CIF = crossprod(IF)` of `ws_ew_setup()` |
| E16 [139] | $\hat R(t_j)$, trapezoid and rectangle (4.3) | `R/internal_helpers.R`, `ws_cumulative_mean()`; called in `R/rmst_calculations.R`, `cumulative_rmst_bands()`: `rmst_time` |
| E17 [147] | $S_0$, $\bar A$, $d\hat\Lambda$ (4.4) | `R/internal_helpers.R`, `ws_score_psi()`: `S0`, `zb` (`S1 / S0`), `dL` |
| E18 [151] | $\psi_i(\beta)$ (4.4) | `ws_score_psi()`: `psi[o] <- w * (D * Kd[grp] * (Z - zb[grp]) - eb * (Z * cK[grp] - cKz[grp]))` |
| E19 [155] | $U_w(\beta)$, $I_w(\beta)$ (4.4) | `R/weightedCox.R`, `score_calculation()`: `score`, `i_bhat`, from the weighted counts built in `cox_score_rhogamma()` |
| E20 [162] | Lemma 2 | statement. Checked numerically in `gates_appendix.R`: central differences of `cox_score_rhogamma()` in one subject's weight, 10 subjects, largest relative difference 4.1e-07 |
| E21 [168] | $B(\beta)$, $\psi_i^{ew}(\beta)$ (4.4) | `R/estimated_weights.R`, `ws_score_ew()`: `psi + IF %*% colSums((psi / a) * dw)` |
| E22 [174] | z, standard error, interval (4.4) | `R/weightedCox.R`, `cox_rhogamma()`: `z.score <- u.zero / sqrt(sig2_score)` with `sig2_score <- sum(psi_zero^2)`; `sig_bhat_asy <- sqrt(sum(psi_bhat^2) / i_bhat^2)`; `ci_cox()`. Log-rank z also in `wlr_cumulative()` and `df_counting()` (`z.score`) |
| table (4.4) | $\kappa(s)$ by scheme | `R/weighted_helpers.R`, `wt.rg.S()`: `wt`; through `get_validated_weights()` |
| E23 [197] | $r_i^{ew}$ (4.4) | `R/df_counting.R`, `df_counting()`: `ew_adj()` |
| E24 [201] | $V^{ew}(t)$, $B(t)$ (4.4) | `R/estimated_weights.R`, `ws_lr_cumvar_ew()`; called in `wlr_cumulative()`. Variance and covariance of the difference: `wlr_dhat_estimates()`: `sig2_dhat <- sum(d^2)`, `cov_wlr_dhat <- sum(-psi * d)` |
| E25 [209] | the three variances (5) | `ws_km_ew_sig2()` (`S^2 * sumsq`) and `KM_diff()` (`sig2_surv1`, `sig2_surv0`); `KM_diff()`: `sig2_dhat`; `cox_rhogamma()`: `sig2_score`, and `sum(ws_score_ew(...)^2)` in `wlr_cumulative()` and `df_counting()` |
| E26 [217] | decomposition; $C$, $M$ (5) | `ws_km_ew_group()`: `sumsq <- V + 2 * rowSums(A * Dg) + rowSums((Dg %*% ew$CIF) * Dg)`; `ws_lr_cumvar_ew()`: `V + 2 * rowSums(B * C) + rowSums((B %*% crossprod(IF)) * B)`. For the score at one $\beta$ the package computes the left-hand side directly |
| E27 [223] | $E\{\partial m / \partial\gamma\} = -E\{m s\}$ (5) | statement (heuristic argument); not computed |
| E28 [227] | residual of the projection (5) | statement; not computed |
| E29 [231] | variance of the residual (5) | statement; not computed. Its sample counterpart $\sum c_i^2 - C^\top M^{-1} C$ is labelled background and is shown in section 8 from article code |
| E30 [243] | resampled curve and difference (6) | `KM_diff()`: `IFG`, `arm_star()`, `surv0_star`, `surv1_star`, `dhat_star` |
| E31 [250] | band critical value and band (6) | `KM_diff()`: `Zdhat_star`, `sups`, `c_alpha_band`, `sb_lower`, `sb_upper` |
| E32 [256] | RMST limits and band (6) | `cumulative_rmst_bands()`: `cum_rmst_stars`, `sig2_rmst_time`, `rmst_time_lower/upper` (1.96), second `KM_diff()` call (`seedstart = 99999`), `c_alpha_band` (0.95 quantile), `rmst_time_sb_lower/upper` |
| E33 [262] | resampled scores, $\hat\beta^\dagger$ (6) | `R/weightedCox.R`, `cox_rhogamma_resample()`: `score_star`, `score_star_null`, `bhat_center_star`, `sig_bhat_star`; `cox_rhogamma()`: `bhat_debiased`, `hr_ci_star` |
| E34 [287] | $\text{ESS}_a$ (7) | `df_counting()`: `ess_arm()`, `ans$ess`; also `print.ps_weights()` |

Symbols defined in the text without a display: $Y_i$, $N_i$ (data; `outer(U, at_points, ">=")` in `calculate_risk_event_counts()`); $G_i$ (`KM_diff()`: `G0.draws`, `G1.draws`; `cox_rhogamma_resample()`: `G`); $\hat\beta_w$ (`find_cox_root()`); $\mathcal L_j$ (`ws_cumulative_mean()`); $c_i(t)$ (`ws_lr_cumvar()` / `ws_lr_cumvar_ew()`, not stored per subject); $\gamma^*$, $\beta_w^*$, $S^a$, $\Delta$, $\Lambda^a$, $y_a$, $m_i$, $P_\gamma$, $\Pi$ (population quantities, not computed).

Notation. The article follows `weightedkm_methodology.Rmd`: $T_i$, $\Delta_i$, $A_i$, $w_i$, $\bar Y_a^w$, $d\bar N_a^w$, $\hat S^a$, $\psi_i(t)$, $\Phi_i^a(t)$, $D_a(t)$, $\hat\Delta(t)$, $G_i$. Where the three articles differ: `weightedcox_methodology.Rmd` writes $X_i$ for the time, $Z_i$ for the arm, $a_i$ for the subject weight and $w(t)$ for the time weight; the new article writes $T_i$, $A_i$, $w_i$ and, to keep $w$ for the subject weight, $\kappa(s)$ for the time weight (stated in its section 4.4). Symbols new to the vignettes: $\gamma$ for the coefficients of the propensity model (`weightedkm_methodology.Rmd` writes $\hat\beta^T X_i$ in its recipe; $\beta$ is the log hazard ratio here), $\xi_i(t)$, $E_a(t)$, $M$, $\mathcal L_j$, $\psi_i^0$.

## 3. Section 8 quantities

Restricted Rotterdam population (829 women, 291 treated; ESS ratio 0.930 treated, 0.984 control), stabilized weights, `weight_variance = "estimated"`. From `gates_output.txt`.

| quantity | by hand | package | difference |
|---|---|---|---|
| SE of the survival difference at 60 months | 0.035166439891200 | 0.035166439891200 (`sqrt(KM_diff(at_points = 60)$sig2_dhat)`) | -1.39e-17 |
| SE of the log hazard ratio, FH(0,0) | 0.084025519440320 | 0.084025519440320 (`cox_rhogamma()$fit$sig_bhat_asy`) | 1.39e-17 |
| log-rank z, FH(0,0) (shown in addition) | 4.868643752370565 | 4.868643752370569 (`cox_rhogamma()$z.score`) | -3.55e-15 |

- How the by-hand values are built: `time`, `delta`, `z`, `w_hat` and `ew_parts` (`IF`, `dw`) of the `df_counting()` result; `ew_parts` equals `pw$IF` and `pw$dw` with the rows in the order of `time` (largest difference 0). The per-subject contributions $\psi_i(t)$ and $\psi_i(\beta)$ are computed in the article from the displayed formulas with base R; no `:::`.
- The score contributions computed from the formula equal those the package uses, obtained from the exported `cox_rhogamma_resample()` with unit multipliers: largest difference 3.1e-15 (shown in the article).
- Also shown in the article: the three terms of the variance and the sample projected variance.

| variance of | known | cross term | quadratic term | estimated (sum) | projected (background) |
|---|---|---|---|---|---|
| survival difference at 60 months | 0.00139487 | -0.000325549 | 0.000167357 | 0.00123668 | 0.00123547 |
| score at the estimate, FH(0,0) | 118.593 | -36.897 | 18.797 | 100.493 | 100.374 |

## 4. How each statement is labelled

**Proved in the article.**
- The derivative of the weights, $dw_i$ (section 3).
- The two expressions for $\psi_i(t)$ are equal (4.1).
- Lemma 1: $\partial \hat S^a(t) / \partial w_i = -\hat S^a(t) \psi_i(t)$.
- $\sum_{i \in a} w_i \psi_i(t) = 0$, and the Kaplan–Meier quantities do not depend on the arm proportion of stabilized weights (4.1).
- $\partial \hat S^a / \partial \gamma = -\hat S^a D_a$ (chain rule, 4.1).
- The cross term of the two arms (4.2).
- The compensator terms of the score sum to zero, so the two expressions for $U_w$ are equal (4.4).
- Lemma 2: $\partial U_w(\beta) / \partial w_j = -\psi_j(\beta) / w_j$ with the time weights held fixed.
- The decomposition of the variance into the known part, a cross term and a quadratic term (5).
- Conditionally on the data the resampled difference is Gaussian with mean zero, and its variance equals the analytic variance (6).

**Cited.**
- Linear representation of $\hat\gamma$ (M-estimation): Lunceford and Davidian (2004), section 3.1.
- Martingale representation of the weighted estimator with a known propensity score: Deng and Wang (2025); Fleming and Harrington (1991) for the unweighted case.
- Robust score test and standard error of `coxph()` under known weights: Lin and Wei (1989).
- Variance with an estimated propensity score equal to the variance with the known one minus a quadratic form, for inverse-probability-weighted means: Lunceford and Davidian (2004), section 3.1, equations (13)–(16) (read in the paper).
- Efficiency of weighting by an estimated propensity score: Hirano, Imbens and Ridder (2003) (statement taken from the abstract).
- Multiplier resampling for the unweighted estimators: Lin, Wei and Ying (1993); Dobler, Beyersmann and Pauly (2017).
- Origin of the IPTW Kaplan–Meier estimator and log-rank test: Xie and Liu (2005) (statement taken from the abstract).

**Heuristic arguments (labelled as such in the article; not proved).**
1. The representation of $\hat\gamma$ still holds, around the limit of the fit, when the propensity model is wrong (3).
2. The first-order expansion of $\hat S^a(t; \hat\gamma)$ giving $\Phi_i^a$; the remainder is not bounded (4.1). The same for the survival difference and the RMST (4.2, 4.3).
3. The expansion of the score in $\gamma$ giving $\psi_i^{ew}$, and the expansion of $\hat\beta_w$ giving its standard error (4.4).
4. Each sum of squared corrected contributions estimates the variance of its estimator (5).
5. The projection argument: with a correct model the large-sample variance under estimated weights is no larger than the `"known"` formula estimates (5). Steps not proved: interchange of differentiation and expectation for the survival contributions; convergence of the derivative and of the estimated contributions.
6. Under a wrong model the linearization still describes the variability around the limit, and the ordering of the two variances is not implied (5).
7. Validity of the multiplier resampling with corrected contributions (conditional multiplier central limit theorem with estimated contributions), under the four conditions listed in section 6.

**Stated as not derived.** The effect of ignoring the dependence of the time weights on the estimated subject weights (4.4), and of ignoring the estimation of the arm proportion in the log-rank and Cox score (3).

## 5. `dev/tasks/TASK_km_fixes_stage2e.md` §2 against the code

No disagreement in the method was found: `IF`, `dw`, both weight types, $D_g$, $\Phi_i$, the variance of a group's curve and of the difference, the resampled curves, $B(\beta)$, $\psi^{ew}$, the single sum over strata, the cumulative variance, the covariance with the survival difference and the `coxph()`-residual route are coded as §2 states them.

Details present in the code and not in §2:
- `ws_km_ew_group()` returns `pmax(sumsq, 0)` and `KM_diff()` returns `pmax(sig2_surv1 + sig2_surv0 - 2 * surv1 * surv0 * cross, 0)`: a variance that rounding makes negative is set to 0.
- The variance of a curve is set to 0 where the curve is 0 (`ws_km_ew_sig2()`, `KM_diff()`).
- $D_g$ is computed as `ws_km_ij_draws(parts, rep(1, n_g), dw_g)`; the implementation note of §2 (marked there "not a contract") writes `ws_km_ij_draws(parts, W, dw_g / W)`. The two are the same quantity.

On `draws > 0`: §2 does not mention it. The task document lists "`draws` is 0" among the validation rules (W2) and gives no reason; `dev/audit/REPORT_km_fixes_stage2e.md` gives none either. The article states what the source does (the `draws > 0` route in `calculate_risk_event_counts()` draws multipliers for one group and uses that group's contributions only; `ws_ew_setup()` stops unless `draws` is 0) and that the record gives no further reason.

## 6. References added

Added to `vignettes/weightedsurv_references.bib`:

1. Lunceford, J. K. and Davidian, M. (2004). Stratification and weighting via the propensity score in estimation of causal treatment effects: a comparative study. *Statistics in Medicine* 23(19), 2937–2960. doi:10.1002/sim.1903. Key `LuncefordDavidian2004`. Details confirmed from the title page of the paper; section 3.1 read.
2. Hirano, K., Imbens, G. W. and Ridder, G. (2003). Efficient estimation of average treatment effects using the estimated propensity score. *Econometrica* 71(4), 1161–1189. doi:10.1111/1468-0262.00442. Key `HiranoImbensRidder2003`. Details confirmed from the Econometric Society and EconPapers records; the paper itself was not read.
3. Xie, J. and Liu, C. (2005). Adjusted Kaplan–Meier estimator and log-rank test with inverse probability of treatment weighting for survival data. *Statistics in Medicine* 24(20), 3089–3110. doi:10.1002/sim.2174. Key `XieLiu2005`. Details confirmed from the PubMed record; the paper itself was not read.

Reused from the file: `LW89` (Lin and Wei 1989, JASA 84, 1074–1078), `LinWeiYing1993` (Lin, Wei and Ying 1993, Biometrika 80(3), 557–572), `DengWang2025`, `FH1991`, `DBP2017`.

Not cited: Leon, Cai and Wei (2009), named in the help page of `df_counting()`; a search confirmed a 2005 working-paper version only, not the journal details.

## 7. Omitted from section 3 of the task, and why

- **Section 5, "variance as the sum of squared corrected contributions, for each estimator".** Given for the curve, the difference and the score. For the RMST the package has no such sum: its variance is the sample variance over resampling draws. The article says so (4.3, 5, 6).
- **Section 6, why `draws > 0` is not available.** The record gives no reason (section 5 above); the article reports the source and says the record gives none.
- **Section 7, the evaluation table.** It is a selection of rows (cells 31, 32, 35, 37, 38, 39) and columns (survival difference at month 18, log hazard ratio FH(0,0), null rejection FH(0,0)) of "Table M7", with six rows ("known", "estimated" on the first 500 replicates of cells 38 and 39) taken from the table "Bootstrap benchmark on the same first 500 replicates" of the same report section, so that the bootstrap row stands beside values from the same replicates. The "Regime" column is a description of each cell written for the article, not a column of the source. The statement that the shortfall is largest late in follow-up rests on `dev/audit/stage2e/m7_summary.csv` (coverage of the survival difference by month, cells 35–39), which the report's table does not show.
- **Section 8.** Nothing omitted. The Kaplan–Meier contributions were not cross-checked through an exported function (the score contributions were).
- **Pointers.** Each is one sentence appended to an existing paragraph; in `weightedsurv_examples.Rmd` it stands at the end of the last paragraph of the section "Weights treated as known or as estimated".
- **pkgdown.** The index entry was added and `check_pkgdown()` run; the site was not built.
