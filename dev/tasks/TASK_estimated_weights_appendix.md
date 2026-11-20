# TASK: estimated-weights statistical appendix — weightedsurv `km-fixes-0.2.0`

Written 2026-11-20. Executor: Claude Code (CC). Approved by Larry.

## 0. Rules

- Branch `km-fixes-0.2.0`; HEAD must contain `a15ff33`. If not, stop.
- First action: copy this file to `dev/tasks/TASK_estimated_weights_appendix.md` and commit it.
- **Scope: `vignettes/`, `_pkgdown.yml` (article index entry only), `dev/`.** **No changes under `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE`, `README.md`.**
- Never push. Never merge.
- Gates (section 5) are stop-on-failure: if one fails, stop, write the report, do not ask.
- Verify from source at the current HEAD. **Where the Stage 2e task document and the code disagree, the code is right**; record the disagreement in the report.
- No simulation, no bootstrap, no refitting loop in this task or in the article. The only compute is knitting.
- Report findings only: no recommendation, no follow-up task proposals.

## 1. Deliverable

A new standalone article: `vignettes/articles/estimated_weights_methodology.Rmd`

Title: "Estimated propensity-score weights: statistical details and properties".

Purpose: one place that gives the derivation and properties behind `weight_variance = "estimated"`, for the weighted Kaplan-Meier curve, the survival difference, RMST, and the weighted log-rank / Cox score. The KM, Cox and RMST methodology articles and the main vignette each get **one sentence** pointing to it (section 3).

Style: same notation, header and reference style as the three existing methodology articles. Mathematical, self-contained, no narrative filler.

## 2. Sources (all in the repo)

- `dev/tasks/TASK_km_fixes_stage2e.md` §2 (method and reference code for estimated weights)
- `dev/audit/REPORT_km_fixes_stage2{c,d,e}.md` and the Stage 2c / 2d task documents in `dev/tasks/`
- `R/` source: `ps_weights()`, `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()`, `plotKM.band_subgroups()`, `cox_rhogamma()`, `ws_score_psi()`, `score_calculation()`, and whatever they call
- The estimated-weights passages CC wrote in the previous pass (`weightedkm_methodology.Rmd`, `weightedcox_methodology.Rmd`, `weightedrmst_methodology.Rmd`)

## 3. Content of the article

Each section states what the package computes. Every displayed formula must correspond to code (section 4, formula–code map).

**1. What the two options mean.**
`"known"`: weights treated as fixed constants. `"estimated"`: the sampling variability of the fitted propensity model is carried into the variance by a linearization (influence-function correction). State plainly what resampling does and does not do: the multipliers are applied to the corrected per-subject contributions; **the propensity model is not refitted within a draw**. No point estimate differs between the two options.

**2. Setup and notation.**
Data, treatment indicator, covariates, counting and at-risk processes, the two weight types (`"ipw"`, `"stabilized"`), the estimand the weights target as implemented.

**3. The propensity model.**
The model `ps_weights()` fits (as coded); its score and information; the per-subject influence function for the coefficients (`IF`); the derivative of each subject's weight with respect to the coefficients (`dw`) for both weight types. For stabilized weights, state whether the estimation of the marginal treatment probability is included in the correction, as the code does it.

**4. Linearization, estimator by estimator.**
- Weighted Kaplan-Meier in one arm: per-subject contribution ψ_i(t) under known weights (Stage 2c); the added term under estimated weights; the corrected contribution Φ_i(t).
- Survival difference: the corrected contribution of the difference; why the two arms are no longer independent (one propensity model feeds both); the variance is not the sum of the arm variances.
- RMST and cumulative RMST difference: integral of the above.
- Weighted log-rank / Cox score: ψ_i(β) under known weights (Stage 2d); the corrected ψ^ew_i; z-score, standard error of log HR, the asymptotic interval. How Fleming–Harrington and the other time-weight schemes enter, and which are covered.

**5. Variance.**
Variance as the sum of squared corrected contributions, for each estimator. Relation to the `"known"` variance. Why `"known"` tends to be conservative when the propensity model is correct: give the argument (projection of the known-weight contribution onto the propensity score space) and label it as an argument or a proof, whichever it is.

**6. Multiplier resampling with corrected contributions.**
The resampled process for the survival difference and for RMST; the standardized supremum and the band as coded; the resampled score in `cox_rhogamma()`. What justifies it (conditional multiplier central limit theorem applied to the linear representation, with estimated contributions replacing true ones) and what conditions that needs. Why `draws > 0` for the pointwise variance in `df_counting()` is not available under `"estimated"`, from the source / Stage 2e record.

**7. Conditions and properties.**
Correct propensity model; overlap; sample size relative to weight concentration (effective sample size as the diagnostic, no threshold). Then the evaluation evidence: **one comparative table** copied from `dev/audit/REPORT_km_fixes_stage2e.md` (coverage of the survival difference and log HR, null rejection, by regime, known vs estimated, with the refitting-bootstrap row), labelled with its source, followed by a short plain-language reading. Do not recompute it. Do not present either option as generally better on the strength of one regime: with even weights `"estimated"` is near nominal and `"known"` conservative; with concentrated weights both fall short, mostly late in follow-up, and `"estimated"` agrees with the refitting bootstrap, so the shortfall comes from the weights.

**8. Numerical illustration (evaluated, fast).**
On the restricted Rotterdam population of the main vignette (`nodes > 0 & meno == 1 & chemo == 0`, model `treat ~ age + size + grade3 + nodes + pgr + er`, `type = "stabilized"`): rebuild the `"estimated"` standard error by hand from the returned pieces (`pw$IF`, `pw$dw`, and the per-subject contributions) for (a) the survival difference at one time point and (b) log HR under FH(0,0), and show each beside the package value. Use exported functions and returned objects only (no `:::`). If the pieces needed are not returned by exported functions, omit the part that cannot be done, say so in the article in one sentence, and record it in the report. Also show the ESS table (`df_counting()$ess`). All numbers inline or in chunks.

**9. Limits and summary.**
What is covered and what is not (model class, weight types, time-weight schemes, `draws`, stratified results), as the code has it. Pointers to the three methodology articles and the main vignette.

**References.** Reuse entries already in the vignettes where they exist. Candidates for the known results, to include only if CC can confirm the full bibliographic details (do not cite from memory): Lunceford & Davidian (2004, Statistics in Medicine); Hirano, Imbens & Ridder (2003, Econometrica); Xie & Liu (2005, Statistics in Medicine); Lin & Wei (1989, JASA); Lin, Wei & Ying (1993, Biometrika). List every reference added in the report.

**Pointers (one sentence each, nothing else changed in these files):**
- `vignettes/weightedsurv_examples.Rmd`, section "Weights treated as known or as estimated"
- `vignettes/articles/weightedkm_methodology.Rmd`, estimated-weights passage
- `vignettes/articles/weightedcox_methodology.Rmd`, section "Estimated weights"
- `vignettes/articles/weightedrmst_methodology.Rmd`, weighted-setting passage

**pkgdown.** If `_pkgdown.yml` has an articles index, add the new article to it beside the other methodology articles. Change nothing else in that file.

## 4. Standards

- **Formula–code map** (goes in the report): one row per displayed formula or defined symbol in the article → file, function, object in `R/`. A formula with no counterpart in the code is either removed or explicitly labelled in the article as background not computed by the package.
- **Argument vs proof.** Every asymptotic claim is labelled in the text as proved here, cited, or a heuristic argument. List the heuristic ones in the report.
- No hard-coded result numerals in prose, except the copied evaluation table (labelled with its source).
- Notation consistent with the three existing methodology articles; where they differ among themselves, follow the KM article and note it in the report.

## 5. Gates (stop-on-failure; record PASS/FAIL with evidence)

| # | Check |
|---|---|
| G1 | `git diff --name-only a15ff33..HEAD` lists nothing outside `vignettes/`, `_pkgdown.yml`, `dev/`. `git diff a15ff33..HEAD -- _pkgdown.yml` is at most the index entry. |
| G2 | `rmarkdown::render()` with `options(warn = 2)`: `OK` for the new article and the four files that received a pointer. |
| G3 | New article render time under 60 seconds. `grep` of the new file finds no `replicate(`, `boot`, `for (b in`, `mclapply`, `future` loop that refits a model. |
| G4 | Section 8 by-hand values agree with the package values: absolute difference below 1e-8 for each quantity shown. (If a part was omitted under the exported-only rule, record which and why; the gate applies to what is shown.) |
| G5 | The pointer is present in all four files: `grep -l "estimated_weights_methodology" ` returns the four. The diff of each of those four files is one added or changed sentence. |
| G6 | `pkgdown::check_pkgdown()` passes (or, if there is no `_pkgdown.yml` articles index, record that). |
| G7 | `R CMD check --as-cran`: 0 errors, 0 warnings, NOTEs no worse than at `a15ff33` in the same session. Tests `PASS 1691`, `FAIL 0`. |
| G8 | Formula–code map complete: every displayed equation in the article has a row; no row marked "not found" unless the article labels that formula as background. |
| G9 | Numerals in the article's prose: each is inline-computed, a design constant, a formula constant, or in the copied evaluation table. |

## 6. Report

`dev/audit/REPORT_estimated_weights_appendix.md`:

1. Gate table.
2. Formula–code map.
3. One table: section 8 quantities, by-hand value / package value / difference.
4. Statements labelled as heuristic arguments; statements cited; statements proved in the article.
5. Disagreements between `dev/tasks/TASK_km_fixes_stage2e.md` §2 and the code, if any.
6. References added, with full details.
7. Anything omitted from section 3 of this task, and why.

After the final commit, copy the report and the rendered HTML of the new article to `~/Downloads/` and print the paths. Do not push.
