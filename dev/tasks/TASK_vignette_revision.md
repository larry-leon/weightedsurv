# TASK: vignette revision pass — weightedsurv `km-fixes-0.2.0`

Written 2026-11-19. Executor: Claude Code (CC). Approved by Larry (items 1–8 and OPEN-1..4 defaults).

## 0. Rules

- Branch `km-fixes-0.2.0`; HEAD must contain `ac76f82`. If not, stop.
- First action: copy this file to `dev/tasks/TASK_vignette_revision.md` and commit it.
- **Scope: `vignettes/` and `README.md` only** (plus `dev/` for the task doc, report and scratch). **No changes under `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE`.**
- Never push. Never merge.
- Gates (section 6) are stop-on-failure: if one fails, stop, write the report, do not ask.
- Verify everything from source at the current HEAD. Line numbers below are from an earlier inventory: locate each passage by content and confirm before editing.
- No point estimate moves from package changes; unweighted results stay byte-identical. Numbers change only because the example population changes (item 1) or the variance option changes (item 2).
- Compute every quoted number inline in the Rmd (inline R or a chunk). No hard-coded numerical results in prose.
- Report findings only: no recommendation, no follow-up task proposals.

## 1. Files

| File | Role |
|---|---|
| `vignettes/weightedsurv_examples.Rmd` | weighted examples (main) |
| `vignettes/articles/martingale_resampling_applications.Rmd` | weighted examples |
| `vignettes/articles/weightedcox_methodology.Rmd` | weighted examples + method text |
| `vignettes/articles/weightedkm_methodology.Rmd` | method text |
| `vignettes/articles/weightedrmst_methodology.Rmd` | method text |
| `README.md` (no `README.Rmd`) | validation statement |

## 2. API facts (verified at `ac76f82`; re-confirm from `R/`)

- `ps_weights(formula = NULL, data = NULL, type = c("ipw", "stabilized"), model = NULL)`. **Default is `"ipw"`; pass `type = "stabilized"`** to reproduce the vignettes' `sw.weights`. Returns class `"ps_weights"` (`weights`, `type`, `treat`, `ps`, `IF`, `dw`, `model`, `call`). `data` must be exactly the rows analysed. `print()` shows n, weight range, effective sample size by arm.
- `weight_variance = c("known", "estimated")` and `weight_model` are arguments of `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()`, `plotKM.band_subgroups()`.
- `get_dfcounting(...)` forwards `...` to `df_counting()`.
- `cox_rhogamma()`, `wlr_cumulative()`, `wlr_dhat_estimates()` inherit from the `df_counting()` result.
- `df_counting()$ess`: data frame `group`, `n`, `ess`, `ratio`; treated row first.
- `draws > 0` is not available with `"estimated"` in `df_counting()`.
- `plotKM.band_subgroups()` default `draws.band` is 1000; weighted `logrank_results` is the `coxph` robust score test.

Pattern:

```r
pw <- ps_weights(treat ~ age + size + grade3 + nodes + pgr + er, data = d, type = "stabilized")
d$sw.weights <- pw$weights
dc <- get_dfcounting(df = d, ..., weight.name = "sw.weights",
                     weight_variance = "estimated", weight_model = pw)
dc$ess
```

## 3. Items

### Item 1 — Rotterdam example: restricted population

In all three files with weighted examples, restrict the Rotterdam weighted example to **postmenopausal women without chemotherapy**: `meno == 1 & chemo == 0` among `nodes > 0` (538 control / 291 treated; n = 829).

- First show the effective sample size on **all node-positive women** (1,207 / 339) with the current model `treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er`, as the stated reason for restricting (treated ESS ratio about 0.31, largest weight about 12.9). Then restrict.
- Restricted-population model: `treat ~ age + size + grade3 + nodes + pgr + er`.
- Keep the unrestricted node-positive data under a separate name (e.g. `df_rotterdam_np`) — the cross-trial example needs it (item 3).
- One sentence of context is allowed: GBSG trial HR 0.69 (0.54–0.89).

### Item 2 — `weight_variance = "estimated"`

In the restricted Rotterdam example:

- Build weights with `ps_weights(..., type = "stabilized")`; replace the hand-rolled propensity code.
- Use `weight_variance = "estimated", weight_model = pw` for the log-rank test, hazard ratio, survival difference and RMST.
- Show `"known"` beside `"estimated"` **once**, as one small comparative table (z, HR CI, RMST CI), with a short reading.
- Keep `"known"` for FH(0,1) and for the cross-trial example.
- Wording about the two options (use this substance, do not present either as generally better): with evenly spread weights, `"estimated"` is near nominal and `"known"` is conservative; with concentrated weights both fall short, mostly late in follow-up, and the shortfall comes from the weights. Point to `dev/audit/REPORT_km_fixes_stage2e.md` for the evaluation.
- `check.seKM = TRUE, draws = 5000` chunk (examples Rmd, about line 712): `draws > 0` is not available with `"estimated"`. Keep that chunk on `"known"` and say so in one sentence.
- Weighted Cox table (examples Rmd about 628–692): `cox_rhogamma(..., draws = 1000)` with CIs from `hr_ci_star`. Under `"estimated"` use the asymptotic CI; keep the resampling CI only in a `"known"` row. Review the comment near line 632 (recommends asymptotic SEs with case weights) and make it accurate.

### Item 3 — Cross-trial example (GBSG treated vs Rotterdam node-positive control; 1,207 / 246)

Keep as a stated hard case.

- Controls come from the **unrestricted** node-positive Rotterdam data (item 1 rename).
- `weight_variance = "known"`.
- Show the ESS table (item 5).
- Add a caveat paragraph: weighting does not balance these trials (grade 3 is 20% of GBSG treated vs 74% of Rotterdam control, grade coded 1/2/3 in GBSG and 2/3 in Rotterdam; postmenopausal 76% vs 51%); treated ESS ratio about 0.30; few treated at risk late.
- Read results to **60 months**; correct the "comparable horizon" text (about 515–517) accordingly.

### Item 4 — Cumulative log-rank figure

Keep. Start the curve at the 10th event in each arm (about 9.4 months, computed from the data, not hard-coded); say so in the caption; set the y-axis to fit the plotted range (replace `ymin <- -1; ymax <- 5`). Correct the surrounding text (about 498–512), including the "flat for time points ≥ 80" sentence and the `abline(v = 80)`, to match what the figure shows and the 60-month reading.

### Item 5 — Effective sample size

After each set of weights, a two-row table by arm from `df_counting()$ess` (`group`, `n`, `ess`, `ratio`) and one sentence on reading it. No threshold.

### Item 6 — Bands

Leave the bands as built.

- Remove `qtau = 0.025` from every band call (it has no effect without `modify_tau = TRUE`). Known locations: examples Rmd 356, 379, 403, 531, 552, 613; martingale article 293, 336, 379; RMST article 286. Search all files.
- Correct text that says the band covers a trimmed / quantile range (examples 341–343 "qtau quartiles"; KM article 248; martingale 303, 419, 448): the band is over the whole range as called.
- KM article 252–255: text gives the unstandardised sup and band; the code uses the standardised version. Make the text match the code.

### Item 7 — Inaccurate text and quoted numbers

Read each passage at HEAD against the current `R/` source and correct it. Compute quoted numbers inline.

`weightedsurv_examples.Rmd`:
- 15–16 (OPEN-4): the studies compare **hormonal therapy vs none**, not "hormon replacement therapy in comparison to chemotherapy"; remove "draft vignette".
- 29–49 (OPEN-4): source files and helpers named — check each against the current `R/`; fix or remove.
- 148–150 (OPEN-4): install comment naming `weightedKMplots` and `github_install` — replace with the current install instructions.
- 703–713: `check.seKM` text.
- 719–768 (`adjustedCurves` comparison; OPEN-1): make the `adjustedCurves` call use the **same event variable** as the weightedsurv analysis (currently it uses `status`; the analysis uses `event <- status2` with `tte` from `rfstime`). State in the data-prep text which event/time definition is used. **Do not change the coding.** State that `adjustedCurves` `iptw_km` standard errors still differ, and that `survfit` weighted KM standard errors agree.
- 773–776, 793: comments naming `survdiff` as a source — weighted `logrank_results` is now the `coxph` robust score test.

`articles/weightedkm_methodology.Rmd`:
- 171 (OPEN-2): keep the truncation sentence as a description of practice; add that the package never truncates.
- 199 (OPEN-3): replace "additional variance … typically small" with a statement of what `weight_variance` does and a pointer to the evaluation.
- 203–207: weighted Greenwood formula and "no adjustment for estimated weights" → robust weighted KM variance as implemented (Stage 2c), plus the estimated-weights option.
- 209 ("preferred"), 220 ("(independent) variance estimates"), 234–242 and 251 (dN-based resampling description) → influence-based resampling as implemented.
- 259, 276, 326 (validation against `adjustedCurves` and `survfit`): correct per the current facts above.

`articles/martingale_resampling_applications.Rmd`: weights 174–185; weighted calls 196, 241, 274, 316, 421; band calls 276, 321, 362, 423; text 303, 419, 448; 426 and 436 use `"ps_weights"` as a column name in an unevaluated snippet — rename the column (e.g. `sw.weights`) so it does not collide with the function name.

`articles/weightedcox_methodology.Rmd`: 101–103, 105, 112–114 (score variance and resampling formulas) → score-residual variance as implemented (Stage 2d); weighted call 315; weights 329–335; 383.

`articles/weightedrmst_methodology.Rmd`: 112 says delta-method variance; the package uses the resampling variance — correct. 219–223, 225; 281 unevaluated call.

`README.md` 231: validation statement — correct per the current facts.

**Method text.** Write replacement text for: robust weighted KM variance; influence-based resampling; score-residual variance for the weighted log-rank / Cox; and one short estimated-weights section (what `ps_weights()` returns, what `"estimated"` adds, the limitation with `draws`). Sources, all in the repo: `dev/tasks/TASK_km_fixes_stage2e.md` §2, `dev/audit/REPORT_km_fixes_stage2{c,d,e}.md`, the Stage 2c and 2d task documents in `dev/tasks/`, and the `R/` source. Formulas must match the code; where text and code disagree, the code is right.

### Item 8 — Truncated weights

Remove the truncated-weight lines and any use of them (examples 306, 452; martingale 188; search for `trunc`).

## 4. Reference values (CC recomputes; these are targets for the gates, not text to paste)

Restricted Rotterdam (postmenopausal, no chemo):

| Quantity | Value |
|---|---|
| n; control / treated | 829; 538 / 291 |
| Treated ESS ratio (control) | 0.93 (0.98) |
| Largest weight | 2.7 |
| Weighted HR, known CI | 0.66 (0.56–0.80) |
| Weighted HR, estimated CI | 0.56–0.78 |
| z, known → estimated | 4.45 → 4.87 |
| Unweighted HR | 0.72 (0.60–0.86) |
| Worst imbalance after weighting | 0.05 (from 0.26) |

All node-positive (shown only as the reason): ESS treated 339 → 105.2 (0.310); control 1,207 → 1,162.9 (0.963); largest weight 12.9.

Cross-trial: ESS treated 246 → 73.0 (0.297); control 1,098.8 (0.910); largest weight 13.0. HR 0.6601 (0.4445, 0.9803) known; z 2.110. Cut at 60 months: HR 0.62 (0.41–0.95), z 2.30. 10th treated event at 9.4 months.

Stored values that changed with Stage 2e (context for any text quoting them): weighted `cox_results$score` 34.70 → 12.06; `logrank_results$chisq` 4.771 → 12.07 (Rotterdam, old population), 13.34 → 4.438 and 15.43 → 4.453 (cross-trial).

## 5. Order of work

1. Copy and commit this task doc.
2. Before any edit: knit all vignettes at HEAD and save every figure to `dev/scratch/vignette_before/` (the "before" set).
3. Read every listed passage at HEAD; record in the report any line that is not where the inventory says or does not say what the inventory says.
4. Items 1, 2, 5, 8 in `weightedsurv_examples.Rmd`; then 3, 4; then 6, 7.
5. Same items in the martingale and Cox articles.
6. Method text in the KM, Cox, RMST articles; `README.md`.
7. Gates (section 6).
8. Before/after figure sheet; report; commit.

Commit per file or per item, messages prefixed `vignettes:`.

## 6. Gates (stop-on-failure; each is a machine check, record PASS/FAIL with the evidence)

| # | Check |
|---|---|
| G1 | `git diff --name-only ac76f82..HEAD` lists nothing outside `vignettes/`, `README.md`, `dev/`. |
| G2 | All vignettes and articles knit without error or warning-as-error. |
| G3 | `R CMD check` (as-cran) clean: 0 errors, 0 warnings; notes no worse than at `ac76f82`. |
| G4 | `grep -rn "qtau = 0.025" vignettes/ README.md` returns nothing. |
| G5 | `grep -rn "sw.weights_trunc\|_trunc" vignettes/` returns nothing. |
| G6 | `ps_weights(` and `weight_variance = "estimated"` each present in all three example files. |
| G7 | Restricted population, computed in the knitted vignette: n = 829; 538 / 291; treated ESS ratio in [0.92, 0.94]; largest weight in [2.6, 2.8]; weighted HR in [0.655, 0.675]; z known in [4.40, 4.50], estimated in [4.82, 4.92]. |
| G8 | Cross-trial: 1,207 / 246; treated ESS ratio in [0.29, 0.30]; HR known in [0.655, 0.665]; `weight_variance` is `"known"`. |
| G9 | Unweighted results: unweighted `df_counting()` outputs on the unrestricted node-positive data identical (`identical()` or `all.equal(tolerance = 0)`) before and after. Test suite still 1,691 passing. |
| G10 | No `"ps_weights"` used as a column name anywhere in `vignettes/`. |
| G11 | No hard-coded result numbers in changed prose: list every numeral in changed text lines and mark each as inline-computed, a design constant, or a literature value (GBSG 0.69 (0.54–0.89)). |
| G12 | Every figure regenerated; before/after sheet built. |

## 7. Report

`dev/audit/REPORT_vignette_revision.md`:

1. Gate table (G1–G12, PASS/FAIL, evidence).
2. One comparative table: restricted Rotterdam and cross-trial, columns = reference value (section 4) / vignette value, rows = quantities in section 4.
3. Per-file change list: passage, what it said, what it says now, item number.
4. Inventory discrepancies from step 5.3.
5. Anything in `vignettes/` or `README.md` found inaccurate but not covered by an item: listed, not changed.

Before/after figure sheet: one PDF or HTML, each figure before (from `dev/scratch/vignette_before/`) beside after, labelled by file and chunk.

After the final commit, copy the report and the figure sheet to `~/Downloads/` and print the paths. Do not push.
