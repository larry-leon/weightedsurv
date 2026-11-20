# weightedsurv 0.2.0: vignette revision pass on branch `km-fixes-0.2.0`

Task: `dev/tasks/TASK_vignette_revision.md`. Start: `ac76f82` (in HEAD). Scope kept: only `vignettes/`, `README.md` and `dev/` changed; nothing under `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE`. No push, no merge.

**Top line.**

- **Gates.** G1–G12 pass. The sixth article, `articles/weightedcox_wald-vs-logrank_simulations.Rmd`, is not in the task's file list and is byte-identical to `ac76f82`; it knits without error at full size (72 minutes; SIMS-NOTE under the gate table).
- **Items.** Items 1–8 are done in every file the task names. Choices the task left open, and two passages changed outside the inventory, are listed under "Execution notes".
- **Numbers.** Every reference value of section 4 of the task is reproduced by the knitted vignette (section 2 below). One displays differently at three decimals: the control ESS ratio on all node-positive women is 0.96350 and prints as 0.964 (task: 0.963).
- **Unchanged results.** The unweighted `df_counting()` result on the unrestricted node-positive data is `identical()` before and after. The eight figures of the main vignette that show unweighted or cross-trial analyses (other than the redrawn cumulative log-rank figure) are byte-identical image files before and after.
- **Deliverables.** This report; the figure sheet `dev/audit/vignette_revision/figures_before_after.html` (25 figure pairs, and the 3 figures of the simulations article); the gate output in one file, `dev/audit/vignette_revision/gates_output.txt`; scripts and figure sets in `dev/scratch/`.

Commits (each prefixed `vignettes:`), in order: `20262b2` task document; `b168581` `weightedsurv_examples.Rmd`; `8cd7f4c` martingale article; `52827de` Cox methodology article; `25c7c53` KM methodology article; `25e5ee8` RMST methodology article; `a15ff33` `README.md`; `79f2330` this report, the gate scripts, the figure sets and the figure sheet; then one commit that adds the result of the full-size knit of the simulations article to this report and to the figure sheet. `a15ff33` is the last commit that touches `vignettes/` or `README.md`.

Execution notes (findings, not recommendations):

- **Environment.** R 4.5.2, knitr 1.52, rmarkdown 2.32, survival 3.8.11, adjustedCurves 0.11.4. The package was installed once from the working tree into a temporary library (`R/` is identical at `ac76f82` and at HEAD) and both figure sets were knitted against that build.
- **Before set.** Knitted at `ac76f82` before any edit (`dev/scratch/knit_vignettes.R`, `knitr::knit`, png at 96 dpi): 18 + 6 + 1 figures (main vignette, martingale article, Cox article); the KM and RMST articles have no figure. Stored in `dev/scratch/vignette_before/`. The objects of each knit (`env.rds`, up to 128 MB) are used by the gate script and are not committed (`dev/scratch/.gitignore`).
- **Known beside estimated, once.** The comparative table (z, hazard-ratio interval, RMST interval) is in the main vignette only (section "Weights treated as known or as estimated"). The martingale article points to it; the Cox article shows FH(0,0) under `"estimated"` and FH(0,1) under `"known"` without a second table.
- **Cross-trial weights.** The cross-trial weight column is still computed by hand (truncated-weight lines removed), not by `ps_weights()`: `ps_weights()` differs from that column by 4.1e-08 (Stage 2e report, D2), and the cross-trial numbers were to stay as they were. `weight_variance = "known"` is passed explicitly in the cross-trial calls.
- **Weighted Cox table, schemes the task does not name.** FH(0,0) has two rows: `"estimated"` with the asymptotic interval, `"known"` with the resampling interval. FH(0,1) is `"known"` with the resampling interval, as the task says. FH_exp2, MB(12) and zero_one(6) are not named in the task; they are shown under `"known"` with the resampling interval, as FH(0,1).
- **Bands left as built.** No band call changed other than the removal of `qtau = 0.025` and, in the restricted Rotterdam example, the two `weight_variance` arguments. The 60-month reading of the cross-trial example is done in text (inline values at 60 months from the plotted fits) and in one table with follow-up cut at 60 months; the plots still run to the smaller of the two arms' last event times.
- **Checks section of the main vignette.** The two `check.seKM` chunks used to overwrite `dfcount_rotterdam_wtd`; they now write `dfcount_rotterdam_known` and `dfcount_rotterdam_draws`, so the `"estimated"` object is not overwritten. The standard-error comparison with `adjustedCurves` uses the analytic `"known"` object (at `ac76f82` it used the `draws = 5000` object, because that chunk ran last).
- **`adjustedCurves` call.** It now uses `ev_time = "tte"`, `event = "event"` (the variables of the weightedsurv analysis; `tte` equals `time_months`), the restricted population, and the same model formula as `ps_weights()` (`grade3`; the call at `ac76f82` had `grade`, which gives the same fit because Rotterdam grade takes the values 2 and 3). No coding was changed.
- **GBSG context value.** The sentence of context is computed inline from the vignette's own GBSG object (`cox_rhogamma()`, FH(0,0), asymptotic interval) and prints as 0.69 (0.54, 0.89), the value the task gives.
- **Two passages changed outside the inventory.** (a) KM article, lead sentence of "Practical Weight Computation" (line 165 at `ac76f82`): it said the stabilized weights are computed "as implemented in the `weightedsurv_examples` vignette"; item 2 replaced that hand-written code by `ps_weights()`, so the sentence now names `ps_weights()`. (b) RMST article, line 178 at `ac76f82`: it said the band interval is determined "from the observed quantiles of the event times (via the `qtau` parameter)"; this is the class of text item 6 corrects, at a location item 6 does not list.
- **Additions.** A `ps_weights()` bullet in the capability list of the main vignette and in the function list of the KM article; a two-column table of standardized differences before and after weighting in the main vignette (the "worst imbalance" row of section 4 of the task is computed from it); chunk labels on the figure chunks that changed; a small formatting helper `fmt()` for the inline numbers.
- **R CMD check NOTE count.** The check of this task shows 3 NOTEs where the stored Stage 2e log shows 2. The third is "checking for future file timestamps ... unable to verify current time". The same check run on `ac76f82` in the same session shows the same 3 NOTEs (G3).

---

## 1. Gates

| # | Result | Evidence |
|---|---|---|
| G1 | PASS | `git diff --name-only ac76f82..HEAD` lists files under `vignettes/` (5), `README.md` (1) and `dev/` only; `grep -v '^vignettes/\|^README.md$\|^dev/'` on that list returns nothing. |
| G2 | PASS | `dev/scratch/render_vignettes.R`: `rmarkdown::render()` with `options(warn = 2)` prints `OK` for `weightedsurv_examples.Rmd` (17 s), `martingale_resampling_applications.Rmd`, `weightedcox_methodology.Rmd`, `weightedkm_methodology.Rmd`, `weightedrmst_methodology.Rmd`. `knitr::knit` of the same five: no error, 0 warnings escaping the chunks. The vignette rebuild inside `R CMD check` is OK (16 s). Sixth article: `knitr::knit` of the file as written, OK (SIMS-NOTE). |
| G3 | PASS | `R CMD check --as-cran --timings` on the tarball built from a clean `git clone` at `a15ff33` (`RSTUDIO_PANDOC` set to the quarto tools directory, `DISPLAY` unset): `Status: 3 NOTEs`, 0 errors, 0 warnings. NOTEs: future file timestamps ("unable to verify current time"); top-level files (pandoc not on PATH for README/NEWS); HTML version of the manual (old HTML Tidy). The identical procedure on `ac76f82` in the same session: `Status: 3 NOTEs`, the same three. The stored log of Stage 2e (`dev/audit/stage2e/00check_2e.log`) has the last two only. Logs: `dev/audit/vignette_revision/00check_a15ff33.log`, `00check_ac76f82_same_session.log`. |
| G4 | PASS | `grep -rn "qtau = 0.025" vignettes/ README.md` returns nothing (exit status 1). |
| G5 | PASS | `grep -rn "sw.weights_trunc\|_trunc" vignettes/` returns nothing (exit status 1). |
| G6 | PASS | Occurrences of `ps_weights(` / `weight_variance = "estimated"`: `weightedsurv_examples.Rmd` 5 / 10; `martingale_resampling_applications.Rmd` 5 / 9; `weightedcox_methodology.Rmd` 3 / 3. |
| G7 | PASS | `dev/scratch/gates.R`, from the objects of the knitted files, the same in all three: n = 829; 538 / 291; treated ESS ratio 0.9301; largest weight 2.7090; weighted HR 0.6648; z known 4.4491, estimated 4.8686; the vignette's weighted object has `weight_variance = "estimated"`. The table printed by the main vignette shows z 4.45 (known) and 4.87 (estimated). |
| G8 | PASS | Same script: 1,207 / 246; treated ESS ratio 0.2967; HR known 0.6601 (0.4445, 0.9803); `weight_variance` is `"known"` in the `df_counting()` result, the `cox_rhogamma()` result, the band fit, the RMST result and the cut-at-60 result. |
| G9 | PASS | (a) Same script: the node-positive data are `identical()` before and after (1,546 rows, 29 columns); the unweighted `df_counting()` result on them (`get_dfcounting(..., by.risk = 24)`, the call of `ac76f82`) is `identical()` to the object of the knit at `ac76f82` in all 47 elements (the `coxph` object compared by coefficients and variance at tolerance 0); z = 2.184188. The vignette's GBSG object `dfcount_gbsg` is `identical()` before and after as well. (b) Test suite inside `R CMD check` at `a15ff33`: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]` (`dev/audit/vignette_revision/testthat_a15ff33.Rout`). |
| G10 | PASS | `grep -rn '"ps_weights"' vignettes/` returns nothing. Every other occurrence of `ps_weights` in `vignettes/` is the function (`ps_weights(` or `` `ps_weights` `` in a function list). The snippet column is now `sw.weights`. |
| G11 | PASS | `dev/scratch/g11_numerals.py` → `dev/audit/vignette_revision/g11_numerals.txt` (exit status 0): of the 545 lines added or changed in the six files, the prose lines hold 109 numerals written literally (inline R removed first); 0 unclassified; none is a result. The 11 comment lines added inside R chunks that contain a numeral are listed in the same file; none holds a result. Summary below. |
| G12 | PASS | 25 figures regenerated at HEAD (`dev/scratch/vignette_after/`): 18 main vignette, 6 martingale article, 1 Cox article; the KM and RMST articles draw none; plus the 3 figures of the unchanged simulations article (one knit, SIMS-NOTE). Sheet: `dev/audit/vignette_revision/figures_before_after.html`, 25 pairs and those 3, each labelled by file and chunk (before and after labels) and marked identical or different. Byte-identical image files: 8 of 18 (main vignette: flowchart, weight functions, the two GBSG K-M figures, the GBSG / cross-trial K-M figure, the GBSG / cross-trial band and RMST figure, the GBSG subgroup figure, the GBSG er > 0 figure), 1 of 6 (martingale article: GBSG K-M), 1 of 1 (Cox article: weight functions). |

**SIMS-NOTE (`articles/weightedcox_wald-vs-logrank_simulations.Rmd`).** The file is unchanged by this task (`git diff ac76f82 HEAD -- <file>` is empty) and passes no subject weights. As written it runs 6 scenarios × 5,000 simulated trials on 12 workers inside the knit (`eval = TRUE`; the pkgdown page in `docs/` was built with 10 trials per scenario).
- A copy with two lines changed (`n_sim = 10`, `num_workers = 2`), rendered with `rmarkdown::render()` and `options(warn = 2)`: `OK`, 68 s, three figures.
- The file as written, knitted with `knitr::knit` (`dev/scratch/knit_vignettes.R`, started when the "before" set was made): `OK`, no error, 0 warnings escaping the chunks, 3 figures, 4,314 s. Its own output: "Completed: 30000 | Failed: 0", "Elapsed: 71.8 minutes". Because the file is byte-identical at the two commits it was knitted once; the figure sheet shows its three figures with the same image on both sides. Simulation results in this knit depend on the run and were not compared with anything.

**G11, classes of the 109 numerals** (full list with line numbers and context in `g11_numerals.txt`):

| Where | Count | Class |
|---|---|---|
| in math (`$...$`, `$$...$$`) | 65 | formula constants: indices, exponents, the 10⁻¹⁰ tolerance of the formulas |
| in code spans | 5 | design constants: argument values or names (`draws > 0`, `1 - qtau`) |
| plain text | 10 | design constants: names of weighting schemes (FH(0,0), FH(0,1), MB(12), zero_one(6)) |
| plain text | 4 | design constant: "ratio near 1" / "close to 1" (reference value of the ESS ratio) |
| plain text | 3 | design constant: confidence level (95%, 0.95) |
| plain text | 2 | design constant: "first 20" resampled curves (fixed in `KM_diff()`) |
| plain text | 3 | design constants: `draws` values (1000 in the chunk that follows; 5000 and 500 in a sentence of `ac76f82`) |
| plain text | 2 | design constants: at-risk time points 48 and 72 months (the counts are inline) |
| plain text | 2 | design constant: treatment coded 0/1 |
| plain text | 1 each | "10 events" (rule of item 4); "Grade 3" (category name); "`draws` must be 0" |
| plain text | 2 | wording of `ac76f82` kept in a rewritten paragraph: "a term is 0 when", "a weight of 2 is not two copies" |
| plain text | 5 | references to numbered steps of the text |
| plain text | 3 | list-item numbers of the algorithm in the KM article |
| inline R | not counted | every result quoted in changed prose: counts, ESS and ratios, largest weights, hazard ratios and intervals (the GBSG value included), z, standardized differences, percentages, months, standard-error differences and ratios |

By file: main vignette 25, martingale article 5, Cox article 34, KM article 42, RMST article 3, README 0.

## 2. Reference values (section 4 of the task) and vignette values

Vignette values are those computed in the knitted main vignette (objects read by `dev/scratch/gates.R`; the display in the text or table is given where it differs in rounding).

| Quantity | Reference | Vignette |
|---|---|---|
| **Restricted Rotterdam (postmenopausal, no chemotherapy)** | | |
| n; control / treated | 829; 538 / 291 | 829; 538 / 291 |
| Treated ESS ratio (control) | 0.93 (0.98) | 0.930 (0.984) |
| Largest weight | 2.7 | 2.709 |
| Weighted HR, known CI | 0.66 (0.56–0.80) | 0.6648 (0.5559, 0.7951); shown 0.66 (0.56, 0.80) |
| Weighted HR, estimated CI | 0.56–0.78 | (0.5639, 0.7838); shown 0.66 (0.56, 0.78) |
| z, known → estimated | 4.45 → 4.87 | 4.4491 → 4.8686; shown 4.45, 4.87 |
| Unweighted HR | 0.72 (0.60–0.86) | 0.7167 (0.6007, 0.8552); shown 0.72 (0.60, 0.86) |
| Worst imbalance after weighting | 0.05 (from 0.26) | 0.053 (from 0.259); shown 0.05, 0.26 |
| **All node-positive (shown as the reason)** | | |
| ESS treated | 339 → 105.2 (0.310) | 339 → 105.2 (0.310) |
| ESS control | 1,207 → 1,162.9 (0.963) | 1,207 → 1,162.9 (0.96350; the table shows 0.964) |
| Largest weight | 12.9 | 12.88; shown 12.9 |
| **Cross-trial** | | |
| n control / treated | 1,207 / 246 | 1,207 / 246 |
| ESS treated | 246 → 73.0 (0.297) | 246 → 73.0 (0.297); text: 73 of 246, ratio 0.30 |
| ESS control | 1,098.8 (0.910) | 1,098.8 (0.910) |
| Largest weight | 13.0 | 12.98; shown 13.0 |
| HR, known | 0.6601 (0.4445, 0.9803) | 0.6601 (0.4445, 0.9803); shown 0.66 (0.44, 0.98) |
| z | 2.110 | 2.110; shown 2.11 |
| Cut at 60 months: HR | 0.62 (0.41–0.95) | 0.62 (0.41, 0.95) |
| Cut at 60 months: z | 2.30 | 2.30 |
| 10th treated event | 9.4 months | 9.40 months; shown 9.4 |
| **Stored values that changed with Stage 2e** | weighted `cox_results$score`, `logrank_results$chisq` | not quoted or printed anywhere in the vignettes, at `ac76f82` or now |

Other values the revised main vignette displays, for the record: RMST difference at 118.0 months 15.6 (9.2, 22.0) known and 15.6 (9.7, 21.4) estimated; cross-trial at 60 months, survival difference 0.19 (0.08, 0.31) and cumulative RMST difference 5.3 (0.2, 10.5) months; treated patients at risk in the cross-trial data 104, 58 and 18 at 48, 60 and 72 months; cumulative log-rank z 2.30 at 60 months and 2.11 at the end; GBSG treated follow-up ends at 87.4 months, Rotterdam control at 230.9; grade 3 in 20% (GBSG treated) and 74% (Rotterdam control); postmenopausal 76% and 51%; weighted K-M standard errors against `survfit()`: largest absolute difference 8e-17 (control), 3e-17 (treated); against `adjustedCurves` `iptw_km`: curves within 7e-16, standard-error ratio 1.00–1.06 (control) and 0.90–1.05 (treated).

## 3. Changes by file

Line numbers are those at `ac76f82`.

### `vignettes/weightedsurv_examples.Rmd`

| Passage | What it said | What it says now | Item |
|---|---|---|---|
| 15–16 | "In this draft vignette … both evaluating hormon replacement therapy in comparison to chemotherapy" | "In this vignette … both comparing hormonal therapy with no hormonal therapy" | 7 |
| 26–50 capability list | plain names of files and helpers; "two or more groups"; `get_weights`/`extract_and_calc_weights` compute weights "for each subject and time point"; "Maggir-Burman" (also line 58) | each function and file checked against `R/` and written as code with its file; two groups (treatment 0/1); time weights at each time point; "Magirr-Burman"; a `ps_weights()` bullet added | 7 |
| 148–150 | "# to install weightedKMplots / library(devtools) / github_install(...)" | `install.packages("weightedsurv")` and `devtools::install_github("larry-leon/weightedsurv")`, as comments | 7 |
| setup chunk | — | helper `fmt()` for inline numbers | — |
| 263–266 | `df_rotterdam <- subset(gbsg_validate, nodes > 0)` | `df_rotterdam_np` holds the node-positive data; a paragraph states the time and event definition (`tte` from `rfstime`, `event = status2`) and the number of rows in which `status2` differs from `status` (inline) | 1, 7 |
| 268–294 | baseline table on `df_rotterdam` | the same table on `df_rotterdam_np` (unchanged content) | 1 |
| 297–308 | hand-written propensity code on all node-positive women; truncated weights at 306 | new section "Effective sample size with all node-positive women": `ps_weights(type = "stabilized")` with the earlier model, ESS table, a paragraph with the treated ESS ratio and largest weight (inline) as the reason for restricting; new section "Restricted population": `df_rotterdam <- subset(df_rotterdam_np, meno == 1 & chemo == 0)`, `ps_weights()` with `treat ~ age + size + grade3 + nodes + pgr + er`, printed object, n by arm (inline), table of standardized differences before and after weighting | 1, 2, 5, 8 |
| 309–316 | weighted `get_dfcounting()` with `weight.name` only | adds `weight_variance = "estimated", weight_model = pw`; ESS table and one sentence; weighted and unweighted hazard ratios and the GBSG value (all inline) | 2, 5 |
| 341–343 | "for timepoints t within "qtau" quartiles of the event times (i.e., for qtau = 2.5% …)" | the first graph shows the first 20 standardized resampled curves; definition of the band and its critical value; the band is over the whole time range of the call, no quantile range trimmed | 6 |
| 346–363, 369–417 | Rotterdam band and RMST calls with `qtau = 0.025`, weights as known | `qtau` removed (also from the GBSG call); `weight_variance = "estimated", weight_model = pw` in the Rotterdam band and RMST calls; chunk labels | 2, 6 |
| new, after 417 | — | section "Weights treated as known or as estimated": one table (z, hazard ratio and interval, RMST difference and interval) under both options, a reading, the wording of item 2 on the two options and the pointer to `dev/audit/REPORT_km_fixes_stage2e.md` | 2 |
| 420–453 | controls from `df_rotterdam`; truncated weights at 452 | opening paragraph (hard case; controls from `df_rotterdam_np`; weights as known); controls from `df_rotterdam_np`; truncated-weight lines removed | 3, 8 |
| 476–481 | weighted call | `weight_variance = "known"` explicit; ESS table and one sentence; caveat list (grade 3 shares and grade coding, postmenopausal shares, treated ESS ratio, treated at risk at 48, 60, 72 months), all inline; results read to 60 months | 3, 5 |
| 498–512 | "the log-rank test restricts comparisons to where both groups have follow-up"; "Notice how it is flat for time points ≈ ≥ 80"; `ymin <- -1; ymax <- 5`; `abline(v = 80)` | the statistic accumulates only while both groups are at risk; the curve is drawn from the 10th event in each arm (time computed, in text and caption); y-axis set from the plotted range; dashed line at 60 months with the z there; the curve is constant after the treated follow-up ends (time and value inline) | 4 |
| 515–517 | "With such dramatic difference in follow-up / It may make more sense to compare K-M differences and RMSTs over comparable 'horizon'" | "Reading the cross-trial comparison at a common horizon": follow-up of the two arms, where the plots end, why they are read at 60 months; a table of z and hazard ratio on all follow-up and cut at 60 months; after the figure, the survival difference and RMST difference at 60 months (inline) | 3 |
| 520–567 | band and RMST calls with `qtau = 0.025` | `qtau` removed; `weight_variance = "known"` explicit in the cross-trial calls | 3, 6 |
| 613 | `qtau = 0.025` in the er > 0 call | removed | 6 |
| 624–694 | five `cox_rhogamma(draws = 1000)` fits on the (known) weighted object, intervals from `hr_ci_star`; comment "with case-weights … recommend asymptotic versions of SEs" | a helper builds the rows; columns `weights` and `interval` added; FH(0,0) estimated / asymptotic, FH(0,0) known / resampling, the other four schemes known / resampling; the comment now says what the two standard errors are and where `hr_ci_asy` and `hr_ci_star` are centred; one introductory paragraph | 2, 7 |
| 697–714 | two headings, two chunks, the note on `draws = 5000` | the checks are stated to be on `"known"`; what `check.seKM` compares; the largest absolute difference from `survfit()` by arm (inline); `draws > 0` is not available with `"estimated"`, so that check stays on `"known"`; objects renamed | 2, 7 |
| 716–768 | `adjustedsurv(..., ev_time = "time_months", event = "status", ...)`, model with `grade`, on all node-positive women | `ev_time = "tte"`, `event = "event"`, restricted population, the model of `ps_weights()`; a sentence before and a paragraph after: curves agree, `iptw_km` standard errors still differ (ratio ranges inline), `survfit()` standard errors agree | 7 |
| 773–776, 793 | comments naming `wlr_estimates`, `survdiff`, `z_score_calculations` | the three sources as they are in the result (`lr^2 / sig2_lr`, `logrank_results$chisq`, `z.score^2`); `logrank_results` is `survdiff` without subject weights and the `coxph` robust score test with subject weights that are not all equal | 7 |

### `vignettes/articles/martingale_resampling_applications.Rmd`

| Passage | What it said | What it says now | Item |
|---|---|---|---|
| setup | — | helper `fmt()` | — |
| 150–191 | node-positive `df_rotterdam`; hand-written propensity code; truncated weights at 188 | `df_rotterdam_np`; `ps_weights()` on all node-positive women with the ESS table and the reason (inline); restricted `df_rotterdam`, `ps_weights()` with the restricted model, printed object, n by arm (inline), what `pw` is used for, pointer to the comparison in the main vignette | 1, 2, 5, 8 |
| 196–236 | weighted `df_counting()` with `weight.name` | adds `weight_variance = "estimated", weight_model = pw`; ESS table and one sentence after the figure | 2, 5 |
| 276–295, 362–399 | weighted band and RMST calls with `qtau = 0.025` | `qtau` removed; the two `weight_variance` arguments added | 2, 6 |
| 321–338 | GBSG band call with `qtau = 0.025` | `qtau` removed | 6 |
| 303 | "The grey traces show individual resampled processes (W†) …" | the first figure shows the first 20 standardized resampled survival-difference curves; how the critical value and the band are formed; the supremum is over the whole range of the call; influence-based curves with subject weights | 6, 7 |
| 419 | "The resampling distribution is modified accordingly, with the at-risk process replaced by the weighted version" | with weights not all equal the resampled curve is −Ŝ_j(t) Σ G_i w_i ψ_i(t), from the influence contributions of the robust variance; second paragraph on `ps_weights()` and `weight_variance = "estimated"` | 7 |
| 421–439 (unevaluated) | `weight.name = "ps_weights"` twice | a `ps_weights()` call, column `sw.weights`, the two `weight_variance` arguments in both steps | 2, 7 |
| 448 | "The same resampling framework applies naturally to propensity-score weighted analyses" | same multipliers, curves from per-subject influence contributions, estimation of the weights can be carried into the bands | 7 |

### `vignettes/articles/weightedcox_methodology.Rmd`

| Passage | What it said | What it says now | Item |
|---|---|---|---|
| 97–107 | one score-variance formula (Fleming–Harrington form) for all cases; K "the FH weight kernel" | the same formula, for no subject weights or all equal, with K written out as in `score_calculation()`; for subject weights not all equal, the score-residual variance Σ ψ_i(β)² with ψ_i, S₀, z̄ and dΛ̂ as in `ws_score_psi()`; z-score and standard error from it; equality with `coxph(robust = TRUE, ties = "breslow")` for FH(0,0); pointer to estimated weights | 7 |
| 110–118 | resampled score written with "dG_i N_i(t)"; "yields resampling-based confidence intervals"; "extends naturally to constructing simultaneous confidence intervals" | the events-only perturbation G_i dN_i(t) without subject weights; −Σ G_i ψ_i(β̂) with subject weights not all equal; what `cox_rhogamma()` takes from the draws (resampling SE, de-biased estimate, `hr_ci_star`); the same multipliers serve several schemes | 7 |
| new section | — | "Estimated weights": what `ps_weights()` returns; ψ^ew_i = ψ_i + IF_i′B; what changes and what does not; limits (model, time weights, `draws`), concentrated weights, pointer to the evaluation | 2, 7 |
| 311–360 | node-positive Rotterdam, hand-written propensity code, two `cox_rhogamma(draws = 1000)` fits with `hr_ci_star` | ESS table on all node-positive women with the reason (inline); restricted population with `ps_weights()`; ESS table; FH(0,0) under `"estimated"` with the asymptotic interval, FH(0,1) under `"known"` with the resampling interval; the printed lines name the option and the interval | 1, 2, 5 |
| 383 | c_max "is computed via the same synthetic-martingale resampling using common G_i draws" | "can be computed"; `cox_rhogamma()` does not return c_max; it returns `fit_resamples$Z_bstar` per scheme, on the same multipliers | 7 |

### `vignettes/articles/weightedkm_methodology.Rmd`

| Passage | What it said | What it says now | Item |
|---|---|---|---|
| 165 | weights "computed as stabilized IPTW weights … as implemented in the `weightedsurv_examples` vignette" | `ps_weights()` computes inverse probability (default) or stabilized weights; the recipe is for stabilized weights | 2 (consequence) |
| 171 | "(Optional) Truncate extreme weights …" | kept; added that this describes practice and that the package never truncates or rescales weights | 7 (OPEN-2) |
| 199 | "the additional variance contribution from the estimated propensity score is typically small …" | what `weight_variance` does under `"known"` and `"estimated"`; no estimate differs; pointer to the evaluation | 7 (OPEN-3) |
| 201–209 | weighted Greenwood formula; "no adjustment for estimated weights"; resampling "preferred" | Greenwood without weights or with all weights equal; the robust variance Ŝ² Σ w_i² ψ_i(t)² with ψ_i as in the help pages; equality with `survfit(robust = TRUE)`; estimated weights: what `ps_weights()` returns, IF_i, dw_i, D_a, Φ_i, the variance; the resampling variance and that it is not available with `"estimated"`; concentrated weights | 7 |
| 216–220 | standard error "from the (independent) variance estimates" | sum of the arm variances without weights or with weights as known; with `"estimated"` not the sum, with the formula | 7 |
| 228–242 | resampled process written as the log-rank W† with K(t)/Ȳ | resampled deviation of each curve: events-only perturbation without weights; −Ŝ_j Σ G_i w_i ψ_i with weights not all equal; Φ under `"estimated"`; the resampled difference | 7 |
| 244–255 | difference "within the quantile-trimmed range (controlled by the `qtau` parameter)"; sup of the unstandardized difference; band Δ̂ ± c | standard error computed in step 1; the grid and the whole range as called, `qtau` only with `modify_tau = TRUE`; sup of the standardized process; band Δ̂ ± c · se | 6 |
| 257–259 | weighted increments "and the variance proxy accounts for the propensity-score weights" | influence-based curves and robust standard error; `"estimated"`; the arguments to pass | 7 |
| 316 list | — | `ps_weights` added to the function list | — |
| 326 | "including cross-validation of weighted SEs against the `adjustedCurves` and `survfit` packages" | the weighted K-M standard errors agree with `survfit(weights = )`; the curves agree with `adjustedCurves` `iptw_km`, whose standard errors differ | 7 |

### `vignettes/articles/weightedrmst_methodology.Rmd`

| Passage | What it said | What it says now | Item |
|---|---|---|---|
| 112 | the delta-method formula "is the standard variance estimator for RMST implemented in the `survRM2` package and in `weightedsurv`" | the formula is the standard analytic estimator (as in `survRM2`); `weightedsurv` does not use it: the variance is the empirical variance of the integrated resampled curves | 7 |
| 178 | interval determined "from the observed quantiles of the event times (via the `qtau` parameter)" | the supremum is over the whole time grid of the fit; no quantile range trimmed; the critical value is the 0.95 quantile over a second set of `draws_sb` curves | 6 |
| 217–225 | weighted W† with "the kernel function appropriate for the integrated RMST functional" | the resampled RMST difference is the integral of the resampled survival difference; events-only without weights; influence-based with weights not all equal; `"estimated"` must be given to the fit and to `cumulative_rmst_bands()` | 7 |
| 261–269 list | — | `weight_variance`, `weight_model` added to the argument list | 2 |
| 281–296 (unevaluated) | `draws.band = 1000, qtau = 0.025` | `qtau` removed; a `ps_weights()` call and the two `weight_variance` arguments in both steps | 2, 6 |

### `README.md`

| Passage | What it said | What it says now | Item |
|---|---|---|---|
| 231 | "Validation against the `survival` and `adjustedCurves` packages" | "Checks against the `survival` package (weighted Kaplan-Meier standard errors agree with `survfit()`) and a comparison with the `adjustedCurves` package (the curves agree; its `iptw_km` standard errors differ)" | 7 |

## 4. Inventory discrepancies (step 5.3)

Every listed passage was read at `ac76f82`. Lines that are not where the inventory says, or do not say what it says:

| Inventory entry | At `ac76f82` |
|---|---|
| Item 6: martingale 303, 419, 448 say the band covers a trimmed / quantile range | None of the three mentions `qtau`, a quantile or a trimmed range. 303 is the "grey traces … (W†)" sentence, 419 the weighted-resampling sentence, 448 a summary bullet. They are inaccurate in the way item 7 lists them, and were corrected there. |
| Item 7, KM article: 259, 276, 326 are the validation against `adjustedCurves` and `survfit` | Only 326 is. 259 is the paragraph "Extension to Weighted Setting" (corrected with the resampling text). 276 is the sentence on RMST bands "by the same martingale resampling approach, applied to the integrated survival difference process"; it agrees with the code and is unchanged. |
| Item 7, examples 703–713: "`check.seKM` text" | The only prose there is two headings and the note on `draws = 5000` (709). The chunks are at 702–705 and 711–714. |
| Item 7, examples 719–768: `adjustedCurves` call uses `status` | Confirmed (726). Not in the inventory: the call's model uses `grade`, the weightedsurv model `grade3` (same fit in these data); the section starts at 716. |
| Item 7, examples 773–776, 793: comments naming `survdiff` as a source | Confirmed. The same comments also name `wlr_estimates` and `z_score_calculations`, which are not functions of the package. |
| Item 7, martingale: weighted calls 196, 241 | 196 and 241 are chunk headers. The weighted `df_counting()` call is at 207–216; the chunk at 241 only plots that object. |
| Item 4: "flat for time points ≥ 80" | Confirmed at 500, with `abline(v = 80)` at 510. From the data: the last treated event is at 77.9 months and the treated follow-up ends at 87.4; z is 2.108 at 80 months and 2.110 from 87.4 months on. |
| Item 6: band calls with `qtau = 0.025` | All ten listed locations confirmed; no other occurrence in `vignettes/` or `README.md`. |
| Item 6: text on a quantile range | One location not listed: RMST article 178 (changed, see section 3). |
| Section 2 (API facts) | All confirmed in `R/`. In addition: `cox_rhogamma(draws > 0)` is available when the `df_counting()` result was made with `"estimated"`; the `draws` restriction concerns the pointwise variance of `df_counting()`, `KM_diff()` and `plotKM.band_subgroups()`. |
| Section 4: control ESS ratio, all node-positive, 0.963 | The value is 0.96350; at three decimals it prints as 0.964. |

## 5. Found inaccurate, not covered by an item: listed, not changed

`vignettes/weightedsurv_examples.Rmd`
- "### cumulative_rmst_bands … Uses resampling (bootstrap or similar)": the function uses multiplier (martingale) resampling; there is no bootstrap.
- Chunk header of the weighted Cox table: `messages = FALSE` is not a knitr option (`message`); the attach messages of `data.table` are printed in the output.
- Overview, "propensity-score adjustment appears to suggest a stronger impact consistent with the randomized trial": written for the earlier population. In the restricted population the hazard ratio is 0.66 weighted and 0.72 unweighted (GBSG 0.69). The sentence was left as it was.

`vignettes/articles/martingale_resampling_applications.Rmd`
- Line 272: the band is written as P(sup |(Ŝ₁ − Ŝ₀)(t) − (S₁ − S₀)(t)| ≤ c_α), the unstandardized supremum; the code takes the supremum of the standardized process and the band is the estimate ± c · se (the same point item 6 corrects in the KM article at 252–255).
- Lines 69 and 76–85 (and Cox article 371): K(t) carries the factor √{(n₀ + n₁)/(n₀ n₁)}; the score in `score_calculation()` has no such factor (the z-score is the same).
- Lines 263 and 408: "broadly consistent with the randomized GBSG trial, supporting the validity of the observational comparison" and "consistent with the GBSG trial findings" are interpretive sentences written for the all node-positive population; not re-examined, left as they were.

`vignettes/articles/weightedcox_methodology.Rmd`
- Summary bullet "Simultaneous inference: Support for simultaneous confidence intervals …": no function returns a simultaneous interval or c_max (the passage at 383 now says what is returned).
- Schemper weights are written F̄(t⁻)/Ḡ(t); the code uses the left-continuous censoring curve, F̄(t⁻)/Ḡ(t⁻). Xu–O'Quigley weights are written with Ȳ(t)/n; the code divides by Ȳ(t).

`vignettes/articles/weightedkm_methodology.Rmd`
- Summary: "with careful attention to numerical validation against established packages": the only comparisons in the vignettes are with `survfit()` (agreement) and `adjustedCurves` (standard errors differ).

`vignettes/articles/weightedrmst_methodology.Rmd`
- Lines 164–176 give pointwise intervals and bands for a general α; `cumulative_rmst_bands()` has no `alpha` argument and uses 1.96 and the 0.95 quantile.
- Lines 162–178 describe intervals and bands for a single RMST curve; the function computes them for the difference of two arms only.
- Line 277, workflow step 1 "Prepare data via `df_counting` or `get_dfcounting`": `plotKM.band_subgroups()` and `cumulative_rmst_bands()` take the data frame, not a `df_counting()` result.

`README.md`
- "Dependencies": lists Imports survival, stats, graphics and Suggests gt, ggplot2 (with an unclosed parenthesis); `DESCRIPTION` has Imports ggplot2, graphics, grDevices, rlang, stats, survival, utils and nine Suggests.
- Weighting-scheme table: Magirr–Burman "Modest downweighting after cutoff time" (the weight 1/max(S(t), S(t*)) rises until t* and is constant afterwards); Schemper "Adjust for informative censoring" (the Cox article describes it as giving an average hazard ratio approximately free of the censoring distribution).
- The README does not mention `ps_weights()`, `weight_variance` or `df_counting()$ess`.

`vignettes/articles/weightedcox_wald-vs-logrank_simulations.Rmd`
- Not read for accuracy (not in the task's file list). Its `run-simulations` chunk has `eval = TRUE` and the comment "This chunk is NOT evaluated during vignette build".
