# REPORT: Stage 3 release prep — weightedsurv 0.2.0, branch `km-fixes-0.2.0`

Written 2026-11-20. Task: `dev/tasks/TASK_stage3_release_prep.md`.

**Status: stopped at gate G7 (FAIL).** The task was stopped after step 2 of the order of work (text
rider). Steps 3 to 6 (NEWS and `cran-comments.md`, `R CMD check`, tarball, URLs, forestsearch, site
build) were not run. The approved compute (full-size simulations article) was not used. Nothing was
pushed, merged, tagged or submitted. The forestsearch repository was not touched.

## 1. Gate table

| # | Result | Evidence |
|---|---|---|
| G1 | PASS at `94f0085` (steps 1-2 only) | `git diff --name-only e95728a..HEAD`: `README.md`, `vignettes/weightedsurv_examples.Rmd`, `vignettes/articles/martingale_resampling_applications.Rmd`, and files under `dev/` (those of `1a5cf48`, the task document, `dev/audit/stage3/`). `git diff e95728a..HEAD -- R/ man/ tests/ NAMESPACE DESCRIPTION` is empty. |
| G2 | PASS | `dev/audit/stage3/g2_changed_lines.txt`: main vignette, 2 prose lines removed and 2 added, 1 chunk header removed and 1 added (the R5 header); martingale article, 4 prose lines removed and 6 added; no changed line inside a chunk. `dev/audit/stage3/g2_compare.txt`: the main vignette knitted from the file of `e95728a` and from the file after the rider, each in a fresh session against the package built from HEAD; 87 objects in each knit environment, same names. `identical()` on the objects as saved: 57 of 87; the other 30 hold environments (closures, formulas, model objects) and are not `identical()` between two knits of the same file either (control: a second knit of the `e95728a` file gives the same 57 of 87 and the same 30 names). With environments removed (closures and language objects replaced by their deparse), `identical()`: 87 of 87, before against after. |
| G3 | NOT RUN | Stopped at G7. |
| G4 | NOT RUN | Stopped at G7. State found, not changed: `DESCRIPTION` has `Version: 0.2.0` and no `Date` field; the first heading of `NEWS.md` is `# weightedsurv 0.2.0`; `cran-comments.md` does not exist. |
| G5 | NOT RUN | Stopped at G7. |
| G6 | NOT RUN | Stopped at G7. |
| G7 | **FAIL** | `dev/audit/stage3/g7_greps.txt`. `grep -rn "supporting the validity\|bootstrap or similar\|messages = FALSE" vignettes/` returns one line: `vignettes/.Rhistory:2:knitr::opts_chunk$set(echo = TRUE, warning = FALSE, messages = FALSE)`. `vignettes/.Rhistory` is a tracked file (commit `a24f97d`, 2025-09-25; 38 lines; the same line is in it at `e95728a`). It is not among the allowed paths of section 0, so it was not changed. No `.Rmd` file under `vignettes/` matches. The other two parts of G7 hold: `README.md` contains `ps_weights(` (3 lines); `dev/audit/stage3/g7_readme_deps.txt` shows the README lists `identical()` to `DESCRIPTION` Imports (7 entries) and Suggests (9 entries), same order, and 95 opening against 95 closing parentheses in `README.md`. |
| G8 | PASS (rider prose) | Numerals in the changed prose: `95%` (confidence level, design constant, R1b and R1c), `FH(0,0)` (scheme name, R1b), `$1-\alpha$` (R2), and the arguments of inline R (`rho = 0`, `gamma = 0`, digits `2` and `1`). Every result value is inline-computed. In `README.md`: the version bound `(>= 3.1.7)` from `DESCRIPTION`, and the numerals of the copied vignette code (`by.risk=24`, `digits = c(0, 0, 1, 3)`). |
| G9 | NOT RUN | Stopped at G7 (not for lack of urlchecker; it is installed). |
| G10 | NOT RUN | Stopped at G7. |

## 2. Rider

Commits `2d60636` (vignettes) and `a9c66c6` (README). Rendered values below are from a render of
the martingale article against the package built from HEAD (1000 draws, as written).

| # | Passage | What it said | What it says now |
|---|---|---|---|
| R1a | `weightedsurv_examples.Rmd`, Overview | "wherein for comparisons of the non-randomized therapy, in contrast to un-adjusted comparisons, propensity-score adjustment [@cole2004adjusted] appears to suggest a stronger impact consistent with the randomized trial." | "for the non-randomized comparison the hazard ratio is estimated with and without propensity-score weighting [@cole2004adjusted] in a restricted Rotterdam population, and shown next to the hazard ratio of the GBSG trial data, each with its confidence interval. The three values are given in the section "Restricted population: postmenopausal women without chemotherapy"." The Overview precedes the chunks that compute the values, so it points to the section; no number is quoted. |
| R1b | `martingale_resampling_applications.Rmd`, after chunk `compare-km` | "The propensity-score weighted Rotterdam analysis shows a treatment effect pattern broadly consistent with the randomized GBSG trial, supporting the validity of the observational comparison after confounding adjustment." | Three inline-computed hazard ratios (FH(0,0), `cox_rhogamma()$hr_ci_asy`, the expression of the main vignette's `hr_text()` written inline) and two inline-computed position words. Rendered: "Hazard ratios (FH(0,0), with asymptotic 95% intervals): in the restricted Rotterdam population the propensity-score weighted estimate is 0.66 (0.56, 0.78) and the unweighted estimate is 0.72 (0.60, 0.86); in the GBSG trial data the estimate is 0.69 (0.54, 0.89). The weighted Rotterdam estimate is below the unweighted one and below the GBSG estimate. The two datasets differ in population and follow-up, and weighting adjusts only for the covariates of the propensity model." |
| R1c | same article, text under chunk `rmst-combined` | "The weighted RMST curve shows a cumulative treatment benefit that grows steadily over follow-up, consistent with the GBSG trial findings." | Inline values from `get_bands_rotterdam$rmst_maxtau_ci` and `get_bands_gbsg$rmst_maxtau_ci`. Rendered: "The weighted cumulative RMST difference at the last time point of its grid is 15.6 months (pointwise 95% limits 9.7, 21.4) at 118.0 months; in the GBSG trial data (top row) it is 6.8 months (pointwise 95% limits 2.2, 11.4) at 77.9 months. The two values are taken at different time points and are not on a common horizon." |
| R2 | same article, "Rotterdam Weighted Survival Difference Bands" | Coverage statement with the unstandardised supremum, `sup |(S1hat - S0hat)(t) - (S1 - S0)(t)| <= c_alpha`. | Standardised supremum: `sup |Dhat(t) - D(t)| / se[Dhat(t)] <= c_{1-alpha}`, with `c_{1-alpha}` the (1 - alpha) quantile of the resampled supremum of the standardised process, and the band `Dhat(t) ± c_{1-alpha} · se[Dhat(t)]`, "as in `KM_diff()`: wider where the standard error is larger". Code read for this: `R/km_wlr_calculations_helpers.R` lines 681-695 (`Zdhat_star <- dhat_star / sqrt(sig2_dhat)`, `quantile(sups, 1 - alpha)`, `dhat ± c_alpha_band * sqrt(sig2_dhat)`). |
| R3 | `README.md`, "Dependencies" | Imports: survival, stats, graphics. Suggests: "gt (for formatted tables" (unclosed parenthesis), "ggplot2 (for weight scheme visualization)". | Imports: ggplot2, graphics, grDevices, rlang, stats, survival, utils. Suggests: adjustedCurves, data.table, DiagrammeR, dplyr, gt, knitr, pammtools, rmarkdown, testthat (>= 3.1.7). |
| R4 | `README.md`, "Key Features" and "Usage" | No mention of `ps_weights()`, `weight_variance` or `ess`. | One feature bullet, and a subsection "Propensity-Score Weights" before "Baseline Characteristics Table" with one snippet. The snippet's lines are copied from two adjacent chunks of the main vignette that run: `ps-weights-restricted` (the `ps_weights()` call and the weight column) and the unnamed chunk that follows (the weighted `get_dfcounting()` call with `weight_variance = "estimated", weight_model = pw`, and the `$ess` table); no single chunk holds all three. Link: `https://larry-leon.github.io/weightedsurv/articles/articles/estimated_weights_methodology.html`. |
| R5 | `weightedsurv_examples.Rmd` | "Uses resampling (bootstrap or similar) to estimate confidence bands for the RMST curves."; chunk header `{r, messages = FALSE, warning = FALSE}` | "Uses multiplier (martingale) resampling to estimate confidence bands for the RMST curves."; `{r, message = FALSE, warning = FALSE}`. |

## 3. `NEWS.md`

Not changed (step 3 not reached). A `# weightedsurv 0.2.0` section exists (220 lines in the file),
with the groups "Changes in behaviour" and "Bug fixes", written during Stages 2 to 2e. It was not
checked against the git history and the reports in this task.

## 4. `R CMD check`

Not run. Environment found: R 4.5.2, macOS (Darwin 27.2.0); pandoc is available only through
quarto (`/Applications/quarto/bin/tools`, used here as `RSTUDIO_PANDOC` for the knits of G2).

## 5. forestsearch

Not run. The forestsearch repository was not read or changed.

## 6. Site build

Not run. The simulations article was not knitted.

## 7. State of the branch

- HEAD before this report's commit: `94f0085`. Commits since `e95728a`: `1a5cf48` (present at the
  start), `52053a9` (task document), `2d60636` (vignette rider), `a9c66c6` (README rider),
  `94f0085` (G2 and G7 scripts and outputs), and the commit of this report.
- `git status`: clean after the commit of this report.
- Ahead of `origin/km-fixes-0.2.0`: 5 commits, this report's first commit included (`52053a9`, `2d60636`,
  `a9c66c6`, `94f0085`, `e19df63`), plus the commit that corrects this line; `1a5cf48` is on the remote. Not pushed.
- No tag created (existing tag: `v0.1.0`). No merge.

## 8. Found, not covered by this task; listed, not changed

1. `vignettes/.Rhistory` is tracked in git (the cause of the G7 failure). It is not in the tarball
   built from the working tree in this session (`R CMD build` leaves it out).
2. `docs/` is in `.gitignore` and no file under `docs/` is tracked; the task (section 5, step 4)
   says to commit `docs/`. `.gitignore` is not among the allowed paths. The local `docs/` is from
   February 2026.
3. The deployed site (`origin/gh-pages`, "Built site for weightedsurv@0.1.0", 2026-02-12) has no
   page for `estimated_weights_methodology`, so the link added to `README.md` in R4 does not
   resolve until a site that contains the article is deployed. On that branch the articles under
   `vignettes/articles/` are at `articles/articles/<name>.html`; copies of five of them are also at
   `articles/<name>.html`.
4. `cran-comments.md` does not exist in the repository.
5. `DESCRIPTION` already has `Version: 0.2.0` (since Stage 2) and has no `Date` field.
6. The weightedsurv installed in the user library of this machine is 0.1.0; the knits of this task
   used a build of HEAD installed in a temporary library.
7. Not produced because the task stopped: `cran-comments.md`, the tarball from a clean clone, the
   forestsearch comparison, the site.

---

# Resume

Written 2026-11-20. Task: `dev/tasks/TASK_stage3_release_prep_resume.md` (committed first, `8d29c98`),
which amends `dev/tasks/TASK_stage3_release_prep.md`. Steps 1 and 2 were not redone; work resumed at
step 3. The part above is the record of the stop and is unchanged.

**Status: all gates G1 to G10 PASS** (G9 and G6 with the exceptions recorded in the amendments
below). Nothing was pushed, merged, tagged, submitted or deployed by this task. The forestsearch
repository was read only (a `git clone` of it into a temporary directory; no edit, no commit). No
full simulation was run.

## R0. Amendments received during the resume

Three amendments arrived as messages while the resume was running. They are recorded here as given.

1. **G9 amendment (approved by Larry).** The README link to
   `articles/articles/estimated_weights_methodology.html` is expected to return 404 until the site
   is deployed, which is outside this task. It is a known exception; G9 passes if it is the only
   broken URL. If `R CMD check` reports the same URL in a NOTE, it is listed in `cran-comments.md` as
   resolving once the site is deployed. Do not stop: continue with step 5 and step 6, and confirm in
   G6 that the built site has the page at that path.
2. **G5 amendment (approved by Larry).** The forestsearch check uses the current CRAN release
   (0.5.0), from the CRAN source tarball (`download.packages("forestsearch", type = "source")`),
   checked twice in separate temporary libraries: (a) weightedsurv 0.1.0 from CRAN, (b) the
   weightedsurv 0.2.0 tarball of G3. G5 is judged on the CRAN release. The local-clone comparison
   (0.4.1.9003 at `02f88da`) stays in the report as a second table, not as the gate. Done before the
   site build; recorded in `cran-comments.md`.
3. **Site-build amendment (approved by Larry).** The simulations article is not knitted at full size
   inside the site build. The article is restructured (tables and figures from a stored full-run file
   `vignettes/articles/sim_results/wald_vs_logrank_full.rds`; a minimal run of 10 trials per scenario
   on 2 workers evaluated as a labelled illustration; the full-run code shown with `eval = FALSE`,
   with the line that saves the file; the "NOT evaluated" comment made true). The stored file is
   built from the objects of the full-size knit of the vignette revision pass if they are complete,
   otherwise from one run outside the knit. G6 changes: the article renders in under 3 minutes with
   `rmarkdown::render` and inside pkgdown; tables and figures come from the stored file; the stored
   run reports 30,000 completed, 0 failed; the figures agree with those of the full-size knit. G1:
   the article and `vignettes/articles/sim_results/` are allowed paths; the stored file is not in the
   tarball.

The third amendment arrived after `pkgdown::build_site()` had been started with the article as it
was (start 10:00:55, at `d39b8c4`). That build was stopped at 10:03:24, about two minutes after the
article's 12 simulation workers had started; the build and the 12 workers were killed. Its log is
`dev/audit/stage3/g6_build_site_log_stopped_run.txt`. The reference pages and three articles it had
written stayed in `docs/`. The full-size simulation did not complete and was not rerun.

## R1. Gate table (final state)

| # | Result | Evidence |
|---|---|---|
| G1 | PASS | `dev/audit/stage3/g1_paths.txt`. `git diff --name-only e95728a..HEAD` outside `dev/`: `.Rbuildignore`, `NEWS.md`, `README.md`, `cran-comments.md`, `vignettes/.Rhistory` (deletion), `vignettes/articles/martingale_resampling_applications.Rmd`, `vignettes/weightedsurv_examples.Rmd` (the three rider files, committed before the resume), `vignettes/articles/weightedcox_wald-vs-logrank_simulations.Rmd` and `vignettes/articles/sim_results/wald_vs_logrank_full.rds` (site-build amendment). `git diff e95728a..HEAD -- R/ man/ tests/ NAMESPACE DESCRIPTION`: empty. `git ls-files docs`: 0 files. `.gitignore` is unchanged (it already had the line `.Rhistory`). `.Rbuildignore`: one line added, `^cran-comments\.md$`. |
| G2 | PASS | As reported above; not rerun. |
| G3 | PASS | `R CMD check --as-cran` on the tarball built by `R CMD build` from a clean `git clone` at `b8a14ab`: `Status: 2 NOTEs`, 0 errors, 0 warnings; tests `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`. Both NOTEs are listed in `cran-comments.md`. Log: `dev/audit/stage3/g3_00check.log`, `g3_testthat.Rout`. The same check on the tarball built from a clean clone at `d39b8c4`, earlier in the resume: the same status, the same two NOTEs, the same test counts. |
| G4 | PASS | `DESCRIPTION` `Version: 0.2.0`; first heading of `NEWS.md` `# weightedsurv 0.2.0`; `cran-comments.md` exists and names 0.2.0. `dev/audit/stage3/g4_api_diff.R` (exports from `NAMESPACE`, argument lists parsed from `R/`, both against the tag `v0.1.0`) → `g4_api_diff_after.txt`: 46 exports at `v0.1.0`, 47 at HEAD, none removed; added export `ps_weights`; added S3 method `print.ps_weights`; added arguments `weight_variance` and `weight_model` (`df_counting()`, `KM_diff()`, `cumulative_rmst_bands()`, `plotKM.band_subgroups()`) and `mean_rule` (`cumulative_rmst_bands()`); changed defaults `KM_diff(show_resamples)` `TRUE` to `FALSE`, `plotKM.band_subgroups(draws)` 20 to 0 and `(draws.band)` 20 to 1000. Not named in `NEWS.md` at the start of the resume (`g4_api_diff_before.txt`): `print.ps_weights` and the argument `type` of `ps_weights()`. After the added entry: none. |
| G5 | PASS | Judged on the CRAN release (amendment 2). forestsearch 0.5.0, CRAN source tarball: `R CMD check` (a) with weightedsurv 0.1.0 and (b) with weightedsurv 0.2.0: `Status: OK` in both; tests `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]` in both. No error, warning or note in either. Section R5. |
| G6 | PASS | `pkgdown::check_pkgdown()`: "No problems found." (`g6_check_pkgdown.txt`). `pkgdown::build_site()` completed locally at `b8a14ab` in 1.4 minutes (`g6_build_site_log.txt`). Every article of the `_pkgdown.yml` index has an HTML file in `docs/articles/` written by this session's builds (`g6_site_files.txt`). At the path of the README link, `docs/articles/articles/estimated_weights_methodology.html`, the built site has a 334-byte redirect page (meta refresh) to `articles/estimated_weights_methodology.html`, where the article itself is (231,671 bytes). Simulations article: `rmarkdown::render` 31 s; inside pkgdown about 37 s; tables and figures from the stored file; the page states "completed: 30,000; failed: 0" for the stored run; the three figures are byte-identical to those of the full-size knit. Section R6. `git status` shows nothing from `docs/` (ignored; no file tracked). |
| G7 | PASS | `dev/audit/stage3/g7_greps_resume.txt`: `grep -rn --include="*.Rmd" "supporting the validity\|bootstrap or similar\|messages = FALSE" vignettes/` returns nothing; `git ls-files vignettes/ \| grep -i "rhistory\|rdata\|ds_store"` returns nothing. README parts as reported above. |
| G8 | PASS | As reported above; not rerun. |
| G9 | PASS with the known exception (amendment 1) | `urlchecker::url_check()` (`g9_urlcheck.txt`, with pandoc on the `PATH`): "Found 1 invalid URL": `README.md:186`, `https://larry-leon.github.io/weightedsurv/articles/articles/estimated_weights_methodology.html`, 404. It is the only one. `g9_http_status.txt`: that URL 404; `articles/articles/weightedkm_methodology.html` and the site root 200. Run at `d39b8c4`, before the article was restructured; the restructured article adds no URL. |
| G10 | PASS | `dev/audit/stage3/g10_tarball.txt`, tarball from the clean clone at `b8a14ab`: 4,329,234 bytes (under 5 MB), 112 files; no `dev/`, `docs/`, `vignettes/articles/`, `_pkgdown.yml`, `cran-comments.md`, and no `sim_results` (the one `.rds` in it is `build/vignette.rds`, the vignette index written by `R CMD build`). `g10_rbuildignore_match.txt`: `cran-comments.md` matches a `.Rbuildignore` pattern. |

## R2. Amendment A: stray files

- `vignettes/.Rhistory` removed with `git rm` (`3456e12`). `.gitignore` already had the line
  `.Rhistory` (and `.DS_Store`), so it was not changed; the file had been added before that line
  took effect for it.
- Other tracked `.Rhistory`, `.RData` or `.DS_Store` files (`git ls-files`): none.

## R3. `NEWS.md` (amendment D)

No existing entry was rewritten or reordered. Two additions (`d39b8c4`, 12 lines):

1. A new entry after the `ess` entry (Ra), at the end of "Changes in behaviour": the arguments of
   `ps_weights()` (`formula`, `data`, `model`, `type = c("ipw", "stabilized")`) and the
   `print()` method of class `"ps_weights"`. Source: the G4 comparison, `R/ps_weights.R` and
   `man/ps_weights.Rd`.
2. A subsection `## Documentation` with one entry, in the words of amendment D.

Source of each entry (the commit that wrote the entry's closing line, `git blame`; the tag at the end
of each entry is the item name used in the stage reports):

| Entries (tag) | Commit | Report |
|---|---|---|
| B4, B1, B3, B5, B10, B6 ("Changes in behaviour"); B7, B2, (B5, B4, B6), V ("Bug fixes") | `b68df5f` (and `8696d62` for two lines of B4) | `REPORT_km_fixes_stage2.md` |
| B8, B9, B11 | `932ba3a` | `REPORT_km_fixes_stage2b.md` |
| R2 | `7ab4171` | `REPORT_km_fixes_stage2c.md` |
| L2, L3, L4, B12 | `8f820bb` | `REPORT_km_fixes_stage2d.md` |
| (M, R1), L1, (W1-W4), Rd, Rb, Rc, Ra | `c195953` (last edit of each; M and R1 were first written in Stages 2b and 2c, L1 in Stage 2d) | `REPORT_km_fixes_stage2e.md` (and 2b, 2c, 2d for M, R1, L1) |
| `ps_weights()` arguments and print method; Documentation | `d39b8c4` | this task |

The entries written before this task were not checked against the reports line by line in the
resume (amendment D limits the work to the two additions).

The file in full, as committed (232 lines):

    # weightedsurv 0.2.0

    ## Changes in behaviour

    * Input checks (`KM_diff()`, `df_counting()`, `plotKM.band_subgroups()`): missing values in the
      time, event, treatment, weight or strata column, non-finite or negative times, negative weights,
      a missing column (including `strata.name` and `weight.name`) and an empty treatment arm are now
      errors that name the column or arm (and the subgroup in `plotKM.band_subgroups()`). Before, some
      of these dropped into cryptic errors and some gave wrong results without a warning (`NA` survival
      for an `NA` weight, curves from negative weights or negative times). Rows are never dropped
      silently. The event column may be 0/1, logical, or 1/2 as `survival::Surv()` reads it, in all
      three functions (before, `KM_diff()` rejected 1/2). Zero weights are allowed: in `df_counting()`
      those rows no longer cause a false "Discrepancy in KM curve fit" warning, and the Cox model and
      log-rank test are fitted on the positive-weight rows (before, `coxph()` stopped on zero weights).
      On degenerate data (an arm without events, no events at all, one subject per arm)
      `df_counting()` now returns a result; quantities that cannot be computed (for example the hazard
      ratio when an arm has no events) are `NA`, with one warning that says why. In a stratified
      analysis a stratum in which one arm is empty contributes nothing (that arm's stratum curve is
      `NA`, with the same single warning) instead of stopping the call.
      `plotKM.band_subgroups()` asks for `Maxtau` when an arm has no events. `cox_rhogamma()` warns
      when an arm has no events (the estimate is not finite). (B4)

    * `KM_diff()` returns the Kaplan-Meier value at every requested time point, whatever `at_points`
      holds. With `modify_tau = TRUE`, with a user grid that omits event times, and through
      `plotKM.band_subgroups(modify_tau = TRUE)`, the curves, differences, variances, resampled
      processes and bands now equal the Kaplan-Meier estimate at those points; before they could be
      far off (errors above 0.4 and negative survival were seen). `at_points` must now be numeric,
      finite and non-decreasing (unsorted points used to give values outside [0, 1]). With the default
      grid, and in `df_counting()`, nothing changes. (B1)

    * Quantiles (`km_quantile()`, `kmq_calculations()`, `df_counting()$quantile_results`) now equal
      `quantile(survfit(...), probs = 1 - qprob)`: when the curve equals the level on an interval the
      midpoint is returned (before, the left end; for example a median of 5.0 is now 5.5), a plateau at
      the level that runs to the end of follow-up is handled as `survfit()` does, and the confidence
      limits are read the same way when the limit curve is not monotone. `qprob` keeps its meaning, the
      time at which survival drops to `qprob`, and is documented. `df_counting()$quantile_check` now
      asks `survfit()` for the same quantity (`probs = 1 - qprob`, at `conf_level`); a real discrepancy
      is a warning even when `verbose = FALSE`. (B3)

    * `KM_diff()`: the simultaneous-band critical value is the `1 - alpha` quantile (it was always the
      0.95 quantile); results at `alpha = 0.05` are unchanged. `wlr_dhat_estimates()` uses the pooled
      survival at `tzero` (it used the next grid point), so `sig2_dhat` and `cov_wlr_dhat` change
      slightly (by 0-8% on the test datasets; the correlation is unchanged); `tzero` before the first
      time point gives a difference of 0, and after the last time point the values there (before, an
      empty result). Magirr-Burman weights with `mb_tstar` at or beyond the last time point now use
      the survival at the last time point (before, all weights were 1, or infinite where survival was
      0). (B5)

    * The Kaplan-Meier plots (`plot_weighted_km()`, `KM_plot_2sample_weighted_counting()`,
      `plot_km_curves_counting()`) now draw each curve from time zero at S = 1, as a step function
      identical to `survfit()`'s, and the confidence band as a step band that stays at or below 1 and
      closes onto the curve where survival reaches 0. Before, each curve started at its first event
      (the first drop and the line under early censoring marks were missing) and the band joined the
      event times with straight lines and could rise above 1. Everything else on the plot is
      unchanged. (B10)

    * Functions that draw random numbers (`KM_diff()` with `draws` or `draws.band`, `df_counting()`
      and `calculate_risk_event_counts()` with `draws`, `cox_rhogamma()` with `draws`) no longer reset
      the random-number stream of the calling session; results for a given `seedstart` are unchanged.
      (B6)

    * Weighted standard errors: when the weights are not all equal, the standard error of a weighted
      Kaplan-Meier curve is now the robust (infinitesimal-jackknife) standard error,
      `S(t)^2 * sum_i w_i^2 psi_i(t)^2` with `psi_i` subject i's influence contribution, in
      `df_counting()`, `KM_diff()` and `calculate_risk_event_counts()`. It equals the standard error
      of `survival::survfit(weights = w, robust = TRUE)`. **Every weighted standard error changes**,
      and with it everything computed from one: pointwise limits, simultaneous-band limits and
      critical values, the standardised resampled curves, and the confidence limits of weighted
      quantiles. **Unweighted results do not change**, and neither do results with all weights equal
      (both use Greenwood's formula; all weights equal to a constant give the unweighted standard
      error, where in 0.1.0 a constant weight c gave the standard error of the data with every row
      repeated c times). Survival estimates, differences, risk tables and the RMST point estimate do
      not change (the log-rank test and Cox estimate with subject weights change separately, see L1-L4
      below). The standard error does not depend on the scale of the weights
      (weights `w` and `10 * w` give the same answer), so weights are not frequency weights: a weight
      of 2 is not two copies of the row. By default weights are treated as fixed and known; for
      weights estimated by a logistic propensity model see `weight_variance = "estimated"` below. With
      extreme weights (low effective sample size) the standard error can be too small. (M, R1)

    * Weighted resampling: when the weights are not all equal, resampled Kaplan-Meier curves are now
      built from the same per-subject influence contributions, `-S(t) * sum_i G_i w_i psi_i(t)` for
      the same standard-normal multipliers `G_i` as before (before, the multipliers perturbed the
      events only). This applies to `resampling_survival()`, the resampled curves of `KM_diff()`
      (`surv0_star`, `surv1_star`, `dhat_star`, `Zdhat_star`), and the `draws > 0` variance of
      `KM_diff()`, `df_counting()` and `calculate_risk_event_counts()`. **Weighted simultaneous bands
      and weighted RMST bands change** (`cumulative_rmst_bands()`: the RMST variance, pointwise
      limits, critical value and band limits; not the RMST estimate), and the resampling variance now
      agrees with the robust standard error. **Unweighted results, and results with all weights
      equal, do not change**: they keep the events-only perturbation, with the same seeds and draws.
      (R2)

    * Log-rank test with subject weights: when the subject weights are not all equal, the variance of
      the weighted log-rank score is now the sum of squared per-subject score contributions,
      `sum_i psi_i(0)^2` (the formula is in `?df_counting`), in `df_counting()` (unstratified and per
      stratum) and in `wlr_cumulative()` when it is given a `df_counting()` result. For
      Fleming-Harrington (0, 0) weights the z-score squared is the robust score test of
      `survival::coxph(weights = w, robust = TRUE, ties = "breslow")`. **z-scores and p-values of
      subject-weighted log-rank tests change** (`z.score`, `sig2_lr`, `get_score`, `zlogrank_text`,
      `z.score_stratified`, `sig2_lr_stratified`; the cumulative `sig2.score` and `z.score` of
      `wlr_cumulative()`); the score itself (`lr`, `score`, `lr_stratified`) does not. They no longer
      depend on how the weights are scaled: `w`, `10 * w` and `w / sum(w)` give the same z-score
      (before, multiplying the weights by 10 multiplied z by about 3.2, and normalised weights gave
      `NaN`). **Unweighted results, and results with all weights equal, do not change.**
      `wlr_cumulative()` on a bare table of counts (no subject-level data) is unchanged too. By default
      weights are treated as fixed and known; for weights estimated by a logistic propensity model see
      `weight_variance = "estimated"` below. With extreme weights (low effective sample size) the test
      can reject too often. (L1)

    * `cox_rhogamma()` with subject weights: when the weights are not all equal, the standard error of
      the log hazard ratio is now `sqrt(sum_i psi_i(bhat)^2) / i_bhat`, from the same per-subject score
      contributions; for Fleming-Harrington (0, 0) weights it is the robust standard error of
      `survival::coxph(weights = w, robust = TRUE, ties = "breslow")`. **Standard errors, confidence
      intervals, z-scores and p-values of subject-weighted fits change** (`fit$sig_bhat_asy`,
      `hr_ci_asy`, `cox_text_asy`, `fit$sig2_score`, `fit$z.score`, `zlogrank_text`, and `sig2_score`
      and `sig2_beta_asy` of `cox_score_rhogamma(score_only = FALSE)`), and no longer depend on how
      the weights are scaled. The estimate (`bhat`), the scores (`u.beta`, `u.zero`), the information
      and the time weights do not change. **Unweighted fits, and fits with all weights equal, do not
      change.** (L2)

    * `cox_rhogamma(draws > 0)` with subject weights: when the weights are not all equal, the
      resampled scores of `cox_rhogamma_resample()` are now `-sum_i G_i psi_i(bhat)` and
      `-sum_i G_i psi_i(0)`, built from the same per-subject contributions with the same multipliers
      `G_i` as before (before, the multipliers perturbed the events only). **The resampling standard
      error, the de-biased estimate and the resampling interval of subject-weighted fits change**
      (`fit_resamples`, `fit$sig_bhat_star`, `fit$bhat_debiased`, `fit$wald_debiased1`,
      `fit$wald_debiased2`, `hr_ci_star`, `sig2_score_star`, `z.score_debiased`, `cox_text_star`); the
      resampling standard error now agrees with `sig_bhat_asy`. **Unweighted fits, and fits with all
      weights equal, do not change**, with the same seeds and draws. (L3)

    * `wlr_dhat_estimates()` with subject weights: on a `df_counting()` result whose weights are not
      all equal, the variance of the log-rank statistic (as L1), the variance of the survival
      difference at `tzero` and their covariance are now computed from the same per-subject
      contributions; `sig2_dhat` equals `KM_diff(draws = 0)$sig2_dhat` at that time point.
      **`sig2.score`, `z.score`, `sig2_dhat`, `cov_wlr_dhat` and `cor_wlr_dhat` of subject-weighted
      calls change**; `score` and `dhat` do not. **Unweighted results, and results with all weights
      equal, do not change.** (L4)

    * `plotKM.band_subgroups()`: the default of `draws` is now 0 (it was 20). With weights that are not
      all equal and `draws` not supplied, the pointwise standard error, the pointwise limits and the
      standardisation of the simultaneous band now use the analytic (robust) variance instead of a
      variance from 20 resampling draws, so **weighted band plots that do not pass `draws` change**.
      Plots without weights, or with all weights equal, are numerically unaffected. Pass `draws = 20`
      for the earlier behaviour. (B12)

    * `cumulative_rmst_bands()` has a new last argument `mean_rule`. The default, `"trapezoid"`, gives
      the same results as before. `"rectangle"` integrates the survival difference as a step function;
      it then equals the difference of `survfit()`'s restricted means whenever the time grid holds
      every event time (as the default `KM_diff()` grid does; not with `modify_tau = TRUE`). The result
      gains an element `mean_rule`. (B8)

    * `KM_diff()` no longer plots the resampled difference curves by default (`show_resamples = FALSE`);
      use `show_resamples = TRUE` to see them. Results are unchanged. (B9)

    * The two confidence bands of the Kaplan-Meier plots (`plot_weighted_km()`,
      `KM_plot_2sample_weighted_counting()`, `plot_km_curves_counting()`) are now semi-transparent, in
      the same colours, so each stays visible where they overlap. Nothing else on the plot changes. (B11)

    * New option for estimated weights. `ps_weights()` (new export) fits a logistic propensity model, or
      takes a fitted one, and returns inverse-probability or stabilised weights together with the
      influence of the model's coefficients and the derivative of each weight. Passed as `weight_model`
      with `weight_variance = "estimated"` to `df_counting()`, `KM_diff()`, `cumulative_rmst_bands()` or
      `plotKM.band_subgroups()` (two new last arguments), it carries the estimation of the weights into
      the Kaplan-Meier standard errors, the variance of the survival difference (no longer the sum of
      the arm variances), the simultaneous band, the RMST variance and band, the log-rank z-score
      (stratified too), `logrank_results`, `cox_results$score`, `$se`, `$hr_ci` and `$pval`;
      `cox_rhogamma()`, `wlr_cumulative()` and `wlr_dhat_estimates()` take the option from the
      `df_counting()` result they are given. **The default is `weight_variance = "known"`: a call that
      does not ask for `"estimated"` returns what it returned before.** No estimate changes under either
      option. `"estimated"` needs a logistic `glm` fitted on exactly the rows analysed (subgroup analyses
      that reuse full-sample weights must use `"known"`) and `draws = 0`; see the section "Estimated
      weights" of `?df_counting`. Results gain an element `weight_variance`. (W1-W4)

    * Zero weights: whether the weights are "all equal" is now decided on the positive weights. **Calls
      whose weight column contains a zero change**: a 0/1 weight column now gives the variances,
      z-scores and intervals of the positive-weight rows analysed without weights (before, the zeros
      made the weights count as unequal and the robust formulas ran; for example a log-rank z of
      -0.7569 against -0.7554 for the same rows without weights). With zeros among unequal positive
      weights the results are those of the positive-weight rows with the same weights, as before.
      Resampled quantities keep a multiplier per row and are not equal to the subset's draw by draw.
      Calls without a zero weight do not change. (Rd)

    * `df_counting()` with subject weights: `logrank_results` and `cox_results$score` are now robust
      tests of the weighted analysis. **Calls whose positive weights are not all equal change in these
      two outputs only**: `logrank_results` is the robust score test at 0 of the weighted Cox model
      (before, an unweighted `survdiff()` test that ignored the weights), and `cox_results$score` is
      `summary(cox_fit)$robscore` (before, the model-based score test, which reads the weights as
      counts). With subject weights `logrank_results` is the `rho = 0` test whatever the scheme.
      `check_results()` now agrees on weighted Fleming-Harrington (0, 0) results instead of warning.
      New elements `logrank_results$method` (`"survdiff"` or `"coxph_robust_score"`) and
      `cox_results$se` (the standard error of the log hazard ratio behind `hr_ci` and `pval`). Without
      weights, and with all positive weights equal, both outputs are unchanged. (Rb)

    * `plotKM.band_subgroups()`: the default of `draws.band` is now 1000 (it was 20), so the critical
      value of the simultaneous band no longer rests on 20 draws. **Calls that do not pass `draws.band`
      change in the band (critical value, band limits, resampled curves)**; pointwise limits and
      estimates do not. Calls that pass `draws.band` are unchanged; pass `draws.band = 20` for the
      earlier behaviour. (Rc)

    * `df_counting()` returns a new element `ess`: by arm, the number of rows with positive weight, the
      effective sample size `(sum w)^2 / sum(w^2)` and their ratio. A low ratio means a few subjects
      carry most of the weight; standard errors, intervals and tests can then be too optimistic. No
      existing element changes and no warning depends on it. (Ra)

    * `ps_weights()` has the arguments `formula` and `data` (the logistic propensity model is fitted
      by `glm()`), `model` (a fitted logistic `glm`, in place of `formula` and `data`) and
      `type = c("ipw", "stabilized")`. The result has class `"ps_weights"` with a `print()` method
      (`print.ps_weights()`) that shows, by arm, the number of subjects, the smallest and largest
      weight, the effective sample size and its ratio to the number of subjects.

    ## Bug fixes

    * Stratified `df_counting()` no longer depends on `survival` being attached or recent enough to
      find `strata()`. (B7)

    * Weighted Kaplan-Meier curves that end on an event now end at exactly 0. Before, a residue of
      about 1e-17, sometimes negative, made `plot_weighted_km()` stop with "values outside [0,1]" and
      gave "NaNs produced" warnings. Confidence limits are `NA` where survival is 0 (as in `survfit()`),
      and the plotted band no longer receives `NaN`. (B2)

    * `wlr_dhat_estimates()` and `wlr_cumulative()` with `scheme = "schemper"` now check for the
      censoring curve and say what is missing; `cumulative_rmst_bands()` says that it needs a
      `KM_diff()` fit made with `draws.band > 0`; `KM_diff(draws.band < 20, show_resamples = TRUE)` no
      longer fails; `get_validated_weights()` no longer prints to the console. (B5, B4, B6)

    * Internal tidying with no change in results: a local function inside `df_counting()` no longer
      shadows the exported `check_km_curve()`, and unused code was removed. (V)

    ## Documentation

    * The vignettes and articles were revised for 0.2.0: the Rotterdam example now uses postmenopausal
      women without chemotherapy and `ps_weights()` with `weight_variance = "estimated"`; effective
      sample size tables were added; a new article gives the statistical details of estimated weights.

## R4. `R CMD check`

- Environment: R 4.5.2 (2025-10-31), aarch64-apple-darwin20, macOS 27.2.1; pandoc from quarto
  (`/Applications/quarto/bin/tools` on the `PATH` and as `RSTUDIO_PANDOC`); `DISPLAY` unset.
- Tarball: `weightedsurv_0.2.0.tar.gz`, built by `R CMD build` from a clean `git clone` at
  `b8a14ab`; 4,329,234 bytes; sha256 `4835afcb6e9a0bc68186263cb26f4a2f3895d08fa4e79eae8263679e8c30d13b`.
  It has the same 112 file names as the tarball built from the clean clone at `d39b8c4`; two files
  differ in content between the two builds, `DESCRIPTION` (written by `R CMD build`) and
  `inst/doc/weightedsurv_examples.html` (rebuilt vignette).
- `R CMD check --as-cran`: `Status: 2 NOTEs`. 0 errors, 0 warnings.
  1. "checking CRAN incoming feasibility ... NOTE": maintainer line, and "Found the following
     (possibly) invalid URLs": the README link of amendment 1, status 404.
  2. "checking HTML version of manual ... NOTE": "Skipping checking HTML validation: 'tidy' doesn't
     look like recent enough HTML Tidy."
- Tests: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`. Examples and the vignette rebuild: OK.
- `cran-comments.md` (new, `f010a16`): the environment, the result with both NOTEs, the reverse
  dependency result on the CRAN release of forestsearch, and a summary of the changes since 0.1.0,
  including that weighted standard errors and tests change and why. It states that no other
  environment (win-builder, R-hub) was used. It is in `.Rbuildignore`.
- Packages used in vignette chunks against `DESCRIPTION` (recorded, not changed;
  `dev/audit/stage3/vignette_packages.txt`): the main vignette and five of the six articles use only
  packages in Imports or Suggests (the main vignette names devtools in one comment line, not in
  code). The simulations article uses 11 packages that are in neither: doFuture, foreach, future,
  future.callr, gsDesign2, htmltools, parallelly, simtrial, tibble, tictoc, tidyr. The articles are
  not in the tarball (`^vignettes/articles$` in `.Rbuildignore`).

## R5. forestsearch

Reverse dependencies on CRAN, `tools::package_dependencies("weightedsurv", reverse = TRUE, which =
"all")`: `forestsearch` (one package). CRAN has weightedsurv 0.1.0 and forestsearch 0.5.0;
forestsearch lists weightedsurv in Imports.

Both comparisons: plain `R CMD check` on one tarball, twice, with `R_LIBS` set to a temporary
library holding (a) weightedsurv 0.1.0 installed from the CRAN source tarball, (b) weightedsurv
0.2.0 installed from the tarball of the first G3 run (clean clone at `d39b8c4`; `R/`, `NAMESPACE`,
`man/` and `tests/` identical to `b8a14ab`). `dev/audit/stage3/forestsearch/library_versions.txt` shows the version each library
resolves to. Logs in `dev/audit/stage3/forestsearch/`.

**Gate table: CRAN release, forestsearch 0.5.0** (source tarball from CRAN, sha256
`0a4082df35be3903b013316c5b4cebfe251cb8acfaceaf3a382b40f269351da5`)

| | (a) weightedsurv 0.1.0 | (b) weightedsurv 0.2.0 |
|---|---|---|
| Status | OK | OK |
| Errors / warnings / notes | 0 / 0 / 0 | 0 / 0 / 0 |
| Tests | `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]` | `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]` |
| Failing tests | none | none |
| Examples | OK | OK |
| Vignettes (check and rebuild) | OK | OK |
| Difference in the check output, timings aside | | one INFO line: "installed size is 5.5Mb" against 5.4Mb in (a) |

**Second table: local clone, forestsearch 0.4.1.9003** (`~/Documents/GitHub/forestsearch`, branch
`master`, commit `02f88da`, working tree clean before and after; tarball built by `R CMD build` from
a `git clone` of it in a temporary directory, with the 0.2.0 library on the path for the vignette
build). The local clone has a lower version number than the CRAN release.

| | (a) weightedsurv 0.1.0 | (b) weightedsurv 0.2.0 |
|---|---|---|
| Status | 1 ERROR | 1 ERROR |
| Errors / warnings / notes | 1 / 0 / 0 | 1 / 0 / 0 |
| Tests | `[ FAIL 1 | WARN 18 | SKIP 92 | PASS 4206 ]` | `[ FAIL 1 | WARN 18 | SKIP 92 | PASS 4206 ]` |
| Failing test | `test-mrct-region-sims-record.R:53:3`, "record = 'point' is the call without the argument and the pre-change call" | the same test |
| Examples | OK | OK |
| Vignettes (check and rebuild) | OK | OK |
| Difference in the check output, timings aside | | none |

The failing test of the local clone compares two results that differ in the last digits
(for example `0.55000000000000004` against `0.55000000000000016`); it fails identically with both
weightedsurv versions.

**weightedsurv functions that forestsearch calls** (`weightedsurv_calls_in_forestsearch.txt`; the
same three functions and the same 13 call sites in `R/` of the CRAN release and of the local clone;
`NAMESPACE` imports `df_counting` and `plot_weighted_km`; no call in `tests/`; the local clone's
articles `paper_simulations.qmd` and `forestsearch-2024.Rmd` call `df_counting()` and
`plot_weighted_km()` after `library(weightedsurv)`). No call passes `weight.name`,
`weight_variance` or `weight_model`.

| Function | Calls in `R/` | Changed in Stages 2c to 2e |
|---|---|---|
| `df_counting()` | 7, without subject weights | Results: the changes of 2c (robust Kaplan-Meier variance), 2d (score-residual variance) and 2e (`logrank_results`, `cox_results$score`, zero weights) apply when the positive weights are not all equal or a weight is zero; `NEWS.md` states that unweighted results do not change. Defaults: none changed. Arguments: `weight_variance`, `weight_model` added last (2e). New result elements: `ess`, `weight_variance`, `logrank_results$method`, `cox_results$se` (2e). |
| `plot_weighted_km()` | 5, on unweighted `df_counting()` results | No change in 2c to 2e. (Stages 2 and 2b changed the drawing: curves from time zero, step bands, semi-transparent bands.) |
| `plotKM.band_subgroups()` | 1, with `draws.band` supplied and `draws` not supplied, without weights | Defaults: `draws` 20 to 0 (2d, B12) and `draws.band` 20 to 1000 (2e, Rc). forestsearch passes `draws.band`, so Rc does not apply to its call; it does not pass `draws`, and `NEWS.md` (B12) states that plots without weights are numerically unaffected. Arguments: `weight_variance`, `weight_model` added last (2e). |

## R6. Site build

- `pkgdown::check_pkgdown()`: "No problems found." (pkgdown 2.2.1), before the first build and again
  at `b8a14ab`.
- First `build_site()` (at `d39b8c4`, article as it was): stopped after 2 minutes 29 seconds on the
  site-build amendment (section R0). Not completed.
- Second `build_site(preview = FALSE)` at `b8a14ab`: completed, exit status 0, 1.4 minutes
  wall-clock (10:06:47 to 10:08:11). Not deployed. `docs/` is not committed (ignored).
- `docs/` was not emptied before the builds: of its 501 files, 157 are older than this session
  (pkgdown assets and files of the February build that these builds did not rewrite).

Per article (all seven of the `_pkgdown.yml` index rendered; pkgdown does not print render times;
the seconds are the interval between the write time of the file and of the one before it in the
completed build):

| Index entry | File in the built site | Rendered | Seconds |
|---|---|---|---|
| `weightedsurv_examples` | `articles/weightedsurv_examples.html` | yes | 12 |
| `articles/weightedcox_methodology` | `articles/weightedcox_methodology.html` | yes | 5 |
| `articles/weightedkm_methodology` | `articles/weightedkm_methodology.html` | yes | 1 |
| `articles/weightedrmst_methodology` | `articles/weightedrmst_methodology.html` | yes | 1 |
| `articles/estimated_weights_methodology` | `articles/estimated_weights_methodology.html` | yes | 2 |
| `articles/weightedcox_wald-vs-logrank_simulations` | `articles/weightedcox_wald-vs-logrank_simulations.html` | yes | 37 |
| `articles/martingale_resampling_applications` | `articles/martingale_resampling_applications.html` | yes | 4 |

For three of them (estimated weights, martingale, Cox methodology) the log of the completed build
shows "Reading" without a "Writing" line; their files nevertheless carry write times within the
completed build (10:07:02, 10:07:06, 10:07:11). The reason for the missing log line was not
investigated.

**Path of the README link.** pkgdown 2.2.1 writes each article of `vignettes/articles/` to
`articles/<name>.html` and adds a redirect page at `articles/articles/<name>.html` for each of the
six. So the built site has, at `articles/articles/estimated_weights_methodology.html`, a redirect
(meta refresh, canonical link) to
`https://larry-leon.github.io/weightedsurv/articles/estimated_weights_methodology.html`, not the
article itself. On the deployed site the README link returns 404 (the site built for 0.1.0 is the
one deployed); it resolves only after a site that contains the article is deployed.

**Simulations article (site-build amendment).** Commit `b8a14ab`.

- *Stored file.* The objects of the full-size knit of the vignette revision pass were found in
  `dev/scratch/vignette_before/sims/env.rds` (94,922,577 bytes, written 2026-11-20 00:49:37, sha256
  `58f8f9ad…f510f2`; not committed, ignored by `dev/scratch/.gitignore`) and are complete: `res_out`
  with 30,000 rows (5,000 per scenario, 6 scenarios), 80 columns, no missing value; the knitted
  output of that run (`sims.md`) has "Completed: 30000 | Failed: 0" and "Elapsed: 71.8 minutes". The
  stored file was built from them by `dev/audit/stage3/build_sim_results.R` (output
  `build_sim_results.txt`). **No simulation was run.**
- *Content and size.* Per-trial estimates. With all 80 columns the per-trial results are
  13,180,664 bytes compressed (xz), over 5 MB. The stored file holds the 27 per-trial columns that
  the tables and figures use (the scenario, the 15 z-statistics that are not de-biased, and the 11
  estimates of the two tables), all 30,000 rows, the setup and timing elements of `res_out`, and a
  `meta` list: **4,699,236 bytes (xz), under 5 MB** (5,708,351 bytes with gzip). The 53 other
  columns (de-biased z, coverage indicators, `maxcombop`, estimates of schemes the tables do not
  show) are not in the file. The kept columns read back `identical()` to the source columns.
- *Settings stated in the text* (inline R from `meta`): 6 scenarios, named; 5,000 trials per
  scenario; seed 8316951; 200 martingale draws; 12 workers; weightedsurv 0.2.0, source commit
  `ac76f82`; date 2026-11-20; completed 30,000, failed 0. The seed, trial count, counts and
  elapsed time come from the stored objects. The martingale draws and the worker count are not in
  the stored objects; they are the values of the article's `run-simulations` chunk, which was
  byte-identical to `ac76f82` when knitted. The commit is that of the working tree knitted (`R/`
  is identical at `ac76f82` and at HEAD). The date is that of the file; the knit started on
  2026-11-19 at about 23:37 and ended on 2026-11-20 at 00:49.
- *Article.* The section "Simulation Execution" now has three parts: "Minimal run (illustration
  only)", chunk `run-simulations-demo`, evaluated, 10 trials per scenario, 2 workers, result
  `res_demo`, stated as not used for any table or figure; "Full run (code shown, not evaluated)",
  chunk `run-simulations-full`, `eval = FALSE`, the full call as before, then the lines that select
  the columns and `saveRDS(stored, "sim_results/wald_vs_logrank_full.rds", compress = "xz")`, with
  the comment "This chunk is NOT evaluated when the article is built"; "Stored full run", chunk
  `load-results`, which reads the file into `res_out`, followed by the settings list. One sentence
  was added at the top of "Results". No chunk after `load-results` was changed. Nothing before
  "Simulation Execution" was changed.
- *Render time.* `rmarkdown::render`: 31 s (`g6_sims_render.txt`). `knitr::knit` with
  `dev/scratch/knit_vignettes.R`: 32 s, 0 warnings escaping the chunks (`g6_sims_knit.txt`).
  Inside pkgdown: about 37 s. The minimal run printed "Completed: 60 | Failed: 0", 0.44 minutes.
- *Figures and tables against the full-size knit, and how compared.* The restructured article was
  knitted with the same script as the full-size knit (`dev/scratch/knit_vignettes.R`: `knitr::knit`,
  png at 96 dpi) against the weightedsurv 0.2.0 build of G3, and the output compared with
  `dev/scratch/vignette_before/sims/` by `dev/audit/stage3/g6_sims_compare.R`
  (`g6_sims_compare.txt`): the three figure files (`main-3tests-plot`, `all-fh-plot`,
  `simtrial-comparison-plot`) are byte-identical; the data of the three `gt` tables and of the three
  ggplot objects are `identical()`; `n_sim`, `pnull_threshold`, `time_inhours` and the annotation
  tables are `identical()`. The figures of the pkgdown page itself (different device settings) were
  not compared pixel by pixel; the page has the same three figures and the text "completed: 30,000;
  failed: 0".
- *Tarball.* The stored file is not in the tarball (G10).

## R7. State of the branch

- HEAD before the commit of this section: `b8a14ab`. Commits of the resume: `8d29c98` (task
  document), `3456e12` (`vignettes/.Rhistory` removed), `d39b8c4` (`NEWS.md`, `.Rbuildignore`, G4),
  `f010a16` (`cran-comments.md`, outputs of G3, G4, G5, G9, G10), `b8a14ab` (simulations article and
  stored file), and the commit of this section. Commits since `e95728a`: those, and the six listed
  in the first part (`1a5cf48`, `52053a9`, `2d60636`, `a9c66c6`, `94f0085`, `e19df63`) with `6cbcded`.
- `git status`: clean after the commit of this section.
- Remote: `origin/km-fixes-0.2.0` is at `d39b8c4`. Ahead of it: `f010a16`, `b8a14ab` and the commit
  of this section (3 commits). No command of this task pushed; see R8, item 1.
- Tags: `v0.1.0` only. No merge. No submission. No deployment.

## R8. Found, not covered by this task; listed, not changed

1. The reflog of `refs/remotes/origin/km-fixes-0.2.0` in this clone has two entries "update by
   push" during the resume: 09:42:36 (to `6cbcded`) and 09:50:28 (to `d39b8c4`). No command run by
   this task was a push. At the end of the first part of this report the branch was 6 commits ahead
   of the remote.
2. The first `R CMD check --as-cran` of the resume did not finish: it stayed at "checking use of S3
   registration" for 14 minutes, in a subprocess waiting to open the X display named by `DISPLAY`
   (XQuartz socket) while loading tcltk. It was stopped and the check rerun with `DISPLAY` unset
   (63 s). The partial log is `dev/audit/stage3/g3_check_stdout_first_attempt_stalled.txt`; with
   pandoc not on the `PATH` it also had the NOTE "Files 'README.md' or 'NEWS.md' cannot be checked
   without 'pandoc' being installed".
3. The path `articles/articles/<name>.html` of the README link is a redirect page in a site built
   by pkgdown 2.2.1; the articles themselves are at `articles/<name>.html` (section R6).
4. The local forestsearch clone (0.4.1.9003 at `02f88da`) has one failing test under both
   weightedsurv versions (section R5), and a version number below the CRAN release (0.5.0).
5. forestsearch 0.5.0: the check output with weightedsurv 0.2.0 reports an installed size of 5.5Mb
   against 5.4Mb with 0.1.0 (an INFO line, not a NOTE).
6. The simulations article needs 11 packages that are in neither Imports nor Suggests to build
   (section R4); its minimal run needs them at site build. The article's note that the helper
   functions "will be included in a future CRAN release" was not touched.
7. `docs/` holds files from the February build next to those of this session (section R6).
8. `vignettes/.DS_Store` exists on disk; it is not tracked.
9. The weightedsurv installed in the user library of this machine is 0.1.0 (as in the first part).
10. `cran-comments.md` records one test environment; no other was used.
