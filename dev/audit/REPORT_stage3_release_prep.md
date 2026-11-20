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
- Ahead of `origin/km-fixes-0.2.0`: 5 commits before this report's commit, 6 after. Not pushed.
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
